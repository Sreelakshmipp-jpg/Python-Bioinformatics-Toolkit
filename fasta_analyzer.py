#!/usr/bin/env python3
"""
fasta_analyzer.py - FASTA analyzer with sequence validation (standard library only).

Improvements over a naive analyzer
  * Full IUPAC nucleotide alphabet (R, Y, S, W, K, M, B, D, H, V, N) is recognised
    as nucleotide, not mistaken for protein.
  * RNA is analysed as RNA (RNA residue masses, U counted as U, no DNA-only Tm).
  * Invalid characters are detected, reported and excluded (never silently ignored).
  * Molecular weight is reported as a min-max range when ambiguity codes are present.
  * Tm uses nearest-neighbor thermodynamics (SantaLucia 1998) with salt correction
    for short unambiguous DNA; a long-sequence formula is used above 60 nt.
  * ORF search supports several genetic codes and optional alternative start codons,
    and handles ambiguous codons (translated to X when the amino acid is not certain).

Usage
  python fasta_analyzer.py seqs.fasta
  python fasta_analyzer.py seqs.fasta --table 11 --alt-starts --min-orf 90
  python fasta_analyzer.py seqs.fasta --type rna --csv results.csv --strict
"""

import argparse
import csv
import math
import sys
from collections import Counter

# ---------------------------------------------------------------------------
# Alphabets
# ---------------------------------------------------------------------------

IUPAC_EXPAND = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "U",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
IUPAC_DNA = set("ACGTRYSWKMBDHVN")
IUPAC_RNA = set("ACGURYSWKMBDHVN")
PROTEIN_ALPHABET = set("ACDEFGHIKLMNPQRSTVWY") | set("BZXUOJ")

DNA_COMP = str.maketrans("ACGTRYKMBVDHSWN", "TGCAYRMKVBHDSWN")
RNA_COMP = str.maketrans("ACGURYKMBVDHSWN", "UGCAYRMKVBHDSWN")

# Fraction of possible bases that are G or C, per IUPAC symbol (for GC range)
GC_POSSIBLE = {"G", "C", "S", "R", "Y", "K", "M", "B", "D", "H", "V", "N"}  # can be G/C
GC_CERTAIN = {"G", "C", "S"}                                                # always G/C

# ---------------------------------------------------------------------------
# Genetic codes
# ---------------------------------------------------------------------------

GENETIC_CODES = {
    1: {"name": "Standard", "changes": {},
        "starts": {"ATG", "TTG", "CTG"}},
    2: {"name": "Vertebrate mitochondrial",
        "changes": {"AGA": "*", "AGG": "*", "ATA": "M", "TGA": "W"},
        "starts": {"ATT", "ATC", "ATA", "ATG", "GTG"}},
    11: {"name": "Bacterial, archaeal and plant plastid", "changes": {},
         "starts": {"TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"}},
}


def build_codon_table(table_id):
    bases = "TCAG"
    amino = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    table, i = {}, 0
    for b1 in bases:
        for b2 in bases:
            for b3 in bases:
                table[b1 + b2 + b3] = amino[i]
                i += 1
    table.update(GENETIC_CODES[table_id]["changes"])
    return table


def translate_codon(codon, table):
    """Translate a codon; ambiguous codons give an amino acid only if all
    possible expansions agree, otherwise X."""
    aa = table.get(codon)
    if aa is not None:
        return aa
    try:
        options = {
            table[a + b + c]
            for a in IUPAC_EXPAND[codon[0]]
            for b in IUPAC_EXPAND[codon[1]]
            for c in IUPAC_EXPAND[codon[2]]
        }
    except KeyError:
        return "X"
    return options.pop() if len(options) == 1 else "X"


def translate(seq, table):
    return "".join(translate_codon(seq[i:i + 3], table) for i in range(0, len(seq) - 2, 3))


# ---------------------------------------------------------------------------
# Molecular weights (linear, single-stranded, 5'-OH / 3'-OH)
# Residue = nucleoside monophosphate minus water; chain correction = -61.96
# ---------------------------------------------------------------------------

DNA_RESIDUE = {"A": 313.2076, "C": 289.1846, "G": 329.2070, "T": 304.1961}
RNA_RESIDUE = {"A": 329.2052, "C": 305.1800, "G": 345.2046, "U": 306.1673}
END_CORRECTION = 61.96

# ---------------------------------------------------------------------------
# Nearest-neighbor parameters, SantaLucia (1998) unified, DNA/DNA
# key = 5'->3' dinucleotide; value = (dH kcal/mol, dS cal/mol/K)
# ---------------------------------------------------------------------------

NN = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}
INIT_GC = (0.1, -2.8)
INIT_AT = (2.3, 4.1)
SYMMETRY_DS = -1.4
GAS_CONSTANT = 1.987  # cal/mol/K


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_fasta(handle):
    """Yield (header, sequence) tuples. Handles multi-line sequences."""
    header, chunks = None, []
    for line in handle:
        line = line.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(chunks)
            header, chunks = line[1:].strip(), []
        else:
            if header is None:
                raise ValueError("File does not look like FASTA (sequence before first '>' header).")
            chunks.append("".join(line.split()).upper())
    if header is not None:
        yield header, "".join(chunks)


# ---------------------------------------------------------------------------
# Sequence helpers
# ---------------------------------------------------------------------------

def reverse_complement(seq, alphabet="DNA"):
    table = RNA_COMP if alphabet == "RNA" else DNA_COMP
    return seq.translate(table)[::-1]


# Letters that occur in proteins but are never nucleotide symbols
PROTEIN_ONLY = PROTEIN_ALPHABET - IUPAC_DNA - {"U"}
MIN_LEN_CONFIDENT = 20  # below this, an all-IUPAC-letter sequence is genuinely ambiguous


def classify(seq, forced):
    """Return (type, note). type is 'DNA', 'RNA', 'Protein' or 'Mixed';
    note is a warning string when auto-detection is only a best guess."""
    if forced != "auto":
        return {"dna": "DNA", "rna": "RNA", "protein": "Protein"}[forced], None

    chars = set(seq)
    # U never occurs in DNA and is very rare in proteins, so a sequence containing U,
    # no T, and only RNA IUPAC symbols is RNA even when it is short.
    if "U" in chars and "T" not in chars and chars <= IUPAC_RNA:
        return "RNA", None

    letters = [c for c in seq if c.isalpha()]
    if not letters:
        return "DNA", None  # nothing usable; reported later as invalid characters

    core = sum(1 for c in letters if c in "ACGTUN") / len(letters)
    has_protein_only = any(c in PROTEIN_ONLY for c in letters)

    # Nucleotide if mostly A/C/G/T/U/N (real proteins are only ~25% of these), or if the
    # sequence is long enough that having no protein-only letters at all is
    # overwhelmingly unlikely for a protein (heavy-IUPAC-ambiguity DNA/RNA).
    nucleotide = core >= 0.75 or (not has_protein_only and len(letters) >= MIN_LEN_CONFIDENT)

    if nucleotide:
        if "T" in chars and "U" in chars:
            return "Mixed", None
        return ("RNA" if "U" in chars else "DNA"), None

    note = None
    if not has_protein_only:
        note = (f"short sequence made only of IUPAC letters was classified as Protein; "
                f"if it is nucleotide, rerun with --type dna or --type rna")
    return "Protein", note


def gc_percent_strict(counts, unambiguous):
    return 100.0 * (counts["G"] + counts["C"]) / unambiguous if unambiguous else None


def gc_range(seq):
    n = len(seq)
    low = sum(1 for c in seq if c in GC_CERTAIN)
    high = sum(1 for c in seq if c in GC_POSSIBLE)
    return 100.0 * low / n, 100.0 * high / n


def gc_skew(counts):
    g, c = counts["G"], counts["C"]
    return (g - c) / (g + c) if (g + c) else None


def molecular_weight_range(seq, stype):
    """Return (min, max) mass in g/mol. Ambiguity codes widen the range."""
    residues = RNA_RESIDUE if stype == "RNA" else DNA_RESIDUE
    lo = hi = 0.0
    for sym, n in Counter(seq).items():
        expansion = IUPAC_EXPAND[sym]
        if stype == "RNA":
            expansion = expansion.replace("T", "U")  # IUPAC table is written in the DNA alphabet
        weights = [residues[b] for b in expansion]
        lo += n * min(weights)
        hi += n * max(weights)
    return lo - END_CORRECTION, hi - END_CORRECTION


def tm_nearest_neighbor(seq, na_mm, oligo_nm):
    dh = ds = 0.0
    for i in range(len(seq) - 1):
        h, s = NN[seq[i:i + 2]]
        dh += h
        ds += s
    for end in (seq[0], seq[-1]):
        h, s = INIT_GC if end in "GC" else INIT_AT
        dh += h
        ds += s
    self_comp = seq == reverse_complement(seq)
    if self_comp:
        ds += SYMMETRY_DS
    x = 1 if self_comp else 4
    ds += 0.368 * (len(seq) - 1) * math.log(na_mm / 1000.0)  # salt correction
    ct = oligo_nm * 1e-9
    return dh * 1000.0 / (ds + GAS_CONSTANT * math.log(ct / x)) - 273.15


def compute_tm(seq, stype, n_ambig, na_mm, oligo_nm):
    """Return (tm or None, method description)."""
    if stype == "RNA":
        return None, "not computed (RNA duplexes need RNA-specific parameters)"
    if n_ambig:
        return None, "not computed (ambiguous bases present)"
    n = len(seq)
    if n < 8:
        return None, "not computed (shorter than 8 nt)"
    if n <= 60:
        return (tm_nearest_neighbor(seq, na_mm, oligo_nm),
                f"nearest-neighbor, SantaLucia 1998; {na_mm:g} mM Na+, {oligo_nm:g} nM strand")
    gc = 100.0 * (seq.count("G") + seq.count("C")) / n
    tm = 81.5 + 16.6 * math.log10(na_mm / 1000.0) + 0.41 * gc - 675.0 / n
    return tm, f"long-sequence GC formula (approximate); {na_mm:g} mM Na+"


def iter_orfs(seq, table, starts, min_len, both_strands, nested=False):
    """Yield ORFs (start codon ... stop codon; stop included in length).

    Default: one ORF per stop codon per frame, beginning at the first start codon
    (the maximal ORF). With nested=True every start codon in an open region also
    yields its own, shorter ORF ending at the same stop (alternative start sites).
    ORFs without a stop codon are skipped. seq must use the T alphabet."""
    strands = [("+", seq)]
    if both_strands:
        strands.append(("-", reverse_complement(seq, "DNA")))
    length = len(seq)
    for strand, s in strands:
        for frame in range(3):
            i = frame
            while i + 3 <= length:
                if s[i:i + 3] not in starts:
                    i += 3
                    continue
                j = i
                while j + 3 <= length and translate_codon(s[j:j + 3], table) != "*":
                    j += 3
                if j + 3 > length:
                    break  # no stop downstream in this frame, so no later start can end either
                begins = [k for k in range(i, j, 3) if s[k:k + 3] in starts] if nested else [i]
                for k in begins:
                    orf_len = j + 3 - k
                    if orf_len < min_len:
                        continue
                    if strand == "+":
                        start, end = k + 1, j + 3
                    else:
                        start, end = length - (j + 3) + 1, length - k
                    yield {
                        "strand": strand, "frame": frame + 1,
                        "start": start, "end": end, "length": orf_len,
                        "start_codon": s[k:k + 3], "nested": k != i,
                        # initiator codon is read as Met regardless of identity
                        "protein": "M" + translate(s[k + 3:j], table),
                    }
                i = j + 3


def find_longest_orf(seq, table, starts, min_len, both_strands):
    """Longest ORF, or None. Nested starts are always shorter than the maximal ORF
    sharing their stop codon, so they never change this result."""
    return max(iter_orfs(seq, table, starts, min_len, both_strands),
               key=lambda o: o["length"], default=None)


# ---------------------------------------------------------------------------
# Record analysis
# ---------------------------------------------------------------------------

def describe_chars(chars, limit=8):
    shown = sorted(chars)[:limit]
    text = ", ".join(repr(c) for c in shown)
    return text + (", ..." if len(chars) > limit else "")


def analyze_record(header, raw_seq, args):
    r = {"id": header.split()[0] if header else "", "description": header,
         "warnings": [], "invalid": 0}

    seq = raw_seq.upper()
    gaps = seq.count("-") + seq.count(".")
    if gaps:
        seq = seq.replace("-", "").replace(".", "")
        r["warnings"].append(f"{gaps} gap character(s) removed")

    if not seq:
        r.update(type="Empty", length=0)
        r["warnings"].append("empty sequence")
        r["invalid"] = 1
        return r

    stype, note = classify(seq, args.type)
    r["type"] = stype
    if note:
        r["warnings"].append(note)

    if stype == "Mixed":
        r["length"] = len(seq)
        r["warnings"].append("contains both T and U; cannot be DNA or RNA (use --type to force)")
        r["invalid"] = 1
        return r

    # ---- Protein -----------------------------------------------------------
    if stype == "Protein":
        if seq.endswith("*"):
            seq = seq[:-1]
        bad = Counter(c for c in seq if c not in PROTEIN_ALPHABET)
        if bad:
            r["invalid"] = sum(bad.values())
            r["warnings"].append(f"{r['invalid']} invalid protein character(s) removed: {describe_chars(bad)}")
            seq = "".join(c for c in seq if c in PROTEIN_ALPHABET)
        r["length"] = len(seq)
        r["composition"] = dict(Counter(seq).most_common())
        return r

    # ---- Nucleotide ----------------------------------------------------------
    valid = IUPAC_RNA if stype == "RNA" else IUPAC_DNA
    bad = Counter(c for c in seq if c not in valid)
    if bad:
        r["invalid"] = sum(bad.values())
        r["warnings"].append(
            f"{r['invalid']} invalid {stype} character(s) removed: {describe_chars(bad)}")
        seq = "".join(c for c in seq if c in valid)

    r["length"] = len(seq)
    if not seq:
        r["warnings"].append("no valid characters left to analyse")
        return r

    counts = Counter(seq)
    unamb_letters = "ACGU" if stype == "RNA" else "ACGT"
    unambiguous = sum(counts[b] for b in unamb_letters)
    n_ambig = len(seq) - unambiguous
    r["counts"] = {b: counts.get(b, 0) for b in unamb_letters}
    r["ambiguity_counts"] = {k: v for k, v in sorted(counts.items()) if k not in unamb_letters}
    r["n_ambig"] = n_ambig
    r["ambig_pct"] = 100.0 * n_ambig / len(seq)
    if r["ambig_pct"] > 10:
        r["warnings"].append(f"{r['ambig_pct']:.1f}% ambiguous bases; derived values are less reliable")

    r["gc_strict"] = gc_percent_strict(counts, unambiguous)
    r["gc_min"], r["gc_max"] = gc_range(seq)
    r["gc_skew"] = gc_skew(counts)
    r["mw_min"], r["mw_max"] = molecular_weight_range(seq, stype)
    r["tm"], r["tm_method"] = compute_tm(seq, stype, n_ambig, args.na_mm, args.oligo_nm)

    # ORFs: codon logic is identical for RNA, so use the T alphabet internally
    dna_view = seq.replace("U", "T")
    table = build_codon_table(args.table)
    starts = GENETIC_CODES[args.table]["starts"] if args.alt_starts else {"ATG"}
    if args.strands == "auto":
        both = stype == "DNA"
    else:
        both = args.strands == "both"
    r["orf_strands"] = "both strands" if both else "forward strand only"
    r["orf"] = find_longest_orf(dna_view, table, starts, args.min_orf, both)
    if args.list_orfs:
        orfs = list(iter_orfs(dna_view, table, starts, args.min_orf, both, nested=args.nested))
        orfs.sort(key=lambda o: (-o["length"], o["strand"], o["start"]))
        r["orf_list"] = orfs
    r["revcomp"] = reverse_complement(seq, stype)
    return r


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def fmt_range(lo, hi, spec):
    if abs(hi - lo) < 1e-9:
        return format(lo, spec)
    return f"{format(lo, spec)} - {format(hi, spec)}"


def print_report(r, args):
    print("=" * 72)
    print(f"ID          : {r['id']}")
    print(f"Description : {r['description']}")
    print(f"Type        : {r['type']}")
    print(f"Length      : {r['length']}")
    for w in r["warnings"]:
        print(f"WARNING     : {w}")

    if r["type"] == "Protein":
        if r.get("composition"):
            print("Composition : " + ", ".join(f"{k}={v}" for k, v in r["composition"].items()))
        return
    if "counts" not in r:
        return

    counts = "  ".join(f"{k}={v}" for k, v in r["counts"].items())
    print(f"Base counts : {counts}")
    if r["n_ambig"]:
        amb = "  ".join(f"{k}={v}" for k, v in r["ambiguity_counts"].items())
        print(f"Ambiguous   : {r['n_ambig']} ({r['ambig_pct']:.1f}%)  [{amb}]")

    if r["gc_strict"] is not None:
        print(f"GC content  : {r['gc_strict']:.2f}% (unambiguous bases only)")
    if r["n_ambig"]:
        print(f"GC range    : {r['gc_min']:.2f}% - {r['gc_max']:.2f}% (over all positions)")
    if r["gc_skew"] is not None:
        print(f"GC skew     : {r['gc_skew']:+.3f}")

    if r["tm"] is not None:
        print(f"Tm          : {r['tm']:.1f} C  [{r['tm_method']}]")
    else:
        print(f"Tm          : {r['tm_method']}")

    kind = "RNA" if r["type"] == "RNA" else "DNA"
    print(f"Mol. weight : {fmt_range(r['mw_min'], r['mw_max'], '.1f')} g/mol "
          f"(linear single-stranded {kind}, 5'-OH/3'-OH)")

    orf = r["orf"]
    if orf:
        prot = orf["protein"]
        shown = prot if len(prot) <= 60 else prot[:57] + "..."
        print(f"Longest ORF : {orf['length']} nt incl. stop, strand {orf['strand']}, "
              f"frame {orf['frame']}, pos {orf['start']}-{orf['end']}, start {orf['start_codon']}")
        print(f"Protein     : {shown} ({len(prot)} aa)")
    else:
        print(f"Longest ORF : none >= {args.min_orf} nt with start and stop ({r['orf_strands']})")

    if args.list_orfs and "orf_list" in r:
        orfs = r["orf_list"]
        label = "all start codons (nested included)" if args.nested else "maximal ORF per stop codon"
        print(f"ORF list    : {len(orfs)} ORF(s) >= {args.min_orf} nt, {label}; "
              f"showing up to {args.max_orfs} by length")
        for o in orfs[:args.max_orfs]:
            tag = " nested" if o["nested"] else ""
            prot = o["protein"] if len(o["protein"]) <= 30 else o["protein"][:27] + "..."
            print(f"  {o['strand']} f{o['frame']} {o['start']}-{o['end']} "
                  f"{o['length']} nt {o['start_codon']}{tag}  {prot} ({len(o['protein'])} aa)")

    if args.revcomp:
        rc = r["revcomp"]
        print("Rev. comp.  :")
        for i in range(0, len(rc), 60):
            print("  " + rc[i:i + 60])


def print_summary(results, args):
    print("=" * 72)
    print("SUMMARY")
    print(f"Records      : {len(results)}")
    types = Counter(r["type"] for r in results)
    print("Types        : " + ", ".join(f"{k}={v}" for k, v in sorted(types.items())))
    lengths = [r["length"] for r in results]
    print(f"Total length : {sum(lengths)} (min {min(lengths)}, max {max(lengths)}, "
          f"mean {sum(lengths) / len(lengths):.1f})")
    gcs = [r["gc_strict"] for r in results if r.get("gc_strict") is not None]
    if gcs:
        print(f"Mean GC      : {sum(gcs) / len(gcs):.2f}% across {len(gcs)} nucleotide record(s)")
    flagged = sum(1 for r in results if r["warnings"])
    print(f"Warnings     : {flagged} record(s) with warnings")
    print(f"Genetic code : table {args.table} ({GENETIC_CODES[args.table]['name']}); "
          f"start codons: {'table-specific alternatives' if args.alt_starts else 'ATG only'}")
    print("Notes        : Tm assumes a perfectly complementary duplex with no secondary structure.")
    print("               ORFs lacking a stop codon (partial ORFs) are not reported.")


CSV_FIELDS = [
    "id", "type", "length", "invalid_removed", "ambiguous_bases", "ambiguous_pct",
    "A", "C", "G", "T", "U", "gc_percent", "gc_min", "gc_max", "gc_skew",
    "tm", "tm_method", "mw_min", "mw_max",
    "orf_length", "orf_strand", "orf_frame", "orf_start", "orf_end", "orf_start_codon",
    "warnings",
]


def write_csv(results, path):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS, restval="")
        writer.writeheader()
        for r in results:
            row = {"id": r["id"], "type": r["type"], "length": r["length"],
                   "invalid_removed": r["invalid"], "warnings": "; ".join(r["warnings"])}
            if "counts" in r:
                row.update(r["counts"])
                row["ambiguous_bases"] = r["n_ambig"]
                row["ambiguous_pct"] = round(r["ambig_pct"], 2)
                if r["gc_strict"] is not None:
                    row["gc_percent"] = round(r["gc_strict"], 2)
                row["gc_min"] = round(r["gc_min"], 2)
                row["gc_max"] = round(r["gc_max"], 2)
                if r["gc_skew"] is not None:
                    row["gc_skew"] = round(r["gc_skew"], 3)
                if r["tm"] is not None:
                    row["tm"] = round(r["tm"], 1)
                row["tm_method"] = r["tm_method"]
                row["mw_min"] = round(r["mw_min"], 1)
                row["mw_max"] = round(r["mw_max"], 1)
                if r["orf"]:
                    o = r["orf"]
                    row.update(orf_length=o["length"], orf_strand=o["strand"],
                               orf_frame=o["frame"], orf_start=o["start"],
                               orf_end=o["end"], orf_start_codon=o["start_codon"])
            writer.writerow(row)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description="Analyze and validate sequences in a FASTA file.")
    p.add_argument("fasta", help="Path to FASTA file ('-' for stdin)")
    p.add_argument("--type", choices=["auto", "dna", "rna", "protein"], default="auto",
                   help="Force sequence type instead of auto-detection")
    p.add_argument("--table", type=int, choices=sorted(GENETIC_CODES), default=1,
                   help="NCBI genetic code table: 1 standard, 2 vertebrate mito, 11 bacterial")
    p.add_argument("--alt-starts", action="store_true",
                   help="Also accept the table's alternative start codons (default: ATG only)")
    p.add_argument("--strands", choices=["auto", "both", "forward"], default="auto",
                   help="ORF search strands (auto: both for DNA, forward for RNA)")
    p.add_argument("--min-orf", type=int, default=30,
                   help="Minimum ORF length in nt, stop codon included (default: 30)")
    p.add_argument("--list-orfs", action="store_true",
                   help="List all ORFs (not just the longest) for each nucleotide record")
    p.add_argument("--nested", action="store_true",
                   help="With --list-orfs, also list ORFs from internal start codons sharing a stop")
    p.add_argument("--max-orfs", type=int, default=20,
                   help="Maximum ORFs to print per record with --list-orfs (default 20)")
    p.add_argument("--na-mm", type=float, default=50.0, help="Na+ concentration in mM for Tm (default 50)")
    p.add_argument("--oligo-nm", type=float, default=250.0,
                   help="Oligonucleotide strand concentration in nM for Tm (default 250)")
    p.add_argument("--revcomp", action="store_true", help="Print reverse complement of each sequence")
    p.add_argument("--csv", metavar="FILE", help="Write a CSV summary table")
    p.add_argument("--strict", action="store_true",
                   help="Exit with status 2 if any record has invalid or empty sequence data")
    args = p.parse_args()

    if args.na_mm <= 0 or args.oligo_nm <= 0:
        sys.exit("Error: --na-mm and --oligo-nm must be positive.")

    try:
        handle = sys.stdin if args.fasta == "-" else open(args.fasta)
        with handle:
            records = list(parse_fasta(handle))
    except (OSError, ValueError) as exc:
        sys.exit(f"Error: {exc}")
    if not records:
        sys.exit("No FASTA records found.")

    results = [analyze_record(h, s, args) for h, s in records]

    seen = Counter(r["id"] for r in results)
    for r in results:
        if seen[r["id"]] > 1:
            r["warnings"].append("duplicate sequence ID in file")

    for r in results:
        print_report(r, args)
    print_summary(results, args)

    if args.csv:
        write_csv(results, args.csv)
        print(f"\nCSV written to {args.csv}")

    if args.strict and any(r["invalid"] for r in results):
        sys.exit(2)


if __name__ == "__main__":
    main()