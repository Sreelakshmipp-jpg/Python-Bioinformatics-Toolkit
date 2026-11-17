#!/usr/bin/env python3
"""
protein_analyzer.py - Protein sequence analyzer (standard library only).

Per record it reports:
  * Validation (invalid characters, gaps, internal stop codons, nucleotide look-alikes)
  * Amino acid composition, including non-standard/ambiguous symbols (U, O, B, Z, J, X)
  * Molecular weight (average mass); a min-max range when ambiguity codes are present
  * Isoelectric point and net charge at a chosen pH (Bjellqvist pK values)
  * Extinction coefficient at 280 nm (Pace et al. 1995), reduced and oxidised forms
  * GRAVY (Kyte-Doolittle), aliphatic index (Ikai 1980), aromaticity
  * Instability index (Guruprasad et al. 1990)
  * Helix / turn / sheet residue-class fractions (crude composition measure, not a prediction)
  * Optional hydropathy scan for hydrophobic segments (Kyte-Doolittle sliding window)

Usage
  python protein_analyzer.py proteins.fasta
  python protein_analyzer.py proteins.fasta --ph 5.5 --hydropathy --csv results.csv
"""

import argparse
import csv
import sys
from collections import Counter

# ---------------------------------------------------------------------------
# Reference data
# ---------------------------------------------------------------------------

STANDARD = "ACDEFGHIKLMNPQRSTVWY"
NONSTANDARD = "UOBZJX"
VALID = set(STANDARD) | set(NONSTANDARD)

# Average residue masses (g/mol; amino acid minus water), as used by ExPASy.
RESIDUE_MASS = {
    "A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886, "C": 103.1388,
    "E": 129.1155, "Q": 128.1307, "G": 57.0519, "H": 137.1411, "I": 113.1594,
    "L": 113.1594, "K": 128.1741, "M": 131.1926, "F": 147.1766, "P": 97.1167,
    "S": 87.0782, "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326,
    "U": 150.0379,   # selenocysteine
    "O": 237.3018,   # pyrrolysine
}
WATER = 18.01524

# Possible identities of ambiguity codes, for the molecular weight range.
AMBIGUITY = {"B": "DN", "Z": "EQ", "J": "IL", "X": STANDARD}

# Kyte-Doolittle hydropathy
KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# pK values: Bjellqvist et al. (1993, 1994) - the method used by the ExPASy tool.
POSITIVE_PK = {"Nterm": 7.5, "K": 10.0, "R": 12.0, "H": 5.98}
NEGATIVE_PK = {"Cterm": 3.55, "D": 4.05, "E": 4.45, "C": 9.0, "Y": 10.0}
CTERM_PK = {"D": 4.55, "E": 4.75}
NTERM_PK = {"A": 7.59, "M": 7.0, "S": 6.93, "P": 8.36, "T": 6.82, "V": 7.44, "E": 7.7}

# Extinction coefficient contributions at 280 nm in water (M^-1 cm^-1), Pace et al. 1995
EXT_W, EXT_Y, EXT_CYSTINE = 5500, 1490, 125

# Residue classes for the crude secondary-structure composition (as in Biopython)
HELIX_SET, TURN_SET, SHEET_SET = set("VIYFWL"), set("NPGS"), set("EMAL")

# Instability index dipeptide weights (Guruprasad et al. 1990).
# Row = first residue, column = second residue.
DIWV_TEXT = """
          A      C      D      E      F      G      H      I      K      L      M      N      P      Q      R      S      T      V      W      Y
A      1.00  44.94  -7.49   1.00   1.00   1.00  -7.49   1.00   1.00   1.00   1.00   1.00  20.26   1.00   1.00   1.00   1.00   1.00   1.00   1.00
C      1.00   1.00  20.26   1.00   1.00   1.00  33.60   1.00   1.00  20.26  33.60   1.00  20.26  -6.54   1.00   1.00  33.60  -6.54  24.68   1.00
D      1.00   1.00   1.00   1.00  -6.54   1.00   1.00   1.00  -7.49   1.00   1.00   1.00   1.00   1.00  -6.54  20.26 -14.03   1.00   1.00   1.00
E      1.00  44.94  20.26  33.60   1.00   1.00  -6.54  20.26   1.00   1.00   1.00   1.00  20.26  20.26   1.00  20.26   1.00   1.00 -14.03   1.00
F      1.00   1.00  13.34   1.00   1.00   1.00   1.00   1.00 -14.03   1.00   1.00   1.00  20.26   1.00   1.00   1.00   1.00   1.00   1.00  33.60
G     -7.49   1.00   1.00  -6.54   1.00  13.34   1.00  -7.49  -7.49   1.00   1.00  -7.49   1.00   1.00   1.00   1.00  -7.49   1.00  13.34  -7.49
H      1.00   1.00   1.00   1.00  -9.37  -9.37   1.00  44.94  24.68   1.00   1.00  24.68  -1.88   1.00   1.00   1.00  -6.54   1.00  -1.88  44.94
I      1.00   1.00   1.00  44.94   1.00   1.00  13.34   1.00  -7.49  20.26   1.00   1.00  -1.88   1.00   1.00   1.00   1.00  -7.49   1.00   1.00
K      1.00   1.00   1.00   1.00   1.00  -7.49   1.00  -7.49   1.00  -7.49  33.60   1.00  -6.54  24.64  33.60   1.00   1.00  -7.49   1.00   1.00
L      1.00   1.00   1.00   1.00   1.00   1.00   1.00   1.00  -7.49   1.00   1.00   1.00  20.26  33.60  20.26   1.00   1.00   1.00  24.68   1.00
M     13.34   1.00   1.00   1.00   1.00   1.00  58.28   1.00   1.00   1.00  -1.88   1.00  44.94  -6.54  -6.54  44.94  -1.88   1.00   1.00  24.68
N      1.00  -1.88   1.00   1.00 -14.03 -14.03   1.00  44.94  24.68   1.00   1.00   1.00  -1.88  -6.54   1.00   1.00  -7.49   1.00  -9.37   1.00
P     20.26  -6.54  -6.54  18.38  20.26   1.00   1.00   1.00   1.00   1.00  -6.54   1.00  20.26  20.26  -6.54  20.26   1.00  20.26  -1.88   1.00
Q      1.00  -6.54  20.26  20.26  -6.54   1.00   1.00   1.00   1.00   1.00   1.00   1.00  20.26  20.26   1.00  44.94   1.00  -6.54   1.00  -6.54
R      1.00   1.00   1.00   1.00   1.00  -7.49  20.26   1.00   1.00   1.00   1.00  13.34  20.26  20.26  58.28  44.94   1.00   1.00  58.28  -6.54
S      1.00  33.60   1.00  20.26   1.00   1.00   1.00   1.00   1.00   1.00   1.00   1.00  44.94  20.26  20.26  20.26   1.00   1.00   1.00   1.00
T      1.00   1.00   1.00  20.26  13.34  -7.49   1.00   1.00   1.00   1.00   1.00 -14.03   1.00  -6.54   1.00   1.00   1.00   1.00 -14.03   1.00
V      1.00   1.00 -14.03   1.00   1.00  -7.49   1.00   1.00  -1.88   1.00   1.00   1.00  20.26   1.00   1.00   1.00  -7.49   1.00   1.00  -6.54
W    -14.03   1.00   1.00   1.00   1.00  -9.37  24.68   1.00   1.00  13.34  24.68  13.34   1.00   1.00   1.00   1.00 -14.03  -7.49   1.00   1.00
Y     24.68   1.00  24.68  -6.54   1.00  -7.49  13.34   1.00   1.00   1.00  44.94   1.00  13.34   1.00 -15.91   1.00  -7.49   1.00  -9.37  13.34
"""


def load_diwv():
    lines = DIWV_TEXT.strip().splitlines()
    cols = lines[0].split()
    table = {}
    for line in lines[1:]:
        parts = line.split()
        table[parts[0]] = {c: float(v) for c, v in zip(cols, parts[1:])}
    return table


DIWV = load_diwv()


# ---------------------------------------------------------------------------
# Parsing and validation
# ---------------------------------------------------------------------------

def parse_fasta(handle):
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


def looks_like_nucleotide(seq):
    letters = [c for c in seq if c.isalpha()]
    return len(letters) >= 20 and all(c in "ACGTUN" for c in letters)


def describe_chars(chars, limit=8):
    shown = sorted(chars)[:limit]
    return ", ".join(repr(c) for c in shown) + (", ..." if len(chars) > limit else "")


# ---------------------------------------------------------------------------
# Calculations
# ---------------------------------------------------------------------------

def molecular_weight_range(counts):
    """(min, max) average mass in g/mol. Ambiguity codes widen the range."""
    lo = hi = WATER
    for sym, n in counts.items():
        if sym in AMBIGUITY:
            masses = [RESIDUE_MASS[b] for b in AMBIGUITY[sym]]
        else:
            masses = [RESIDUE_MASS[sym]]
        lo += n * min(masses)
        hi += n * max(masses)
    return lo, hi


def partial_charge(pk, ph):
    return 1.0 / (10.0 ** (ph - pk) + 1.0)


def net_charge(counts, first, last, ph):
    """Net charge at pH (Henderson-Hasselbalch, Bjellqvist pK values).
    Residue-specific pK values apply to the N-terminal and C-terminal residues.
    Ambiguous and non-standard residues are ignored."""
    pos = dict(POSITIVE_PK)
    neg = dict(NEGATIVE_PK)
    pos["Nterm"] = NTERM_PK.get(first, POSITIVE_PK["Nterm"])
    neg["Cterm"] = CTERM_PK.get(last, NEGATIVE_PK["Cterm"])
    positive = partial_charge(pos["Nterm"], ph)
    negative = 1.0 / (10.0 ** (neg["Cterm"] - ph) + 1.0)
    for aa in ("K", "R", "H"):
        positive += counts.get(aa, 0) * partial_charge(pos[aa], ph)
    for aa in ("D", "E", "C", "Y"):
        negative += counts.get(aa, 0) / (10.0 ** (neg[aa] - ph) + 1.0)
    return positive - negative


def isoelectric_point(counts, first, last):
    """pH at which the net charge is zero (bisection on 0-14)."""
    lo, hi = 0.0, 14.0
    while hi - lo > 1e-5:
        mid = (lo + hi) / 2.0
        if net_charge(counts, first, last, mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def instability_index(seq):
    """Guruprasad et al. 1990. Dipeptides with non-standard symbols are skipped
    (approximate when such symbols are present)."""
    total = 0.0
    for i in range(len(seq) - 1):
        a, b = seq[i], seq[i + 1]
        if a in DIWV and b in DIWV[a]:
            total += DIWV[a][b]
    return 10.0 * total / len(seq)


def gravy(seq):
    """Mean Kyte-Doolittle hydropathy over residues that have a defined value."""
    vals = [KD[c] for c in seq if c in KD]
    return (sum(vals) / len(vals)) if vals else None, len(seq) - len(vals)


def aliphatic_index(counts, length):
    pct = lambda aa: 100.0 * counts.get(aa, 0) / length
    return pct("A") + 2.9 * pct("V") + 3.9 * (pct("I") + pct("L"))


def hydropathy_regions(seq, window, threshold):
    """Regions where the Kyte-Doolittle window mean is >= threshold.
    Overlapping qualifying windows are merged. Returns list of (start, end, max_mean), 1-based."""
    n = len(seq)
    if n < window:
        return []
    regions = []
    current = None
    for s in range(0, n - window + 1):
        vals = [KD[c] for c in seq[s:s + window] if c in KD]
        mean = (sum(vals) / len(vals)) if vals else None
        if mean is not None and mean >= threshold:
            if current is None:
                current = [s + 1, s + window, mean]
            else:
                current[1] = s + window
                current[2] = max(current[2], mean)
        elif current is not None:
            regions.append(tuple(current))
            current = None
    if current is not None:
        regions.append(tuple(current))
    return regions


# ---------------------------------------------------------------------------
# Record analysis
# ---------------------------------------------------------------------------

def analyze_record(header, raw, args):
    r = {"id": header.split()[0] if header else "", "description": header,
         "warnings": [], "invalid": 0, "skipped": False}
    seq = raw.upper()

    gaps = seq.count("-") + seq.count(".")
    if gaps:
        seq = seq.replace("-", "").replace(".", "")
        r["warnings"].append(f"{gaps} gap character(s) removed")
    if seq.endswith("*"):
        seq = seq[:-1]
        r["warnings"].append("terminal stop symbol '*' removed")

    if not seq:
        r.update(length=0, skipped=True, invalid=1)
        r["warnings"].append("empty sequence")
        return r

    if looks_like_nucleotide(seq) and not args.force:
        r.update(length=len(seq), skipped=True, invalid=1)
        r["warnings"].append("looks like a nucleotide sequence; protein properties not computed "
                             "(use fasta_analyzer.py, or --force to override)")
        return r

    bad = Counter(c for c in seq if c not in VALID)
    if bad:
        r["invalid"] = sum(bad.values())
        r["warnings"].append(
            f"{r['invalid']} invalid character(s) removed: {describe_chars(bad)} "
            f"(removal joins neighbouring residues; dipeptide-based values may be affected)")
        seq = "".join(c for c in seq if c in VALID)
    if not seq:
        r.update(length=0, skipped=True)
        r["warnings"].append("no valid residues left to analyse")
        return r

    n = len(seq)
    counts = Counter(seq)
    r["length"] = n
    r["counts"] = {aa: counts.get(aa, 0) for aa in STANDARD}
    r["nonstandard"] = {k: counts[k] for k in NONSTANDARD if counts.get(k)}
    n_special = sum(r["nonstandard"].values())
    if n_special:
        pct = 100.0 * n_special / n
        r["warnings"].append(
            f"{n_special} non-standard/ambiguous residue(s) ({pct:.1f}%): mass is given as a range; "
            f"charge, pI, GRAVY and instability index ignore them")
    if "C" in counts:
        r["cys"] = counts["C"]

    r["mw_min"], r["mw_max"] = molecular_weight_range(counts)
    first, last = seq[0], seq[-1]
    r["charge"] = net_charge(counts, first, last, args.ph)
    r["pi"] = isoelectric_point(counts, first, last)
    r["gravy"], r["gravy_skipped"] = gravy(seq)
    r["aliphatic"] = aliphatic_index(counts, n)
    r["aromaticity"] = (counts.get("F", 0) + counts.get("W", 0) + counts.get("Y", 0)) / n
    r["instability"] = instability_index(seq)

    w, y, c = counts.get("W", 0), counts.get("Y", 0), counts.get("C", 0)
    reduced = w * EXT_W + y * EXT_Y
    oxidised = reduced + (c // 2) * EXT_CYSTINE
    r["ext_reduced"], r["ext_oxidised"] = reduced, oxidised
    mw_mid = (r["mw_min"] + r["mw_max"]) / 2.0
    r["abs01_reduced"] = reduced / mw_mid
    r["abs01_oxidised"] = oxidised / mw_mid

    r["helix"] = sum(counts.get(a, 0) for a in HELIX_SET) / n
    r["turn"] = sum(counts.get(a, 0) for a in TURN_SET) / n
    r["sheet"] = sum(counts.get(a, 0) for a in SHEET_SET) / n

    if args.hydropathy:
        r["hydro_regions"] = hydropathy_regions(seq, args.window, args.hydro_threshold)
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
    print(f"Length      : {r['length']} aa")
    for w in r["warnings"]:
        print(f"WARNING     : {w}")
    if r["skipped"]:
        return

    print("Composition :")
    items = [f"{aa}={r['counts'][aa]:<3d}({100.0 * r['counts'][aa] / r['length']:4.1f}%)" for aa in STANDARD]
    for i in range(0, len(items), 5):
        print("  " + "  ".join(items[i:i + 5]))
    if r["nonstandard"]:
        print("  Other     : " + ", ".join(f"{k}={v}" for k, v in r["nonstandard"].items()))

    print(f"Mol. weight : {fmt_range(r['mw_min'], r['mw_max'], ',.2f')} Da (average mass)")
    print(f"pI          : {r['pi']:.2f}  (Bjellqvist pK values; Cys treated as free thiols)")
    print(f"Net charge  : {r['charge']:+.2f} at pH {args.ph:g}")
    c = r["counts"]
    print(f"Charged     : negative D+E={c['D'] + c['E']}, positive K+R={c['K'] + c['R']}, His={c['H']}")
    lower = " (lower bound: ambiguous residues could add absorbance)" if r["nonstandard"] else ""
    print(f"Ext. coeff. : {r['ext_reduced']} M-1cm-1 (Cys reduced), "
          f"{r['ext_oxidised']} M-1cm-1 (all Cys as cystines) at 280 nm{lower}")
    print(f"Abs 0.1%    : {r['abs01_reduced']:.3f} (reduced), {r['abs01_oxidised']:.3f} (oxidised) "
          f"for 1 g/L")
    if r["gravy"] is not None:
        extra = f" ({r['gravy_skipped']} residue(s) without a value ignored)" if r["gravy_skipped"] else ""
        print(f"GRAVY       : {r['gravy']:+.3f}{extra}")
    print(f"Aliphatic   : {r['aliphatic']:.2f}")
    print(f"Aromaticity : {r['aromaticity']:.3f}  (F+W+Y fraction)")
    label = "predicted stable (< 40)" if r["instability"] < 40 else "predicted unstable (>= 40)"
    print(f"Instability : {r['instability']:.2f}  {label}")
    print(f"Residue-class fractions (composition only, not a structure prediction): "
          f"helix {r['helix']:.2f}, turn {r['turn']:.2f}, sheet {r['sheet']:.2f}")

    if args.hydropathy:
        regs = r.get("hydro_regions", [])
        if r["length"] < args.window:
            print(f"Hydropathy  : sequence shorter than the {args.window}-residue window")
        elif not regs:
            print(f"Hydropathy  : no {args.window}-residue window with mean >= {args.hydro_threshold:g}")
        else:
            print(f"Hydropathy  : {len(regs)} hydrophobic segment(s) "
                  f"(window {args.window}, mean >= {args.hydro_threshold:g})")
            for s, e, m in regs[:10]:
                print(f"  {s}-{e} ({e - s + 1} aa), max window mean {m:.2f}")
            if len(regs) > 10:
                print(f"  ... and {len(regs) - 10} more")


def print_summary(results):
    ok = [r for r in results if not r["skipped"]]
    print("=" * 72)
    print("SUMMARY")
    print(f"Records      : {len(results)} ({len(ok)} analysed, {len(results) - len(ok)} skipped)")
    if ok:
        lengths = [r["length"] for r in ok]
        print(f"Length       : min {min(lengths)}, max {max(lengths)}, mean {sum(lengths) / len(lengths):.1f} aa")
        mws = [(r["mw_min"] + r["mw_max"]) / 2 for r in ok]
        print(f"Mol. weight  : min {min(mws):,.1f}, max {max(mws):,.1f} Da")
    print(f"Warnings     : {sum(1 for r in results if r['warnings'])} record(s) with warnings")
    print("Notes        : pI and charge assume standard pK values and ignore modifications and")
    print("               disulfides; other tools use different pK sets and can differ by several tenths.")
    print("               The instability index is an empirical in vitro estimate, not a measurement.")


CSV_FIELDS = (["id", "length", "invalid_removed", "mw_min", "mw_max", "pi", "charge", "ph",
               "gravy", "aliphatic_index", "aromaticity", "instability_index",
               "ext_reduced", "ext_oxidised", "helix_fraction", "turn_fraction", "sheet_fraction"]
              + list(STANDARD) + list(NONSTANDARD) + ["warnings"])


def write_csv(results, path, ph):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS, restval="")
        writer.writeheader()
        for r in results:
            row = {"id": r["id"], "length": r["length"], "invalid_removed": r["invalid"],
                   "warnings": "; ".join(r["warnings"])}
            if not r["skipped"]:
                row.update(mw_min=round(r["mw_min"], 2), mw_max=round(r["mw_max"], 2),
                           pi=round(r["pi"], 2), charge=round(r["charge"], 2), ph=ph,
                           aliphatic_index=round(r["aliphatic"], 2),
                           aromaticity=round(r["aromaticity"], 3),
                           instability_index=round(r["instability"], 2),
                           ext_reduced=r["ext_reduced"], ext_oxidised=r["ext_oxidised"],
                           helix_fraction=round(r["helix"], 3), turn_fraction=round(r["turn"], 3),
                           sheet_fraction=round(r["sheet"], 3))
                if r["gravy"] is not None:
                    row["gravy"] = round(r["gravy"], 3)
                row.update(r["counts"])
                for k in NONSTANDARD:
                    row[k] = r["nonstandard"].get(k, 0)
            writer.writerow(row)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description="Analyze protein sequences in a FASTA file.")
    p.add_argument("fasta", help="Path to protein FASTA file ('-' for stdin)")
    p.add_argument("--ph", type=float, default=7.4, help="pH for net charge (default 7.4)")
    p.add_argument("--hydropathy", action="store_true",
                   help="Scan for hydrophobic segments with a Kyte-Doolittle sliding window")
    p.add_argument("--window", type=int, default=19, help="Hydropathy window size (default 19)")
    p.add_argument("--hydro-threshold", type=float, default=1.6,
                   help="Window mean threshold for hydrophobic segments (default 1.6)")
    p.add_argument("--force", action="store_true",
                   help="Analyse sequences even if they look like nucleotides")
    p.add_argument("--csv", metavar="FILE", help="Write a CSV summary table")
    p.add_argument("--strict", action="store_true",
                   help="Exit with status 2 if any record has invalid, empty or skipped data")
    args = p.parse_args()

    if not 0 <= args.ph <= 14:
        sys.exit("Error: --ph must be between 0 and 14.")
    if args.window < 1:
        sys.exit("Error: --window must be at least 1.")

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
    print_summary(results)

    if args.csv:
        write_csv(results, args.csv, args.ph)
        print(f"\nCSV written to {args.csv}")

    if args.strict and any(r["invalid"] or r["skipped"] for r in results):
        sys.exit(2)


if __name__ == "__main__":
    main()