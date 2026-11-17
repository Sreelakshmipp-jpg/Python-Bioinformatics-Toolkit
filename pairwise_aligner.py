#!/usr/bin/env python3
"""
pairwise_aligner.py - Pairwise sequence alignment (standard library only).

  * Global alignment  (Needleman-Wunsch) and local alignment (Smith-Waterman)
  * Affine gap penalties (Gotoh algorithm): a gap of length k costs
        gap_open + (k - 1) * gap_extend
  * DNA / RNA scoring with match/mismatch (IUPAC ambiguity codes score 0)
  * Protein scoring with BLOSUM62
  * Input from one FASTA file (first two records, or chosen with --ids),
    two FASTA files, raw strings (--seq1/--seq2), or all-vs-all comparison
  * Every alignment is re-scored from the aligned strings as a built-in self-check

Usage
  python pairwise_aligner.py seqs.fasta
  python pairwise_aligner.py seqs.fasta --ids Seq1 Seq3 --mode local
  python pairwise_aligner.py a.fasta b.fasta --gap-open 10 --gap-extend 1
  python pairwise_aligner.py --seq1 ACGTTGCA --seq2 ACGTGCA
  python pairwise_aligner.py seqs.fasta --all-vs-all --out results.txt

Penalties are given as POSITIVE numbers (e.g. --gap-open 5 means a score of -5).
"""

import argparse
import sys

# ---------------------------------------------------------------------------
# Scoring matrices
# ---------------------------------------------------------------------------

BLOSUM62_TEXT = """
      A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A     4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
R    -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
N    -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
D    -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
C     0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
Q    -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
E    -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
G     0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
H    -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
I    -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
L    -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
K    -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
M    -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
F    -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
P    -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
S     1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
T     0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
W    -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
Y    -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
V     0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
B    -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
Z    -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
X     0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
*    -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
"""


def load_blosum62():
    lines = BLOSUM62_TEXT.strip().splitlines()
    cols = lines[0].split()
    matrix = {}
    for line in lines[1:]:
        parts = line.split()
        matrix[parts[0]] = {c: int(v) for c, v in zip(cols, parts[1:])}
    return matrix


NT_SYMBOLS = "ACGTURYSWKMBDHVN"
NT_UNAMBIGUOUS = set("ACGTU")
PROTEIN_ONLY = set("EFILPQZXOJ")   # letters that are never nucleotide symbols
PROTEIN_REMAP = {"U": "X", "O": "X", "J": "X"}  # valid protein letters absent from BLOSUM62


def nucleotide_matrix(match, mismatch):
    """match/mismatch for unambiguous bases; any pair involving an IUPAC
    ambiguity code scores 0 (neutral)."""
    matrix = {}
    for x in NT_SYMBOLS:
        matrix[x] = {}
        for y in NT_SYMBOLS:
            if x in NT_UNAMBIGUOUS and y in NT_UNAMBIGUOUS:
                matrix[x][y] = match if x == y else -mismatch
            else:
                matrix[x][y] = 0.0
    return matrix


# ---------------------------------------------------------------------------
# Input handling
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
            chunks.append("".join(line.split()))
    if header is not None:
        yield header, "".join(chunks)


def read_fasta_file(path):
    try:
        with open(path) as fh:
            records = [(h.split()[0] if h else "", s.upper()) for h, s in parse_fasta(fh)]
    except OSError as exc:
        sys.exit(f"Error: {exc}")
    except ValueError as exc:
        sys.exit(f"Error in {path}: {exc}")
    if not records:
        sys.exit(f"Error: no FASTA records found in {path}")
    return records


def detect_kind(seq):
    """Return 'nucleotide' or 'protein'. Same rules as fasta_analyzer.py."""
    chars = set(seq)
    if "U" in chars and "T" not in chars and chars <= set(NT_SYMBOLS):
        return "nucleotide"
    letters = [c for c in seq if c.isalpha()]
    if not letters:
        return "nucleotide"
    core = sum(1 for c in letters if c in "ACGTUN") / len(letters)
    has_protein_only = any(c in PROTEIN_ONLY for c in letters)
    if core >= 0.75 or (not has_protein_only and len(letters) >= 20):
        return "nucleotide"
    return "protein"


def prepare_sequences(records, forced):
    """Validate and normalise a list of (id, sequence). Returns (records, kind, notes).
    Invalid characters are an error here (removing them would shift positions)."""
    notes = []
    cleaned = []
    for name, raw in records:
        seq = raw.replace("-", "").replace(".", "")
        if len(seq) != len(raw):
            notes.append(f"{name}: gap characters removed before alignment")
        if not seq:
            sys.exit(f"Error: sequence '{name}' is empty.")
        cleaned.append((name, seq))

    if forced == "auto":
        kinds = {detect_kind(s) for _, s in cleaned}
        if len(kinds) > 1:
            sys.exit("Error: sequences look like a mixture of nucleotide and protein. "
                     "Use --type dna, rna or protein to force a type.")
        kind = kinds.pop()
    else:
        kind = "protein" if forced == "protein" else "nucleotide"

    out = []
    if kind == "nucleotide":
        has_u = [("U" in s) for _, s in cleaned]
        has_t = [("T" in s) for _, s in cleaned]
        for (name, s), u, t in zip(cleaned, has_u, has_t):
            if u and t:
                sys.exit(f"Error: sequence '{name}' contains both T and U.")
        if any(has_u) and any(has_t):
            notes.append("RNA and DNA sequences mixed: U treated as T for scoring")
        for name, s in cleaned:
            bad = sorted(set(s) - set(NT_SYMBOLS))
            if bad:
                sys.exit(f"Error: sequence '{name}' has invalid nucleotide character(s): {', '.join(bad)}. "
                         f"Alignment positions would shift if they were removed; please fix the input.")
            if any(has_u) and any(has_t):
                s = s.replace("U", "T")
            out.append((name, s))
    else:
        valid = set(load_blosum62())
        for name, s in cleaned:
            if s.endswith("*"):
                s = s[:-1]
            remapped = sum(s.count(c) for c in PROTEIN_REMAP)
            if remapped:
                notes.append(f"{name}: {remapped} U/O/J residue(s) scored as X (not in BLOSUM62)")
                for old, new in PROTEIN_REMAP.items():
                    s = s.replace(old, new)
            bad = sorted(set(s) - valid)
            if bad:
                sys.exit(f"Error: sequence '{name}' has invalid protein character(s): {', '.join(bad)}.")
            if not s:
                sys.exit(f"Error: sequence '{name}' is empty after cleaning.")
            out.append((name, s))
    return out, kind, notes


# ---------------------------------------------------------------------------
# Alignment (Gotoh affine-gap dynamic programming)
# ---------------------------------------------------------------------------

def align(a, b, matrix, mode, gap_open, gap_ext):
    """Align sequences a and b.

    gap_open / gap_ext are NEGATIVE scores; a gap of length k scores
    gap_open + (k - 1) * gap_ext. mode is 'global' or 'local'.
    Returns a dict, or None if a local alignment has no positive score.

    Three states per cell: M (a[i] aligned to b[j]), X (a[i] against a gap),
    Y (b[j] against a gap). Only two score rows are kept; traceback pointers
    are stored in compact bytearrays (1 = from M, 2 = from X, 3 = from Y, 0 = start).
    """
    m, n = len(a), len(b)
    local = mode == "local"
    NEG = float("-inf")

    pM = [bytearray(n + 1) for _ in range(m + 1)]
    pX = [bytearray(n + 1) for _ in range(m + 1)]
    pY = [bytearray(n + 1) for _ in range(m + 1)]

    prevM = [NEG] * (n + 1)
    prevX = [NEG] * (n + 1)
    prevY = [NEG] * (n + 1)
    if not local:
        prevM[0] = 0.0
        for j in range(1, n + 1):
            if j == 1:
                prevY[j], pY[0][j] = gap_open, 1
            else:
                prevY[j], pY[0][j] = prevY[j - 1] + gap_ext, 3

    best, best_i, best_j = 0.0, 0, 0
    for i in range(1, m + 1):
        curM = [NEG] * (n + 1)
        curX = [NEG] * (n + 1)
        curY = [NEG] * (n + 1)
        if not local:
            if i == 1:
                curX[0], pX[i][0] = gap_open, 1
            else:
                curX[0], pX[i][0] = prevX[0] + gap_ext, 2
        row = matrix[a[i - 1]]
        pMi, pXi, pYi = pM[i], pX[i], pY[i]
        for j in range(1, n + 1):
            # M: diagonal move
            vm, vx, vy = prevM[j - 1], prevX[j - 1], prevY[j - 1]
            if vm >= vx and vm >= vy:
                v, p = vm, 1
            elif vx >= vy:
                v, p = vx, 2
            else:
                v, p = vy, 3
            if local and v <= 0:
                v, p = 0.0, 0
            curM[j] = v + row[b[j - 1]]
            pMi[j] = p
            if local and curM[j] > best:
                best, best_i, best_j = curM[j], i, j

            # X: a[i] against a gap (move down)
            vm, vx, vy = prevM[j] + gap_open, prevX[j] + gap_ext, prevY[j] + gap_open
            if vm >= vx and vm >= vy:
                curX[j], pXi[j] = vm, 1
            elif vx >= vy:
                curX[j], pXi[j] = vx, 2
            else:
                curX[j], pXi[j] = vy, 3

            # Y: b[j] against a gap (move right)
            vm, vx, vy = curM[j - 1] + gap_open, curX[j - 1] + gap_open, curY[j - 1] + gap_ext
            if vm >= vx and vm >= vy:
                curY[j], pYi[j] = vm, 1
            elif vy >= vx:
                curY[j], pYi[j] = vy, 3
            else:
                curY[j], pYi[j] = vx, 2
        prevM, prevX, prevY = curM, curX, curY

    if local:
        if best <= 0:
            return None
        score, i, j, state = best, best_i, best_j, 1
        end_i, end_j = best_i, best_j
    else:
        options = [(prevM[n], 1), (prevX[n], 2), (prevY[n], 3)]
        score, state = max(options, key=lambda t: t[0])  # ties: M, then X, then Y
        i, j = m, n
        end_i, end_j = m, n

    ra, rb = [], []
    while True:
        if state == 1:
            ra.append(a[i - 1])
            rb.append(b[j - 1])
            p = pM[i][j]
            i -= 1
            j -= 1
            if local and p == 0:
                break
            state = p
        elif state == 2:
            ra.append(a[i - 1])
            rb.append("-")
            state = pX[i][j]
            i -= 1
        else:
            ra.append("-")
            rb.append(b[j - 1])
            state = pY[i][j]
            j -= 1
        if not local and i == 0 and j == 0:
            break

    return {
        "score": score,
        "aligned_a": "".join(reversed(ra)),
        "aligned_b": "".join(reversed(rb)),
        "start_a": i + 1, "end_a": end_i,
        "start_b": j + 1, "end_b": end_j,
    }


def rescore(al_a, al_b, matrix, gap_open, gap_ext):
    """Recompute an alignment's score from the aligned strings (independent check)."""
    total, prev = 0.0, None
    for x, y in zip(al_a, al_b):
        if x == "-" or y == "-":
            g = "a" if x == "-" else "b"
            total += gap_open if g != prev else gap_ext
            prev = g
        else:
            total += matrix[x][y]
            prev = None
    return total


def alignment_stats(al_a, al_b, matrix, protein):
    ident = pos = mism = gaps = opens = 0
    prev = None
    for x, y in zip(al_a, al_b):
        if x == "-" or y == "-":
            gaps += 1
            g = "a" if x == "-" else "b"
            if g != prev:
                opens += 1
            prev = g
            continue
        prev = None
        unambiguous = (x != "X" and x != "*") if protein else (x in NT_UNAMBIGUOUS)
        if x == y and unambiguous:
            ident += 1
        elif protein and matrix[x][y] > 0:
            pos += 1
        else:
            mism += 1
    length = len(al_a)
    return {"length": length, "identities": ident, "positives": pos, "mismatches": mism,
            "gaps": gaps, "gap_opens": opens,
            "identity_pct": 100.0 * ident / length if length else 0.0,
            "similarity_pct": 100.0 * (ident + pos) / length if length else 0.0}


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def match_symbol(x, y, matrix, protein):
    if x == "-" or y == "-":
        return " "
    unambiguous = (x != "X" and x != "*") if protein else (x in NT_UNAMBIGUOUS)
    if x == y and unambiguous:
        return "|"
    if protein and matrix[x][y] > 0:
        return ":"
    return "."


def format_alignment(res, name_a, name_b, matrix, protein, width=60):
    al_a, al_b = res["aligned_a"], res["aligned_b"]
    cur_a, cur_b = res["start_a"] - 1, res["start_b"] - 1
    la, lb = name_a[:12], name_b[:12]
    indent = " " * 21
    lines = []
    for k in range(0, len(al_a), width):
        sa, sb = al_a[k:k + width], al_b[k:k + width]
        mid = "".join(match_symbol(x, y, matrix, protein) for x, y in zip(sa, sb))
        na, nb = len(sa) - sa.count("-"), len(sb) - sb.count("-")
        s_a, e_a = (cur_a + 1, cur_a + na) if na else (cur_a, cur_a)
        s_b, e_b = (cur_b + 1, cur_b + nb) if nb else (cur_b, cur_b)
        lines.append(f"{la:<12} {s_a:>7} {sa} {e_a:>7}")
        lines.append(f"{indent}{mid}")
        lines.append(f"{lb:<12} {s_b:>7} {sb} {e_b:>7}")
        lines.append("")
        cur_a += na
        cur_b += nb
    return lines


def fmt(x):
    return f"{x:g}"


def describe_scoring(kind, params):
    if kind == "protein":
        return (f"BLOSUM62, gap open {fmt(-params['gap_open'])}, gap extend {fmt(-params['gap_ext'])}")
    return (f"match +{fmt(params['match'])}, mismatch -{fmt(params['mismatch'])}, "
            f"gap open {fmt(-params['gap_open'])}, gap extend {fmt(-params['gap_ext'])}; "
            f"IUPAC ambiguity codes score 0")


def run_pair(rec_a, rec_b, matrix, protein, args, params):
    res = align(rec_a[1], rec_b[1], matrix, args.mode, params["gap_open"], params["gap_ext"])
    if res is None:
        return None
    check = rescore(res["aligned_a"], res["aligned_b"], matrix, params["gap_open"], params["gap_ext"])
    if abs(check - res["score"]) > 1e-9:
        raise RuntimeError(f"Internal consistency check failed: DP score {res['score']} "
                           f"!= re-scored alignment {check}")
    res["stats"] = alignment_stats(res["aligned_a"], res["aligned_b"], matrix, protein)
    return res


def report_pair(rec_a, rec_b, res, kind, params, matrix, args):
    protein = kind == "protein"
    algo = ("global, Needleman-Wunsch" if args.mode == "global" else "local, Smith-Waterman")
    out = [f"Pairwise alignment ({algo}, affine gaps)",
           f"Sequence 1 : {rec_a[0]} ({len(rec_a[1])} {'aa' if protein else 'nt'})",
           f"Sequence 2 : {rec_b[0]} ({len(rec_b[1])} {'aa' if protein else 'nt'})",
           f"Scoring    : {describe_scoring(kind, params)}",
           "Gap cost   : gap of length k scores open + (k-1) * extend"]
    if res is None:
        out += ["", "No positive-scoring local alignment found.", ""]
        return out
    st = res["stats"]
    out.append(f"Score      : {fmt(res['score'])}")
    if args.mode == "local":
        out.append(f"Aligned    : seq1 {res['start_a']}-{res['end_a']}, seq2 {res['start_b']}-{res['end_b']}")
    out.append(f"Length     : {st['length']} columns")
    out.append(f"Identity   : {st['identities']}/{st['length']} ({st['identity_pct']:.1f}%)")
    if protein:
        out.append(f"Similarity : {st['identities'] + st['positives']}/{st['length']} "
                   f"({st['similarity_pct']:.1f}%)  [identities + positive BLOSUM62 scores]")
    out.append(f"Mismatches : {st['mismatches']}")
    out.append(f"Gaps       : {st['gaps']} positions in {st['gap_opens']} gap(s)")
    legend = "| identity" + (", : similar" if protein else "") + ", . mismatch"
    out.append(f"Legend     : {legend}")
    out.append("")
    out += format_alignment(res, rec_a[0], rec_b[0], matrix, protein)
    return out


def report_all_vs_all(records, matrix, protein, args, params):
    n = len(records)
    ident = [[None] * n for _ in range(n)]
    score = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            res = run_pair(records[i], records[j], matrix, protein, args, params)
            if res is not None:
                ident[i][j] = ident[j][i] = res["stats"]["identity_pct"]
                score[i][j] = score[j][i] = res["score"]
        ident[i][i] = 100.0

    names = [r[0][:10] for r in records]
    algo = "global" if args.mode == "global" else "local"
    out = [f"All-vs-all pairwise comparison ({algo}, affine gaps)",
           f"Scoring: {describe_scoring('protein' if protein else 'nucleotide', params)}",
           "Identity = identical columns / alignment columns (%)", ""]
    header = " " * 11 + "".join(f"{nm:>11}" for nm in names)
    out.append("% identity")
    out.append(header)
    for i in range(n):
        cells = "".join(f"{'-' if ident[i][j] is None else format(ident[i][j], '.1f'):>11}" for j in range(n))
        out.append(f"{names[i]:<11}{cells}")
    out += ["", "Alignment score", header]
    for i in range(n):
        cells = "".join(
            f"{'' if i == j else ('-' if score[i][j] is None else fmt(score[i][j])):>11}" for j in range(n))
        out.append(f"{names[i]:<11}{cells}")
    out.append("")
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def positive(value):
    x = float(value)
    if x < 0:
        raise argparse.ArgumentTypeError("penalties must be given as positive numbers")
    return x


def main():
    p = argparse.ArgumentParser(description="Pairwise sequence alignment (global or local, affine gaps).")
    p.add_argument("fasta", nargs="*", help="One FASTA file (first two records) or two FASTA files (first record of each)")
    p.add_argument("--seq1", help="First sequence as a raw string")
    p.add_argument("--seq2", help="Second sequence as a raw string")
    p.add_argument("--ids", nargs=2, metavar=("ID1", "ID2"), help="Choose two records by ID from a single FASTA file")
    p.add_argument("--all-vs-all", action="store_true", help="Compare every pair of records in the file")
    p.add_argument("--mode", choices=["global", "local"], default="global")
    p.add_argument("--type", choices=["auto", "dna", "rna", "protein"], default="auto")
    p.add_argument("--match", type=positive, help="Match score for DNA/RNA (default 2)")
    p.add_argument("--mismatch", type=positive, help="Mismatch penalty for DNA/RNA (default 3)")
    p.add_argument("--gap-open", type=positive, help="Gap open penalty (default 5 DNA/RNA, 11 protein)")
    p.add_argument("--gap-extend", type=positive, help="Gap extend penalty (default 2 DNA/RNA, 1 protein)")
    p.add_argument("--max-cells", type=int, default=10_000_000,
                   help="Refuse alignments larger than this many DP cells (default 10,000,000)")
    p.add_argument("--out", metavar="FILE", help="Also write the report to FILE")
    args = p.parse_args()

    # ---- choose records ------------------------------------------------------
    if args.seq1 or args.seq2:
        if not (args.seq1 and args.seq2):
            sys.exit("Error: give both --seq1 and --seq2.")
        records = [("seq1", args.seq1.upper()), ("seq2", args.seq2.upper())]
    elif len(args.fasta) == 2:
        records = [read_fasta_file(args.fasta[0])[0], read_fasta_file(args.fasta[1])[0]]
    elif len(args.fasta) == 1:
        records = read_fasta_file(args.fasta[0])
        if args.ids:
            by_id = {}
            for rid, seq in records:
                by_id.setdefault(rid, seq)
            missing = [i for i in args.ids if i not in by_id]
            if missing:
                sys.exit(f"Error: ID(s) not found: {', '.join(missing)}. Available: {', '.join(by_id)}")
            records = [(i, by_id[i]) for i in args.ids]
        elif not args.all_vs_all:
            if len(records) < 2:
                sys.exit("Error: need at least two records to align.")
            if len(records) > 2:
                print(f"Note: file has {len(records)} records; aligning the first two "
                      f"(use --ids or --all-vs-all).", file=sys.stderr)
            records = records[:2]
        elif len(records) < 2:
            sys.exit("Error: need at least two records for --all-vs-all.")
    else:
        p.error("provide one FASTA file, two FASTA files, or --seq1 and --seq2")

    records, kind, notes = prepare_sequences(records, args.type)
    protein = kind == "protein"

    # ---- scoring --------------------------------------------------------------
    if protein:
        if args.match is not None or args.mismatch is not None:
            notes.append("--match/--mismatch ignored for protein (BLOSUM62 is used)")
        matrix = load_blosum62()
        params = {"gap_open": -(args.gap_open if args.gap_open is not None else 11.0),
                  "gap_ext": -(args.gap_extend if args.gap_extend is not None else 1.0)}
    else:
        match = args.match if args.match is not None else 2.0
        mismatch = args.mismatch if args.mismatch is not None else 3.0
        if match <= 0:
            sys.exit("Error: --match must be greater than 0.")
        matrix = nucleotide_matrix(match, mismatch)
        params = {"match": match, "mismatch": mismatch,
                  "gap_open": -(args.gap_open if args.gap_open is not None else 5.0),
                  "gap_ext": -(args.gap_extend if args.gap_extend is not None else 2.0)}

    # ---- size guard -----------------------------------------------------------
    lens = [len(s) for _, s in records]
    worst = max(lens) * sorted(lens)[-2] if len(lens) > 1 else 0
    if worst > args.max_cells:
        sys.exit(f"Error: largest alignment needs {worst:,} DP cells (limit {args.max_cells:,}). "
                 f"Raise --max-cells if you have the memory and patience (pure Python).")

    # ---- run ---------------------------------------------------------------------
    lines = []
    for note in notes:
        lines.append(f"Note: {note}")
    if args.all_vs_all and not (args.seq1 or len(args.fasta) == 2 or args.ids):
        lines += report_all_vs_all(records, matrix, protein, args, params)
    else:
        res = run_pair(records[0], records[1], matrix, protein, args, params)
        lines += report_pair(records[0], records[1], res, kind, params, matrix, args)

    text = "\n".join(lines)
    print(text)
    if args.out:
        try:
            with open(args.out, "w") as fh:
                fh.write(text + "\n")
        except OSError as exc:
            sys.exit(f"Error writing {args.out}: {exc}")
        print(f"Report written to {args.out}")


if __name__ == "__main__":
    main()