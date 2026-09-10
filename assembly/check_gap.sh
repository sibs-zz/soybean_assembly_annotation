#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# 检测 FASTA 中每个 scaffold 上的 gap（连续 N/n）
#
# 用法：
#   bash check_gap.sh genome.fa
#   bash check_gap.sh genome.fa.gz
#
# 输出（TSV）：
#   scaffold  gap_index  start  end  length
# 末尾打印汇总。
# 坐标：1-based，闭区间 [start, end]
# ============================================================

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <genome.fa[.gz]>" >&2
    exit 1
fi

FA="$1"
[[ -s "$FA" ]] || { echo "ERROR: file missing or empty: $FA" >&2; exit 1; }

python3 - "$FA" <<'PY'
import sys
import gzip
import re

path = sys.argv[1]
opener = gzip.open if path.endswith(".gz") else open

gap_re = re.compile(r"[Nn]+")

print("scaffold\tgap_index\tstart\tend\tlength")

n_scaff = 0
n_gap = 0
total_gap_bp = 0
scaff_with_gap = 0
details = []  # for summary

def flush(name, seq):
    global n_scaff, n_gap, total_gap_bp, scaff_with_gap
    if not name:
        return
    n_scaff += 1
    gaps = []
    for m in gap_re.finditer(seq):
        start = m.start() + 1  # 1-based
        end = m.end()          # inclusive
        length = end - start + 1
        gaps.append((start, end, length))
    if gaps:
        scaff_with_gap += 1
    for i, (start, end, length) in enumerate(gaps, 1):
        n_gap += 1
        total_gap_bp += length
        print(f"{name}\t{i}\t{start}\t{end}\t{length}")
    details.append((name, len(seq), gaps))

name = None
chunks = []
with opener(path, "rt") as fh:
    for line in fh:
        if line.startswith(">"):
            flush(name, "".join(chunks))
            name = line[1:].strip().split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    flush(name, "".join(chunks))

print("", file=sys.stderr)
print("========== summary ==========", file=sys.stderr)
print(f"scaffolds          : {n_scaff}", file=sys.stderr)
print(f"scaffolds_with_gap : {scaff_with_gap}", file=sys.stderr)
print(f"gap_count          : {n_gap}", file=sys.stderr)
print(f"total_gap_bp       : {total_gap_bp}", file=sys.stderr)

# per-scaffold short summary on stderr
print("", file=sys.stderr)
print("scaffold\tseq_len\tgap_n\tgap_bp\tmax_gap", file=sys.stderr)
for name, seq_len, gaps in details:
    if not gaps:
        continue
    gap_bp = sum(g[2] for g in gaps)
    max_gap = max(g[2] for g in gaps)
    print(f"{name}\t{seq_len}\t{len(gaps)}\t{gap_bp}\t{max_gap}", file=sys.stderr)
PY
