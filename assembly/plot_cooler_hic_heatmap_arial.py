#!/usr/bin/env python3

import argparse
import numpy as np
import cooler
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

def setup_font():
    """
    Prefer Arial. If Arial is unavailable, use Liberation Sans, then DejaVu Sans.
    The output PDF uses editable TrueType fonts when possible.
    """
    available_fonts = {f.name for f in fm.fontManager.ttflist}

    if "Arial" in available_fonts:
        use_font = "Arial"
    elif "Liberation Sans" in available_fonts:
        use_font = "Liberation Sans"
    else:
        use_font = "DejaVu Sans"

    mpl.rcParams["font.family"] = use_font
    mpl.rcParams["font.sans-serif"] = [use_font]
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    mpl.rcParams["svg.fonttype"] = "none"
    mpl.rcParams["axes.unicode_minus"] = False

    print(f"[INFO] Using matplotlib font: {use_font}")

def fetch_matrix(c, chrom1, chrom2=None, use_balance=True):
    if chrom2 is None:
        chrom2 = chrom1

    try:
        mat = c.matrix(balance=use_balance, sparse=False).fetch(chrom1, chrom2)
    except Exception as e:
        if use_balance:
            print(f"[WARNING] Balanced fetch failed for {chrom1}, {chrom2}: {e}. Fallback to raw.")
            mat = c.matrix(balance=False, sparse=False).fetch(chrom1, chrom2)
        else:
            raise

    return mat

def main():
    setup_font()

    parser = argparse.ArgumentParser()
    parser.add_argument("--cool", required=True)
    parser.add_argument("--chromsizes", required=True)
    parser.add_argument("--out_pdf", required=True)
    parser.add_argument("--out_png", required=True)
    parser.add_argument("--title", default="")
    parser.add_argument("--vmax_percentile", type=float, default=99.5)
    parser.add_argument("--figsize", default="8,8")
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--cmap", default="Reds")
    parser.add_argument("--raw", action="store_true", help="Use raw counts instead of balanced matrix.")
    args = parser.parse_args()

    use_balance = not args.raw

    chroms = []
    with open(args.chromsizes) as f:
        for line in f:
            if not line.strip():
                continue
            chroms.append(line.strip().split()[0])

    c = cooler.Cooler(args.cool)

    bin_counts = []
    for chrom in chroms:
        mat = fetch_matrix(c, chrom, use_balance=use_balance)
        bin_counts.append(mat.shape[0])

    total_bins = sum(bin_counts)
    full = np.zeros((total_bins, total_bins), dtype=float)

    starts = np.cumsum([0] + bin_counts[:-1])
    ends = np.cumsum(bin_counts)

    for i, chr_i in enumerate(chroms):
        for j, chr_j in enumerate(chroms):
            block = fetch_matrix(c, chr_i, chr_j, use_balance=use_balance)
            full[starts[i]:ends[i], starts[j]:ends[j]] = block

    full = np.nan_to_num(full, nan=0.0, posinf=0.0, neginf=0.0)

    if use_balance:
        plotmat = full
        cbar_label = "Normalized Hi-C contact intensity"
    else:
        plotmat = np.log1p(full)
        cbar_label = "log1p(contact count)"

    positive = plotmat[plotmat > 0]

    if positive.size > 0:
        vmax = np.percentile(positive, args.vmax_percentile)
        if vmax <= 0:
            vmax = positive.max()
    else:
        vmax = 1

    w, h = [float(x) for x in args.figsize.split(",")]

    fig, ax = plt.subplots(figsize=(w, h))

    im = ax.imshow(
        plotmat,
        cmap=args.cmap,
        vmin=0,
        vmax=vmax,
        interpolation="nearest",
        origin="upper"
    )

    centers = [(s + e) / 2 for s, e in zip(starts, ends)]

    ax.set_xticks(centers)
    ax.set_xticklabels(chroms, rotation=90, fontsize=9)
    ax.set_yticks(centers)
    ax.set_yticklabels(chroms, fontsize=9)

    for b in ends[:-1]:
        ax.axhline(b - 0.5, color="black", linewidth=0.15, alpha=0.65)
        ax.axvline(b - 0.5, color="black", linewidth=0.15, alpha=0.65)

    if args.title and args.title != "NONE":
        ax.set_title(args.title, fontsize=12)

    ax.set_xlabel("")
    ax.set_ylabel("")

    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label(cbar_label, fontsize=9)
    cbar.ax.tick_params(labelsize=9)

    ax.tick_params(axis="both", which="major", labelsize=9)

    plt.tight_layout()
    plt.savefig(args.out_pdf, dpi=args.dpi)
    plt.savefig(args.out_png, dpi=args.dpi)
    plt.close()

if __name__ == "__main__":
    main()
