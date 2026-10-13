"""Shared figure style for 3_experiments/make_figures.py and compare_report.py.

Palette: the first three categorical slots of a colour-vision-deficiency
validated palette (blue, orange, aqua; worst adjacent CVD Delta E 9.2). Aqua
is below 3:1 contrast on white, so every bar carries a direct value label.
Text is always in neutral ink, never in a series colour.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]   # slot 1, 2, 3
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e4e3df"
SURFACE = "#fcfcfb"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.edgecolor": GRID, "axes.labelcolor": INK_2, "axes.titlecolor": INK,
    "xtick.color": INK_2, "ytick.color": INK_2, "text.color": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "axes.grid.axis": "y", "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.axisbelow": True, "legend.frameon": False, "savefig.dpi": 200,
})


def grouped_bars(ax, labels, series, names, value_fmt="{:.1f}", above=None):
    """Bars of len(series) per label; 'above' = optional text above each group."""
    k = len(series)
    width = 0.8 / k
    gap = 0.02
    for j, (vals, name) in enumerate(zip(series, names)):
        xs = [i - 0.4 + width * (j + 0.5) for i in range(len(labels))]
        bars = ax.bar(xs, vals, width - gap, color=SERIES[j], label=name, linewidth=0)
        for b, v in zip(bars, vals):
            if v == v:   # not NaN
                ax.text(b.get_x() + b.get_width() / 2, b.get_height(), value_fmt.format(v),
                        ha="center", va="bottom", fontsize=6.5, color=INK_2)
    top = max(max(v for v in s if v == v) for s in series)
    if above:
        for i, t in enumerate(above):
            group_top = max(s[i] for s in series if s[i] == s[i])
            if t:
                ax.text(i, group_top + top * 0.07, t, ha="center", va="bottom", fontsize=8,
                        color=INK, fontweight="bold")
    ax.set_ylim(0, top * 1.38)          # empty band on top keeps the legend clear of data
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels)
    ax.legend(loc="upper center", ncol=k, fontsize=8)


def stacked_hbars(ax, labels, parts, names, total_fmt="{:.1f} s"):
    """Horizontal stacked bars, one per label; parts[j][i] = segment j of bar i."""
    left = [0.0] * len(labels)
    ys = list(range(len(labels)))[::-1]
    for j, (vals, name) in enumerate(zip(parts, names)):
        ax.barh(ys, vals, left=left, height=0.6, color=SERIES[j], label=name,
                edgecolor=SURFACE, linewidth=1.5)
        left = [a + b for a, b in zip(left, vals)]
    for y, tot in zip(ys, left):
        ax.text(tot, y, "  " + total_fmt.format(tot), va="center", fontsize=8, color=INK)
    ax.set_yticks(ys)
    ax.set_yticklabels(labels)
    ax.grid(axis="x"); ax.grid(axis="y", visible=False)
    ax.set_xlim(0, max(left) * 1.18)
    ax.legend(loc="lower right", fontsize=8)


def save(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"))
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)
    return path.with_suffix(".png")
