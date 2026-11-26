"""Plot helpers. Colour is never the only carrier of meaning: significance also changes marker fill."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import config as C

BLUE = "#1f5fa8"
ORANGE = "#c8561b"
GREY = "#8a8a8a"


def forest(ax, labels, est, lo, hi, q, is_subset, title, xlabel):
    """One dot per group. Filled = FDR below alpha. Square marker = labeled subset (outside FDR)."""
    y = np.arange(len(labels))
    for i in range(len(labels)):
        sig = (not np.isnan(q[i])) and q[i] < C.FDR_ALPHA
        marker = "s" if is_subset[i] else "o"
        colour = ORANGE if is_subset[i] else BLUE
        ax.plot([lo[i], hi[i]], [y[i], y[i]], color=colour, lw=1.2, alpha=0.8)
        ax.plot(est[i], y[i], marker, color=colour, mfc=colour if sig else "white", mew=1.4, ms=6)
    ax.axvline(0, color=GREY, lw=0.8, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)


def save(fig, name):
    C.FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(C.FIG_DIR / f"{name}.png", dpi=200, bbox_inches="tight")
    fig.savefig(C.FIG_DIR / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
