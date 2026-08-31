"""Plots of the vanishing results, for ``g = 2`` and ``g = 3``.

Weights are drawn in the convention of the rest of the package (dominant =
decreasing), coloured by the degree above which the cohomology is known to
vanish.
"""

import matplotlib.pyplot as plt

COLORS = ["#1b1b1b", "#1f6feb", "#d1242f", "#1a7f37", "#8250df", "#bf8700"]


def plot(engine, title=None, ax=None, figsize=None):
    """Scatter the output of ``engine``, coloured by concentration degree."""
    data = engine.data
    res = engine.convert()

    if data.g == 2:
        if ax is None:
            plt.figure(figsize=figsize or (11, 8))
            ax = plt.gca()
        for k in range(data.d):
            if res[k]:
                x, y = zip(*res[k])
                ax.scatter(x, y, s=14, color=COLORS[k % len(COLORS)],
                           label="concentrated in degrees [0 : %d]  (%d)" % (k, len(res[k])))
        ax.set_xlabel(r"$k_1$")
        ax.set_ylabel(r"$k_2$")
        ax.set_title(title or
                     "Vanishing results for the Siegel threefold, $p = %d$" % data.p)
        ax.legend(loc="upper left", frameon=False)
        ax.grid(alpha=0.15)
        return ax

    if data.g == 3:
        if ax is None:
            fig = plt.figure(figsize=figsize or (13, 10))
            ax = fig.add_subplot(projection="3d")
        for k in range(data.d):
            if res[k]:
                x, y, z = zip(*res[k])
                # depthshade dims distant points, which greys the palette
                # exactly where the picture is densest.
                ax.scatter3D(x, y, z, s=10, alpha=0.75, depthshade=False,
                             color=COLORS[k % len(COLORS)],
                             label="degrees [0 : %d]  (%d)" % (k, len(res[k])))
        # Shrink the box inside the axes, or the k_3 label falls off the canvas.
        ax.set_box_aspect(None, zoom=0.88)
        ax.set_xlabel(r"$k_1$")
        ax.set_ylabel(r"$k_2$")
        ax.set_zlabel(r"$k_3$")
        ax.set_title(title or
                     "Vanishing results for the Siegel variety of genus 3, $p = %d$"
                     % data.p)
        ax.legend(loc="upper left", frameon=False)
        return ax

    raise ValueError("Plotting handles g = 2 and g = 3; got g = %s" % data.g)


def plot_slice(engine, fixed, axis=2, ax=None):
    """For ``g = 3``: the slice of the output at one fixed coordinate.

    ``axis`` is the index held fixed (0, 1 or 2) and ``fixed`` its value.  Useful
    for reading a three-dimensional picture one plane at a time.
    """
    data = engine.data
    if data.g != 3:
        raise ValueError("plot_slice is for g = 3")
    res = engine.convert()
    keep = [i for i in range(3) if i != axis]
    if ax is None:
        plt.figure(figsize=(9, 7))
        ax = plt.gca()
    for k in range(data.d):
        pts = [(w[keep[0]], w[keep[1]]) for w in res[k] if w[axis] == fixed]
        if pts:
            ax.scatter(*zip(*pts), s=22, color=COLORS[k % len(COLORS)],
                       label="degrees [0 : %d]  (%d)" % (k, len(pts)))
    ax.set_xlabel(r"$k_%d$" % (keep[0] + 1))
    ax.set_ylabel(r"$k_%d$" % (keep[1] + 1))
    ax.set_title(r"Slice $k_%d = %d$, $p = %d$" % (axis + 1, fixed, data.p))
    ax.legend(loc="upper left", frameon=False)
    ax.grid(alpha=0.15)
    return ax


def plot_primes(engines, ax=None):
    """Counts per degree, one line per prime, from a dict ``{p: engine}``."""
    if ax is None:
        plt.figure(figsize=(9, 6))
        ax = plt.gca()
    degrees = 0
    for i, (p, eng) in enumerate(sorted(engines.items())):
        counts = [len(eng.Cvan[k]) for k in range(eng.data.d)]
        degrees = max(degrees, len(counts))
        ax.plot(range(len(counts)), counts, marker="o",
                color=COLORS[i % len(COLORS)], label="$p = %d$" % p)
    # k is a degree; fractional ticks would be meaningless.
    ax.set_xticks(range(degrees))
    ax.set_xlabel("degree $k$")
    ax.set_ylabel(r"weights with $H^{>k} = 0$")
    ax.set_title("Output as a function of $p$")
    ax.legend(frameon=False)
    ax.grid(alpha=0.15)
    return ax
