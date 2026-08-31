"""Inspecting a computation: certificates, thresholds, and where results come from.

Nothing here changes a result.  These are tools for asking *why* a weight is in
the output, *which* parabolic produced it, and *what* the threshold of
Definition 6.10' looks like at a given genus.
"""

from sage.all import Integer

from .theorem616p import Theorem616p


def certificate(engine, char, kmin, kmax):
    """Every ``(I_0, e, mu_top)`` at which Theorem 6.16' admits ``char``.

    Returns a list of dicts, ordered by ``e``.  Empty if the weight is not in the
    output.  Each entry records the base point ``mu_top``, the derived weights
    ``chi`` and ``lambda'``, and how many weights the conditions (C1)--(C2)
    required -- a rough measure of how hard that particular certificate was.
    """
    data = engine.data
    w = tuple(Integer(x) for x in char)
    op = Theorem616p(data)
    out = []
    for T in data.parabolics():
        for e in range(data.d):
            ctx = op.context(T, e)
            known = engine.Cvan[e + 1] if e + 1 < data.d else None
            for mu_top in data.mu[data.d - e - 1]:
                req = op._requirement(ctx, w, mu_top)
                if req is None:
                    continue
                if known is not None and not req <= known:
                    continue
                chi = tuple(a - b for a, b in zip(w, mu_top))
                lamp = tuple(a + b for a, b in zip(chi, ctx["tworho"]))
                out.append({
                    "I_0": T, "e": e, "mu_top": mu_top,
                    "chi": chi, "lambda_prime": lamp,
                    "conditions": len(req),
                })
    out.sort(key=lambda r: (r["e"], r["I_0"]))
    return out


def explain(engine, char, kmin, kmax):
    """Print a readable certificate for ``char``."""
    data = engine.data
    w = tuple(Integer(x) for x in char)
    k = engine.concentration(w)
    if k is None:
        print("%s: no vanishing known." % (list(w),))
        return []
    print("%s: H^i(X, nabla(lambda)(-D)) = 0 for all i > %d." % (list(w), k))
    certs = certificate(engine, w, kmin, kmax)
    if not certs:
        print("  (obtained at a degree above the one shown, or loaded from file)")
        return certs
    print("  admitted by Theorem 6.16' at %d parabolic/degree pairs:" % len(certs))
    for c in certs[:12]:
        print("    I_0 = %-10s e = %d   mu_top = %-16s lambda' = %-16s %d conditions"
              % (c["I_0"], c["e"], str(list(c["mu_top"])),
                 str(list(c["lambda_prime"])), c["conditions"]))
    if len(certs) > 12:
        print("    ... and %d more" % (len(certs) - 12))
    return certs


def contributions(engine, kmin, kmax):
    """How many weights each parabolic can produce, degree by degree.

    Counts the weights that Theorem 6.16' admits at ``(I_0, e)`` against the final
    ``Cvan``.  A weight admitted at several parabolics is counted at each, so the
    rows do not sum to the totals -- the point is to see which parabolics carry
    the computation and which are redundant.
    """
    data = engine.data
    op = Theorem616p(data)
    weights = data.buildWeights(kmin, kmax)
    table = {}
    for T in data.parabolics():
        row = []
        for e in range(data.d):
            known = engine.Cvan[e + 1] if e + 1 < data.d else None
            n = 0
            for lam, reqs in op.compile(T, e, weights):
                if known is None or any(r <= known for r in reqs):
                    n += 1
            row.append(n)
        table[T] = row
    return table


def show_contributions(engine, kmin, kmax):
    data = engine.data
    table = contributions(engine, kmin, kmax)
    print("%-14s" % "I_0", end="")
    for e in range(data.d):
        print("%8s" % ("e=%d" % e), end="")
    print()
    for T in data.parabolics():
        print("%-14s" % (T,), end="")
        for n in table[T]:
            print("%8d" % n, end="")
        print()
    return table


def delta_table(data, emax=None):
    """The threshold ``delta_{I_0,e}`` of Definition 6.10' at every parabolic.

    Prints the maximum of ``delta_{I_0,e}`` over the simple roots outside
    ``I_0``; ``0`` means the threshold imposes nothing there.
    """
    emax = data.d - 1 if emax is None else emax
    print("max_alpha delta_{I_0,e}(alpha)   (g = %d)" % data.g)
    print("%-14s" % "I_0", end="")
    for e in range(emax + 1):
        print("%6s" % ("e=%d" % e), end="")
    print()
    out = {}
    for T in data.parabolics():
        row = []
        print("%-14s" % (T,), end="")
        for e in range(emax + 1):
            dl = data.delta(T, e)
            m = max([dl[i] for i in range(1, data.g) if T[i - 1] == 0], default=0)
            row.append(m)
            print("%6d" % m, end="")
        print()
        out[T] = row
    return out


def compare_primes(g, primes, kmin, kmax, check=True):
    """Saturate at several primes and tabulate the counts.

    Useful for seeing where the output stabilises in ``p``: the orbital
    ``p``-closeness condition loosens as ``p`` grows and eventually stops binding.
    """
    from .rootdata import SiegelData
    from .engine import VanishingEngine

    rows = {}
    for p in primes:
        data = SiegelData(g, p, check=check)
        eng = VanishingEngine(data)
        eng.saturate(kmin, kmax)
        rows[p] = [len(eng.Cvan[k]) for k in range(data.d)]
    d = g * (g + 1) // 2
    print("%-8s" % "p", end="")
    for k in range(d):
        print("%9s" % ("[0:%d]" % k), end="")
    print()
    for p in primes:
        print("%-8d" % p, end="")
        for n in rows[p]:
            print("%9d" % n, end="")
        print()
    return rows
