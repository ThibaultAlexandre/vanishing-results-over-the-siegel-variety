#!/usr/bin/env sage -python
"""Reproduce the numerical results of Section 5 of the corrigendum.

    PYTHONPATH=.. sage -python reproduce.py [g ...]

with no argument running every case.  Two things are recomputed:

* **the count tables** (Sec. 5.2), setting the published Theorem 6.16 -- read as
  Definition 6.15 with the coding defect repaired -- against Theorem 6.16'
  of the corrigendum, degree by degree, at ``g = 2, 3, 4``;
* **the lost weight at g = 5** (Prop. C.18), every assertion of the proposition
  checked one at a time.

Each printed number is compared against the value the corrigendum carries, and
the script ends in ``ALL FIGURES REPRODUCED`` or lists what differs.  Nothing is
written to disk and nothing here is imported by the package.

Why a window control.  Every count comes from saturating an operator over a
finite box ``[-W,0]^g``.  The operators ask whether certain weights already lie
in ``C^{e+1}_van``; a weight outside the box never does, so a candidate near the
face of the box is rejected for a reason that has nothing to do with the
statement being tested, and **every count is a lower bound**.  The bias is not
the same for the two operators -- Definition 6.15 quantifies over all
``binom(d,n)`` translated T-weights and 6.15' only over the Littlewood--Richardson
constituents, so they consult different weights at different depths -- so a
difference at one window is not yet a difference of operators.  Each operator is
therefore saturated at two windows ``W1 < W2`` and reported only on the largest
sub-box ``[-B,0]^g`` on which its two runs agree.

Runtime: about two minutes in all, nearly every second of it the ``g = 4``
case; ``g = 2`` and ``g = 3`` take seconds, and the ``g = 5`` check is immediate,
testing one weight and needing no saturation.
"""

import sys
import time

sys.path.insert(0, "..")

from siegel import SiegelData, Theorem616p                       # noqa: E402
from siegel.rootdata import add, pair, sub                      # noqa: E402

from definition615 import Definition615                         # noqa: E402

# (g, p, W1, W2, B): the windows the corrigendum's tables were computed at, and
# the stable sub-box they yield.  W2 was pushed up until the radius stopped
# growing; B is an output of the run, checked below rather than assumed.
CASES = [
    (2, 7, 60, 80, 59),
    (3, 11, 26, 34, 23),
    (4, 31, 18, 24, 13),
]

# Section 5.2, one row per operator, then the two set differences.
EXPECTED = {
    (2, 7): {
        "6.16": [962, 1113, 1162],
        "6.16'": [962, 1114, 1162],
        "gained": [0, 1, 0],
        "lost": [0, 0, 0],
    },
    (3, 11): {
        "6.16": [763, 948, 972, 986, 1026, 1347],
        "6.16'": [763, 950, 978, 993, 1037, 1356],
        "gained": [0, 2, 6, 7, 11, 9],
        "lost": [0, 0, 0, 0, 0, 0],
    },
    (4, 31): {
        "6.16": [99, 130, 134, 139, 141, 141, 141, 141, 142, 490],
        "6.16'": [99, 133, 142, 149, 151, 152, 154, 155, 163, 504],
        "gained": [0, 3, 8, 10, 10, 11, 13, 14, 21, 14],
        "lost": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    },
}

FAILURES = []


def check(label, got, want):
    ok = got == want
    if not ok:
        FAILURES.append("%s: got %s, corrigendum says %s" % (label, got, want))
    return ok


# ----------------------------------------------------------------------
# One driver, both operators.  Any difference in the output is then a
# difference between the operators and not between their drivers.
#
# The conditions of either statement split in two: a part that does not depend
# on what is already known -- (D) and (A), and for Theorem 6.16' the
# Littlewood--Richardson sets as well -- and a membership test against
# ``C^{e+1}_van``.  Each operator evaluates the first once per ``(I_0, e)`` and
# hands back, for every surviving weight, either the set whose membership would
# suffice (Theorem 6.16', whose conditions are a single such set per base point)
# or ``None``, meaning "ask me again each sweep" (the printed statement, whose
# condition is a conjunction with an early exit, cheaper to re-decide than to
# tabulate).  Sweeping then costs a subset test per candidate.
# ----------------------------------------------------------------------
def saturate(data, operator, W, seed=None, max_rounds=64):
    """Iterate ``operator`` over ``[-W,0]^g`` to a fixpoint; return ``Cvan``.

    ``seed`` is a ``Cvan`` known to be contained in the fixpoint -- the run over
    a narrower box, say.  Admission is monotone in what is already known, so a
    weight admitted over a sub-box is admitted here; starting from it saves
    sweeps and cannot change the result.
    """
    Cvan = [set() for _ in range(data.d)]
    if seed is not None:
        for k in range(data.d):
            Cvan[k] |= set(w for w in seed[k] if min(w) >= -W)
    weights = data.buildWeights(-W, 0)

    contexts, compiled = {}, {}
    for T in data.parabolics():
        for e in range(data.d):
            contexts[(T, e)] = operator.context(T, e)
            compiled[(T, e)] = operator.compile(T, e, weights)

    for _ in range(max_rounds):
        grew = False
        for T in data.parabolics():
            for e in range(data.d):
                known = Cvan[e + 1] if e + 1 < data.d else None
                here = Cvan[e]
                gained = []
                for lam, reqs in compiled[(T, e)]:
                    if lam in here:
                        continue
                    if reqs is None:
                        if operator.admits(contexts[(T, e)], lam, Cvan):
                            gained.append(lam)
                    elif known is None or any(req <= known for req in reqs):
                        gained.append(lam)
                for lam in gained:
                    for k in range(e, data.d):
                        Cvan[k].add(lam)
                grew = grew or bool(gained)
        if not grew:
            break
    return Cvan


def restrict(Cvan, B):
    """``Cvan`` cut down to the sub-box ``[-B,0]^g``."""
    return [set(w for w in s if min(w) >= -B) for s in Cvan]


def stability_radius(narrow, wide, W1):
    """Largest ``B <= W1`` on which two saturations agree at every degree."""
    best = -1
    for B in range(W1 + 1):
        if all(a == b for a, b in zip(restrict(narrow, B), restrict(wide, B))):
            best = B
        else:
            break
    return best


def windowed(data, operator, W1, W2):
    """Saturate at both windows; return the wide run and its stability radius."""
    narrow = saturate(data, operator, W1)
    wide = saturate(data, operator, W2, seed=narrow)
    return wide, stability_radius(narrow, wide, W1)


# ----------------------------------------------------------------------
# Sec. 5.2, the count tables
# ----------------------------------------------------------------------
def table(g, p, W1, W2, B_expected):
    print("\n" + "=" * 72)
    print("g = %d, p = %d   windows [-%d,0]^%d and [-%d,0]^%d" % (g, p, W1, g, W2, g))
    print("=" * 72)
    data = SiegelData(g, p)
    d = data.d

    runs, radii = {}, {}
    for label, operator in (("6.16", Definition615(data)),
                            ("6.16'", Theorem616p(data))):
        t = time.time()
        runs[label], radii[label] = windowed(data, operator, W1, W2)
        print("  %-6s radius %2d   %7.1fs" % (label, radii[label], time.time() - t))

    B = min(radii.values())
    print("  common stability radius B = %d" % B)
    check("g=%d p=%d stability radius" % (g, p), B, B_expected)

    S = {label: restrict(runs[label], B) for label in runs}
    counts = {label: [len(s) for s in S[label]] for label in S}
    gained = [len(S["6.16'"][k] - S["6.16"][k]) for k in range(d)]
    lost = [len(S["6.16"][k] - S["6.16'"][k]) for k in range(d)]

    print("\n  counts on the stable sub-box [-%d,0]^%d, by degree:\n" % (B, g))
    print("  %-14s%s" % ("e", "".join("%7d" % k for k in range(d))))
    for label in ("6.16", "6.16'"):
        print("  %-14s%s" % ("Thm. " + label,
                             "".join("%7d" % c for c in counts[label])))
    print("  %-14s%s" % ("gained", "".join("%7d" % x for x in gained)))
    print("  %-14s%s" % ("lost", "".join("%7d" % x for x in lost)))

    want = EXPECTED[(g, p)]
    for label in ("6.16", "6.16'"):
        check("g=%d p=%d row %s" % (g, p, label), counts[label], want[label])
    check("g=%d p=%d gained" % (g, p), gained, want["gained"])
    check("g=%d p=%d lost" % (g, p), lost, want["lost"])


# ----------------------------------------------------------------------
# Sec. 5.4, Proposition C.18: a weight lost at g = 5
# ----------------------------------------------------------------------
def lost_weight_g5():
    print("\n" + "=" * 72)
    print("g = 5, p = 29, e = d - 1 : the lost weight of Proposition C.18")
    print("=" * 72)
    g, p = 5, 29
    data = SiegelData(g, p)
    old, new = Definition615(data), Theorem616p(data)
    e = data.d - 1
    lam = (-5, -6, -7, -8, -11)
    T_empty = tuple(0 for _ in range(g - 1))

    # n = d - e = 1, so the plethysm has a single constituent and no reading
    # of Definition 6.6 enters.
    constituents = sorted(data.mu[data.d - e - 1], reverse=True)
    print("  constituents of Lambda^1 Sym^2 std_5 : %s"
          % [list(x) for x in constituents])
    check("g=5 single constituent at n=1", constituents, [(0, 0, 0, 0, -2)])

    chi = sub(lam, constituents[0])
    pv = [int(pair(chi, data.Delta_L_vec[i - 1])) for i in range(1, g)]
    lamp = add(chi, _tworho(data, T_empty))
    amp = bool(data.ample(T_empty, lamp))
    delta = [int(data.delta(T_empty, e)[i]) for i in range(1, g)]
    print("  chi = lambda - mu = %s   pairings %s   regular: %s"
          % (list(chi), pv, all(v > 0 for v in pv)))
    print("  lambda' = chi + 2rho = %s   Z_0-ample and orbitally p-close: %s"
          % (list(lamp), amp))
    print("  delta_{empty,%d} = %s" % (e, delta))
    check("g=5 witness chi", list(chi), [-5, -6, -7, -8, -9])
    check("g=5 witness regular", all(v > 0 for v in pv), True)
    check("g=5 witness lambda'", list(lamp), [-1, -4, -7, -10, -13])
    check("g=5 witness satisfies (A)", amp, True)
    check("delta_{empty,%d}" % e, delta, [2] * (g - 1))
    check("g=5 witness fails (D) at every simple root",
          all(pv[i] < delta[i] for i in range(g - 1)), True)

    # C^d_van = X*, so at e = d - 1 the degeneration conditions of both
    # definitions are discharged and admission is decided by (D) and (A).
    admitted_old = bool(old.admits(old.context(T_empty, e), lam, None))
    admitted_new = any(new.requirements(new.context(T, e), lam)
                       for T in data.parabolics())
    print("  admitted by Theorem 6.16  at I_0 = empty        : %s" % admitted_old)
    print("  admitted by Theorem 6.16' at some (I_0, mu_top) : %s" % admitted_new)
    check("g=5 witness admitted by Theorem 6.16", admitted_old, True)
    check("g=5 witness lost by Theorem 6.16'", admitted_new, False)

    # The e = 3 candidate of the paragraph after the proposition: blocked by
    # (D) everywhere, its admission by Theorem 6.16 left undecided.
    e3, lam3 = 3, (-7, -11, -13, -16, -18)
    chi3 = sub(lam3, (-2, -5, -5, -6, -6))
    lamp3 = add(chi3, _tworho(data, T_empty))
    blocked = not any(new.requirements(new.context(T, e3), lam3)
                      for T in data.parabolics())
    print("\n  e = 3 candidate %s : chi = %s, lambda' = %s"
          % (list(lam3), list(chi3), list(lamp3)))
    print("  blocked by Theorem 6.16' at every parabolic     : %s" % blocked)
    check("g=5 e=3 candidate chi", list(chi3), [-5, -6, -8, -10, -12])
    check("g=5 e=3 candidate lambda'", list(lamp3), [-1, -4, -8, -12, -16])
    check("g=5 e=3 candidate satisfies (A)", bool(data.ample(T_empty, lamp3)), True)
    check("g=5 e=3 candidate blocked by Theorem 6.16'", blocked, True)


def _tworho(data, T):
    _, S = data.buildRoots(T)
    out = tuple(0 for _ in range(data.g))
    for x in S:
        out = add(out, x)
    return out


# ----------------------------------------------------------------------
if __name__ == "__main__":
    want = [int(x) for x in sys.argv[1:]] or [2, 3, 4, 5]

    # State the labelling before printing a single weight.  Both labellings call
    # a weight dominant when it decreases, so dominance settles nothing; the
    # anchor does.
    _probe = SiegelData(3, 11)
    print("Weights are labelled as in Definition 3.10 of the paper: the twist is")
    print("tau(a_1,...,a_g) = (-a_g,...,-a_1), and Sym^2 std carries the label")
    print("mu[0] = %s, of intrinsic highest weight (2,0,...,0)."
          % (list(_probe.mu[0][0]),))
    assert _probe.mu[0] == [(0, 0, -2)], _probe.mu[0]
    assert all(_probe.changeConvention(_probe.changeConvention(w)) == w
               for w in _probe.mu[2]), "the twist is not an involution"

    for g, p, W1, W2, B in CASES:
        if g in want:
            table(g, p, W1, W2, B)
    if 5 in want:
        lost_weight_g5()

    print("\n" + "=" * 72)
    if FAILURES:
        print("%d FIGURE(S) DIFFER FROM THE CORRIGENDUM" % len(FAILURES))
        for line in FAILURES:
            print("  " + line)
        sys.exit(1)
    print("ALL FIGURES REPRODUCED")
