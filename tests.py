#!/usr/bin/env sage -python
"""Self-checks for the implementation.

    PYTHONPATH=. sage -python tests.py

Everything here is internal: each test either checks the code against an
independent computation inside Sage, or checks a property the mathematics
guarantees.  No reference data is consulted, so the suite is meaningful on a
fresh checkout.
"""

import random
import sys

from sage.all import Integer

from siegel import SiegelData, VanishingEngine
from siegel.rootdata import add
from siegel.psmall import psmall

failures = []


def check(name, ok, detail=""):
    print("[%s] %s%s" % ("PASS" if ok else "FAIL", name,
                         (" -- " + detail) if detail else ""))
    if not ok:
        failures.append(name)


# ----------------------------------------------------------------------
# The hypothesis of Theorem 6.16' is enforced
# ----------------------------------------------------------------------
raised = False
try:
    SiegelData(3, 7)          # 7 < 3^2
except ValueError:
    raised = True
check("p > g^2 is enforced", raised)
check("p > g^2 can be bypassed deliberately",
      SiegelData(3, 7, check=False).p == 7)


# ----------------------------------------------------------------------
# The labelling convention
#
# Which of the two labellings is in force is not visible from dominance --
# both call a weight dominant when it decreases -- so it is pinned here by
# the anchor the paper prints under both labels.
# ----------------------------------------------------------------------
random.seed(20260728)
for g, p in [(2, 5), (3, 11), (4, 17)]:
    data = SiegelData(g, p)
    tau = data.changeConvention

    sample = [tuple(Integer(x) for x in
                    sorted((random.randint(-9, 3) for _ in range(g)), reverse=True))
              for _ in range(50)]
    sample.append(tuple(Integer(0) for _ in range(g)))
    check("g=%d: the twist is an involution" % g,
          all(tau(tau(w)) == w for w in sample))
    check("g=%d: the twist preserves dominance" % g,
          all(data.Ldominant(tau(w)) for w in sample))

    # rho = (g-1, g-3, ..., 1-g)/1, the half-sum of the positive roots of L,
    # doubled to stay integral.  tau is the diagram automorphism, so it fixes it.
    two_rho = tuple(Integer(g - 1 - 2 * i) for i in range(g))
    check("g=%d: the twist fixes 2rho" % g, tau(two_rho) == two_rho)

    # Sym^2 std has intrinsic highest weight (2,0,...,0) and label (0,...,0,-2).
    anchor = tuple(Integer(0) for _ in range(g - 1)) + (Integer(-2),)
    check("g=%d: the anchor Sym^2 std is labelled (0,...,0,-2)" % g,
          data.mu[0] == [anchor],
          "mu[0] = %s" % (data.mu[0],))

data3 = SiegelData(3, 11)
check("g=3, n=3: the two maximal constituents carry the published labels",
      sorted(data3.mu[2]) == [(-1, -1, -4), (0, -3, -3)],
      "mu[2] = %s" % (sorted(data3.mu[2]),))


# ----------------------------------------------------------------------
# Littlewood--Richardson against Sage's Weyl character ring
# ----------------------------------------------------------------------
random.seed(20260728)
for g, p in [(2, 5), (3, 11), (4, 17)]:
    data = SiegelData(g, p)
    bad = 0
    for _ in range(120):
        mu = tuple(Integer(x) for x in
                   sorted((random.randint(-9, 3) for _ in range(g)), reverse=True))
        nu = tuple(Integer(x) for x in
                   sorted((random.randint(-9, 3) for _ in range(g)), reverse=True))
        viaLR = set(data.LR(mu, nu))
        viaWCR = set(tuple(Integer(x[i]) for i in range(g))
                     for x in (data.A(list(mu)) * data.A(list(nu))).monomial_coefficients())
        if viaLR != viaWCR:
            bad += 1
    check("g=%d: LR rule agrees with the Weyl character ring" % g, bad == 0,
          "120 random pairs")


# ----------------------------------------------------------------------
# The good filtration reproduces the character it decomposes
# ----------------------------------------------------------------------
for g, p in [(2, 5), (3, 11), (4, 17), (5, 29)]:
    data = SiegelData(g, p)
    ok = True
    for T in data.parabolics():
        phi_I, S = data.buildRoots(T)
        zero = tuple(Integer(0) for _ in range(g))
        for j in range(len(S) + 1):
            rebuilt = {}
            for nu in data.good_filtration_hw(T, j):
                for w, c in data.nabla_L0_character(T, nu).items():
                    rebuilt[w] = rebuilt.get(w, 0) + c
            original = {}
            for M in data.subsets_of_size(S, j):
                s_M = zero
                for x in M:
                    s_M = add(s_M, x)
                w = tuple(-c for c in s_M)
                original[w] = original.get(w, 0) + 1
            if rebuilt != original:
                ok = False
    check("g=%d: good filtration reproduces the character of Lambda^j(u_0^-)^v" % g, ok)


# ----------------------------------------------------------------------
# The threshold behaves as Definition 6.10' says
# ----------------------------------------------------------------------
for g, p, bound in [(2, 5, 0), (3, 11, 0), (4, 17, 1), (5, 29, 2)]:
    data = SiegelData(g, p)
    worst = 0
    zero_low = True
    for T in data.parabolics():
        for e in range(data.d):
            dl = data.delta(T, e)
            m = max([dl[i] for i in range(1, g) if T[i - 1] == 0], default=0)
            worst = max(worst, m)
            if e <= 1 and m != 0:
                zero_low = False
    check("g=%d: max delta = %d as expected" % (g, bound), worst == bound,
          "got %d" % worst)
    check("g=%d: delta vanishes for e <= 1" % g, zero_low)


# ----------------------------------------------------------------------
# Properties of a saturated computation
# ----------------------------------------------------------------------
data = SiegelData(2, 7)
X = VanishingEngine(data)
X.saturate(-50, 0)

check("Cvan is increasing in the degree",
      all(X.Cvan[k] <= X.Cvan[k + 1] for k in range(data.d - 1)))
check("every result is a dominant weight in the window",
      all(all(w[i] >= w[i + 1] for i in range(data.g - 1))
          and min(w) >= -50 and max(w) <= 0
          for k in range(data.d) for w in X.Cvan[k]))
check("saturation is a fixpoint", not X.sweep(-50, 0))

Y = VanishingEngine(SiegelData(2, 7))
Y.saturate(-50, 0)
check("computation is deterministic",
      all(X.Cvan[k] == Y.Cvan[k] for k in range(data.d)))

# Order of the sweep must not change the fixpoint.
Z = VanishingEngine(SiegelData(2, 7))
for _ in range(12):
    for T in reversed(Z.data.parabolics()):
        for e in reversed(range(Z.data.d)):
            Z.apply(T, e, -50, 0)
check("fixpoint is independent of sweep order",
      all(X.Cvan[k] == Z.Cvan[k] for k in range(data.d)))


# ----------------------------------------------------------------------
# Certificates really certify
# ----------------------------------------------------------------------
from siegel.analysis import certificate

sample = sorted(X.Cvan[0])[:40] + sorted(X.Cvan[data.d - 1])[:40]
ok = True
for w in sample:
    k = X.concentration(w)
    certs = certificate(X, w, -50, 0)
    if not certs:
        continue
    if min(c["e"] for c in certs) != k:
        ok = False
check("certificates match the recorded concentration degree", ok,
      "%d weights" % len(sample))


# ----------------------------------------------------------------------
# Round trip through disk
# ----------------------------------------------------------------------
import tempfile, os

with tempfile.TemporaryDirectory() as tmp:
    X.save(directory=tmp)
    R = VanishingEngine(SiegelData(2, 7))
    R.load(directory=tmp)
    check("save/load round trip", all(R.Cvan[k] == X.Cvan[k] for k in range(data.d)))


# ----------------------------------------------------------------------
# p-small weights
# ----------------------------------------------------------------------
res = psmall(SiegelData(2, 7), 42, write=False)
check("p-small weights are dominant",
      all(w[0] >= w[1] for w in res), "%d weights" % len(res))


print("\n%s" % ("all checks passed" if not failures else "FAILURES: %s" % failures))
sys.exit(1 if failures else 0)
