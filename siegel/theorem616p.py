"""The operator ``g_{I_0,e}`` of Definition 6.15', and Theorem 6.16'.

These are the statements of the **corrigendum**, not of the published paper.
The difference is the threshold ``delta_{I_0,e}`` in (D) below: the published
Definition 6.10 asked only that ``chi`` be dominant.  It vanishes identically
for ``g <= 3`` and for ``e <= 1`` at every genus, so there the two agree.

Fix a standard parabolic ``I_0`` and a cohomological degree ``e``.  Write ``S``
for the positive roots of ``L`` outside ``L_0``, ``r_0 = |S|`` and
``2rho_{I_0} = sum(S)``.  For a constituent ``mu^{d-e}_top`` of Definition 6.6' put

    chi = lambda - mu^{d-e}_top,        lambda' = chi + 2rho_{I_0}.

Then ``lambda`` lies in ``g^{mu_top}_{I_0,e}(C)`` if and only if

(D)  ``chi`` lies in ``X^*(P_0)`` and ``<chi, alpha^v> >= delta_{I_0,e}(alpha)``
     for every simple ``alpha`` outside ``I_0``;
(A)  ``lambda'`` lies in ``C_{ample,I_0}``;
(C1) for every ``k`` with ``1 <= k <= min(e, r_0)``, every ``-s_M`` in
     ``W_{I_0}(r_0-k)`` for which ``lambda' - s_M`` is dominant, and every ``j``,
     ``LR(mu^{d-e+k}_j, lambda' - s_M)`` is contained in ``C``;
(C2) for every ``j``, ``LR(mu^{d-e}_j, chi) \\ {lambda}`` is contained in ``C``.

and ``g_{I_0,e}(C)`` is the union of these over the *maximal* constituents
``mu_top``, selected by :meth:`~siegel.rootdata.SiegelData.maximal_mu`.

**Theorem 6.16'.** Assume ``p > g^2``.  If ``C`` is contained in ``C^{e+1}_van``
then ``g_{I_0,e}(C)`` is contained in ``C^e_van``.

``W_{I_0}(j)`` is the multiset of highest weights of the good filtration of
``Lambda^j(u_0^-)^v``.  Its members are ``L_0``-dominant weights ``-s_M`` with
``M`` a subset of ``S`` of size ``j``, but not every such weight occurs: it is
determined by the character, and computed in
:meth:`~siegel.rootdata.SiegelData.good_filtration_hw`.

Conditions (D), (A) and the Littlewood--Richardson sets of (C1)--(C2) do not
depend on ``C``.  Only the final memberships do.  :meth:`compile` exploits this:
it evaluates everything ``C``-independent once and returns, for each weight, the
set of weights whose membership in ``C`` would suffice.  Deciding the operator
then costs one subset test, which is what makes iterating to a fixpoint cheap.
"""

from .rootdata import add, pair, sub


class Theorem616p:
    """The operator of Definition 6.15', compiled against a fixed ``(I_0, e)``."""

    def __init__(self, data):
        self.data = data

    # ------------------------------------------------------------------
    def context(self, T, e):
        """Everything depending on ``(I_0, e)`` alone."""
        d = self.data
        T = tuple(T)
        phi_I, S = d.buildRoots(T)
        zero = tuple(0 for _ in range(d.g))
        tworho = zero
        for x in S:
            tworho = add(tworho, x)
        r0 = len(S)

        # W_{I_0}(r_0 - k) for the sizes (C1) can ask for.  Stored as the
        # weights -s_M themselves, deduplicated: (C1) quantifies over the
        # distinct constituents, so the multiplicities play no role here.
        W = {}
        for k in range(1, min(e, r0) + 1):
            W[k] = sorted(set(d.good_filtration_hw(T, r0 - k)))

        return {
            "T": T, "e": e, "S": S, "r0": r0,
            "tworho": tworho, "zero": zero,
            "delta": d.delta(T, e),
            "W": W,
        }

    # ------------------------------------------------------------------
    def requirements(self, ctx, lam):
        """Sets of weights witnessing ``lam``, one per admissible base point.

        Returns a list of ``frozenset``s.  ``lam`` is admitted as soon as any one
        of them is contained in ``C``; an empty list means no base point passes
        the ``C``-independent conditions (D) and (A), so ``lam`` can never be
        admitted at this ``(I_0, e)``.  An empty frozenset in the list means
        ``lam`` is admitted unconditionally.
        """
        d = self.data
        e = ctx["e"]
        out = []
        for mu_top in d.maximal_mu(d.d - e):
            req = self._requirement(ctx, lam, mu_top)
            if req is not None:
                out.append(req)
        return out

    def _requirement(self, ctx, lam, mu_top):
        d = self.data
        e = ctx["e"]
        chi = sub(lam, mu_top)

        # (D).  Membership of chi in X^*(P_0) is the vanishing on I_0; together
        # with the threshold outside I_0 it makes chi dominant, which (C2) needs.
        delta = ctx["delta"]
        T = ctx["T"]
        for i in range(1, d.g):
            v = pair(chi, d.Delta_L_vec[i - 1])
            if T[i - 1] != 0:
                if v != 0:
                    return None
            elif v < delta[i]:
                return None

        lamp = add(chi, ctx["tworho"])

        # (A)
        if not d.ample(T, lamp):
            return None

        # C^d_van is everything, so at e + 1 >= d there is nothing left to ask.
        if e + 1 >= d.d:
            return frozenset()

        need = set()

        # (C1)
        for k in range(1, min(e, ctx["r0"]) + 1):
            table = d.mu[d.d - e + k - 1]
            for minus_s_M in ctx["W"][k]:
                nu_M = add(lamp, minus_s_M)
                if not d.Ldominant(nu_M):
                    continue  # the graded piece vanishes outright
                for mu_j in table:
                    need.update(d.LR(mu_j, nu_M))

        # (C2).  lambda occurs once, at j = top, and is excluded there only.
        for mu_j in d.mu[d.d - e - 1]:
            for eta in d.LR(mu_j, chi):
                if eta != lam:
                    need.add(eta)

        return frozenset(need)

    # ------------------------------------------------------------------
    def compile(self, T, e, weights):
        """Precompute the requirement sets for every weight at ``(I_0, e)``.

        Returns a list of ``(lam, [req, ...])``, omitting weights that no base
        point admits.
        """
        ctx = self.context(T, e)
        out = []
        for lam in weights:
            reqs = self.requirements(ctx, lam)
            if reqs:
                out.append((lam, reqs))
        return out

    def admits(self, ctx, lam, Cvan):
        """Decide the operator directly, without compiling."""
        e = ctx["e"]
        known = Cvan[e + 1] if e + 1 < self.data.d else None
        for req in self.requirements(ctx, lam):
            if known is None or req <= known:
                return True
        return False
