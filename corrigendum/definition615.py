"""Definition 6.15 *as printed in the published paper*, for comparison only.

.. warning::

   This module is **not** part of what the repository deposits.  The statement
   it implements does not hold as printed -- that is what the corrigendum is
   about -- and nothing here should be used to derive a vanishing result.  It
   exists for one purpose: to let the comparison tables of the corrigendum,
   which set the published operator against the corrected one, be recomputed
   from source rather than taken on trust.  The deposited statement is
   Theorem 6.16', implemented in ``siegel/theorem616p.py``, and it is the only
   one the paper's conclusions may rest on.

   Keep this folder out of any import path used for real computation, and do
   not promote this module into the ``siegel`` package.

What is implemented here is the printed Definition 6.15 with the coding
defect of the accompanying code repaired -- the corrigendum compares against
the *statement*, not against the shipped implementation, so the slip is
fixed before the comparison is drawn:

1. ``p > g^2`` is enforced (:class:`~siegel.rootdata.SiegelData` does this);
2. the subset ``M = emptyset``, i.e. ``k = r_0``, is not skipped.

Two further points, both faithful to the printed text rather than to the
corrected one:

* the base point.  Definition 6.6 speaks of "the" highest weight of
  ``Lambda^{d-e} Sym^2 std_L``, and there is in general no such thing: the
  constituents are pairwise incomparable, so each is maximal and none of them
  is "the" one.  The printed definition therefore has to be *read*, and the
  reading is stated in the corrigendum: the union over the maximal
  constituents, which is the convention Definition 6.6' uses as well.

  This is the charitable reading -- it gives the printed operator every base
  point it could have meant, so a weight it still fails to admit is one no
  reading admits -- and it is the one that isolates what the tables measure,
  since it leaves the threshold as the only difference between the two
  operators.

  The alternative, fixing a single constituent, was what this module used to
  do, via ``data.mu[n-1][0]``.  That list is built from
  ``monomial_coefficients()`` and carries no defined order, so the printed
  operator's output depended on the enumeration order of the plethysm routine.
  At ``g = 5``, ``n = 12`` the two constituents give a dominant and a
  non-dominant ``chi``, so that accident decided admission;
* the degeneration conditions quantify over **all** ``binom(d, n)`` translated
  ``T``-weights of ``Lambda^n Sym^2 std_L``, not over the Littlewood--Richardson
  constituents of the good filtrations.  This is the difference that makes the
  printed operator ask more questions than Definition 6.15' does.

The threshold ``delta_{I_0,e}`` is absent, as it is from the printed statement:
condition (D) asks only that ``chi`` be dominant.  That absence is what the
comparison measures, in both directions.

The interface matches :class:`siegel.theorem616p.Theorem616p` -- ``context``,
``compile`` and ``admits`` -- so that one driver can run either operator and any
difference in the output is a difference between the operators.
"""

from siegel.rootdata import add, pair, sub


class Definition615:
    """The published operator ``g_{I_0,e}``, compiled against a fixed ``(I_0, e)``."""

    name = "Thm. 6.16"

    def __init__(self, data):
        self.data = data
        self._all_weights = None

    # ------------------------------------------------------------------
    def all_T_weights(self, n):
        """Every ``T``-weight of ``Lambda^n Sym^2 std_L``, in the paper's convention.

        ``SiegelData`` keeps only the constituent highest weights, those being
        what Definition 6.6' needs; the printed Definition 6.15 quantifies over
        the full weight set, which is rebuilt here rather than in the package.
        """
        d = self.data
        if self._all_weights is None:
            Omega1 = d.A([2] + [0] * (d.g - 1))
            self._all_weights = [
                [d.changeConvention(x)
                 for x in Omega1.exterior_power(k).weight_multiplicities()]
                for k in range(1, d.d + 1)
            ]
        return self._all_weights[n - 1]

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
        simple_I0 = [d.Delta_L_vec[i - 1] for i in range(1, d.g) if T[i - 1] != 0]
        return {
            "T": T, "e": e, "S": S, "r0": len(S),
            "tworho": tworho, "zero": zero,
            "simple_I0": simple_I0,
        }

    def _local_at(self, ctx, lam, mu_top):
        """``chi`` at one base point, if ``lam`` passes there the conditions
        independent of ``C``.

        These are (D) as printed -- ``chi`` dominant on ``S`` and trivial on
        ``I_0``, with no threshold -- and (A).  ``None`` if either fails.
        """
        d = self.data
        chi = sub(lam, mu_top)

        for alpha in ctx["S"]:
            if pair(chi, alpha) < 0:
                return None
        for alpha in ctx["simple_I0"]:
            if pair(chi, alpha) != 0:
                return None

        if not d.ample(ctx["T"], add(chi, ctx["tworho"])):
            return None
        return chi

    def _local(self, ctx, lam):
        """Every ``(mu_top, chi)`` at which ``lam`` passes (D) and (A).

        The union runs over the maximal constituents, which is how the printed
        Definition 6.6 is read here; see the module docstring.  An empty list
        means no base point passes, in which case no amount of prior knowledge
        can admit ``lam`` at this ``(I_0, e)``.
        """
        d = self.data
        e = ctx["e"]
        out = []
        for mu_top in d.maximal_mu(d.d - e):
            chi = self._local_at(ctx, lam, mu_top)
            if chi is not None:
                out.append((mu_top, chi))
        return out

    def compile(self, T, e, weights):
        """The weights that (D) and (A) leave standing at ``(I_0, e)``.

        Returned as ``(lam, None)``, the ``None`` asking the driver to call
        :meth:`admits` each sweep: the degeneration conditions are a conjunction
        with an early exit, and re-deciding one is cheaper than tabulating the
        set of every weight it consults.  Above ``e = d - 1`` there is nothing
        left to ask and the pair carries an empty requirement instead.
        """
        ctx = self.context(T, e)
        unconditional = e + 1 >= self.data.d
        out = []
        for lam in weights:
            if self._local(ctx, lam):
                out.append((lam, [frozenset()] if unconditional else None))
        return out

    def admits(self, ctx, lam, Cvan):
        """True if ``lam`` is in the image of the printed operator at ``(I_0, e)``."""
        d = self.data
        e = ctx["e"]

        bases = self._local(ctx, lam)
        if not bases:
            return False

        # C^d_van is everything, so at e + 1 >= d there is nothing left to ask.
        if e + 1 >= d.d:
            return True

        # Admitted as soon as one base point carries lam through.
        return any(self._degeneration(ctx, lam, Cvan, mu_top, chi)
                   for mu_top, chi in bases)

    def _degeneration(self, ctx, lam, Cvan, mu_top, chi):
        """The degeneration conditions of Definition 6.15, at one base point."""
        d = self.data
        e = ctx["e"]
        S = ctx["S"]
        r0 = ctx["r0"]
        known = Cvan[e + 1]

        # Over all T-weights and all subsets M of S of size r_0 - k.
        # Repair (2): k = r_0, i.e. M = emptyset, is included.
        for k in range(e + 1):
            if r0 - k < 0:
                break
            for M in d.subsets_of_size(S, r0 - k):
                s_M = ctx["zero"]
                for x in M:
                    s_M = add(s_M, x)
                table = self.all_T_weights(d.d - e + k)
                for mu_j in table:
                    # At k = 0, M = S the base point contributes lambda itself,
                    # which the printed condition excludes just as (C2) does.
                    if k == 0 and mu_j == mu_top:
                        continue
                    value = add(add(sub(chi, s_M), ctx["tworho"]), mu_j)
                    if d.Ldominant(value) and value not in known:
                        return False
        return True
