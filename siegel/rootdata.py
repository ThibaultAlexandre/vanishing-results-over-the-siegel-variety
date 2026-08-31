"""Root data and representation theory for the Siegel variety.

Weights of ``L = GL_g`` are stored as ``g``-tuples in the convention of the
paper, including the labelling of its Definition 3.10: the label of an
automorphic bundle is the character one induces from, so that ``pi_* L_lambda =
nabla(lambda)`` is label-preserving, and the underlying module has highest weight
the image of the label under the twist ``tau(a_1,...,a_g) = (-a_g,...,-a_1)``
(the action of ``-w_0``).  ``changeConvention`` is that twist, and it is applied
in exactly the two places where the package meets representation theory stated
intrinsically -- the plethysm constituents in ``_build_form_weights`` and the
``p``-small weights.  Everything else here is a character of a parabolic or a
weight in the root lattice and carries no twist.

Both labellings share one notion of dominance, so ``dominant means decreasing``,
``lambda_1 >= ... >= lambda_g``, holds either way and does not by itself say
which is in force: ``(-2,-8)`` is dominant, ``(-8,-2)`` is not.  The anchor that
does say it is ``Sym^2 std``, of intrinsic highest weight ``(2,0,...,0)`` and
label ``(0,...,0,-2)``; the latter is what ``self.mu[0]`` holds.  This is the
convention of the files in ``save/``.

Sage's ``WeylCharacterRing(['A', g-1])`` uses the same decreasing convention, so
the stored tuples can be handed to it directly.  Since ``-w_0`` is a ring
automorphism of the character ring (it is ``V |-> V^*``), decomposing a tensor
product here and decomposing it in Sage's convention correspond term by term.
"""

from sage.all import Integer, RootSystem, WeylCharacterRing, WeylGroup, Primes


def vec(x, g):
    """Coordinates of a weight as a plain tuple of integers."""
    return tuple(Integer(x[i]) for i in range(g))


def pair(char, alpha):
    """The pairing ``<char, alpha>``.  In types A and C roots are identified with
    coroots in the ambient space, so a single routine serves both."""
    total = 0
    for j in range(len(char)):
        total += char[j] * alpha[j]
    return total


def add(u, v):
    return tuple(a + b for a, b in zip(u, v))


def sub(u, v):
    return tuple(a - b for a, b in zip(u, v))


class SiegelData:
    """Root data, differential-form weights and the ampleness test for ``X``.

    ``X`` is a smooth projective toroidal compactification of the Siegel variety
    of genus ``g`` over ``F_p``, of dimension ``d = g(g+1)/2``, with boundary
    ``D``.  Theorem 6.16' assumes ``p > g^2``, which is enforced here; pass
    ``check=False`` only to explore the boundary of that hypothesis deliberately.
    """

    def __init__(self, g, p, check=True):
        self.g = Integer(g)
        self.p = Integer(p)
        self.d = Integer(g * (g + 1) // 2)

        if self.g < 2:
            raise ValueError("The genus g must be greater than 1")
        if self.p not in Primes():
            raise ValueError("The number %s is not prime" % p)
        if check and self.p <= self.g ** 2:
            raise ValueError(
                "Theorem 6.16' requires p > g^2; got p = %s and g^2 = %s"
                % (p, self.g ** 2)
            )

        self.L = RootSystem("A" + str(self.g - 1)).ambient_space()
        self.LG = RootSystem("C" + str(self.g)).ambient_space()
        self.WG = WeylGroup(self.LG)

        self.A = WeylCharacterRing(["A", self.g - 1])
        self.C = WeylCharacterRing(["C", self.g])
        self.phi_L = self.A.positive_roots()
        self.phi_G = self.C.positive_roots()
        self.Delta_L = self.A.simple_roots()

        self.phi_L_vec = [vec(x, self.g) for x in self.phi_L]
        self.Delta_L_vec = [vec(self.Delta_L[i], self.g) for i in range(1, self.g)]

        from sage.all import SymmetricFunctions, ZZ
        self._schur = SymmetricFunctions(ZZ).schur()

        self._build_form_weights()
        self._build_ample_tables()
        self._LR_cache = {}
        self._weights_cache = {}
        self._gl_cache = {}
        self._gf_cache = {}

    # ------------------------------------------------------------------
    # Differential forms
    # ------------------------------------------------------------------
    def _build_form_weights(self):
        """Constituent highest weights of ``Lambda^n Sym^2 std_L`` for ``n = 1..d``.

        These are the ``mu^n_j`` of Definition 6.6': the highest weights occurring
        in the plethysm decomposition -- equivalently, the constituents of a good
        filtration of ``Omega^n(log D)`` -- and *not* all of its ``T``-weights.
        They are pairwise incomparable, so in general there is no single highest
        weight among them; the operator unions over the maximal ones.
        """
        Omega1 = self.A([2] + [0] * (self.g - 1))
        self.mu = []
        self._maximal_mu = {}
        for n in range(1, self.d + 1):
            Omega = Omega1.exterior_power(n)
            self.mu.append([self.changeConvention(x)
                            for x in Omega.monomial_coefficients()])

    def maximal_mu(self, n):
        """The maximal ``mu^n_j``, in the dominance order, among the constituents.

        Definition 6.15' restricts the base point ``mu^{d-e}_{j_0}`` to a
        *maximal* constituent, and unions over those.  Maximality is what lets
        ``lambda = mu_{j_0} + chi`` be moved to the quotient end of the good
        filtration in the last step of the proof (Remark C.17); a non-maximal
        base point would put ``lambda`` strictly below another constituent
        ``mu_j + chi``, where the argument does not run.

        The constituents of ``Lambda^n Sym^2 std`` appear to be pairwise
        incomparable at every ``g`` and ``n``, so this filter is expected to be
        the identity -- it is, by enumeration, for every ``g <= 5``.  It is
        applied rather than assumed, so that the operator matches its
        definition without relying on that claim.  All the constituents share a
        total degree, so ``a <= b`` here is the comparison of partial sums.
        """
        hit = self._maximal_mu.get(n)
        if hit is not None:
            return hit

        cands = self.mu[n - 1]

        def leq(a, b):
            run = 0
            for x, y in zip(a, b):
                run += x - y
                if run > 0:
                    return False
            return True

        out = [a for a in cands if not any(b != a and leq(a, b) for b in cands)]
        self._maximal_mu[n] = out
        return out

    def changeConvention(self, char):
        """The ``-w_0`` twist: reverse and negate."""
        return tuple(-Integer(char[i]) for i in range(self.g - 1, -1, -1))

    # ------------------------------------------------------------------
    # Parabolics
    # ------------------------------------------------------------------
    def buildIprime(self, T):
        """Simple roots of ``I_0``, selected by the flag vector ``T``."""
        return [self.Delta_L[i] for i in range(1, len(T) + 1) if T[i - 1] != 0]

    def buildRoots(self, T):
        """Split the positive roots of ``L`` into those of ``L_0`` and ``S``.

        ``S`` is the complement, whose cardinality is ``r_0``.
        """
        Iprime = self.buildIprime(T)
        phi_I = []
        for i in range(1 << len(Iprime)):
            temp = sum((Iprime[j] for j in range(len(Iprime)) if (i & (1 << j))),
                       self.L.zero())
            if temp != 0:
                if self.LG([temp[k] for k in range(self.g)]) in self.phi_G:
                    phi_I.append(vec(temp, self.g))
        S = [x for x in self.phi_L_vec if x not in phi_I]
        return phi_I, S

    def parabolics(self):
        """All flag vectors ``T``, indexing the standard parabolics ``I_0``."""
        out = []
        for i in range(1 << (self.g - 1)):
            out.append(tuple(1 if (i & (1 << j)) else 0 for j in range(self.g - 1)))
        return out

    def blocks(self, T):
        """For each ``i``, the size of the block of ``{1,...,g}`` containing it."""
        sizes = []
        current = 1
        for j in range(self.g - 1):
            if T[j] != 0:
                current += 1
            else:
                sizes.append(current)
                current = 1
        sizes.append(current)
        out = []
        for s in sizes:
            out.extend([s] * s)
        return out

    def delta(self, T, e):
        """The threshold ``delta_{I_0,e}`` of Definition 6.10', as a dict on ``I``.

        With ``N_e(alpha_i, I_0) = min{e, g - |B_i| - |B_{i+1}|}`` outside ``I_0``
        and ``0`` on ``I_0``, one has ``delta = max{0, N_e - 1}``.  It vanishes
        identically for ``e <= 1`` at every genus and for ``g <= 3`` at every
        ``e``, is at most ``1`` for ``e = 2``, and is independent of ``e`` once
        ``e >= g - 2``.
        """
        block_of = self.blocks(T)
        out = {}
        for i in range(1, self.g):
            if T[i - 1] != 0:
                out[i] = Integer(0)
            else:
                N = min(Integer(e), self.g - block_of[i - 1] - block_of[i])
                out[i] = max(Integer(0), N - 1)
        return out

    # ------------------------------------------------------------------
    # Dominance
    # ------------------------------------------------------------------
    def Ldominant(self, char):
        for a in self.Delta_L_vec:
            if pair(char, a) < 0:
                return False
        return True

    def L0dominant(self, phi_I, char):
        for a in phi_I:
            if pair(char, a) < 0:
                return False
        return True

    # ------------------------------------------------------------------
    # Ampleness
    # ------------------------------------------------------------------
    def _build_ample_tables(self):
        """Group the coroots into Weyl orbits.

        The orbital ``p``-closeness condition asks that
        ``|<chi, w(alpha^v)>| <= (p-1)|<chi, alpha^v>|`` for every ``w`` and every
        positive ``alpha``.  Scanning pairs ``(w, alpha)`` costs ``|W| x |Phi^+|``
        pairings per weight.  But ``max_w |<chi, w(alpha^v)>|`` depends only on
        the Weyl orbit of ``alpha^v``, so it suffices to hold, for each orbit, the
        orbit itself and the positive coroots lying in it: the cost drops to the
        size of the union of the orbits, ``2g^2`` pairings, independent of ``|W|``.
        """
        seen = {}
        orbits = []
        for alpha in self.phi_G:
            cr_amb = alpha.associated_coroot()
            cr = vec(cr_amb, self.g)
            if cr in seen:
                orbits[seen[cr]][1].append(cr)
                continue
            orb = set()
            for w in self.WG:
                orb.add(vec(w.action(cr_amb), self.g))
            idx = len(orbits)
            orbits.append((sorted(orb), [cr]))
            for v in orb:
                seen.setdefault(v, idx)
            seen[cr] = idx
        self._orbits = [(tuple(o), tuple(pos)) for o, pos in orbits]

        self._noncompact = [
            vec(alpha.associated_coroot(), self.g)
            for alpha in self.phi_G
            if vec(alpha, self.g) not in self.phi_L_vec
        ]
        self._ample_cache = {}

    def ample(self, T, char):
        """``True`` if ``L_char`` is ``D``-ample on the flag bundle of type ``I_0``.

        This is the sufficient subset of the ``D``-ample cone furnished by
        Theorem 5.12 and Remark 6.17 -- orbitally ``p``-close and ``Z_0``-ample.
        The ``D``-ample cone itself is not reconstructible from the paper, so
        every count derived from this is a lower bound on the true output.
        """
        T = tuple(T)
        key = (T, char)
        hit = self._ample_cache.get(key)
        if hit is None:
            hit = self._ample(T, char)
            self._ample_cache[key] = hit
        return hit

    def _ample(self, T, char):
        # Z_0-ampleness first: a handful of pairings, and the usual reason to fail.
        for cr in self._noncompact:
            if pair(char, cr) >= 0:
                return False
        for i in range(1, self.g):
            if T[i - 1] == 0:
                if pair(char, self.Delta_L_vec[i - 1]) <= 0:
                    return False
        # Orbital p-closeness, one pass per Weyl orbit of coroots.
        bound = self.p - 1
        for orb, positives in self._orbits:
            best = 0
            for v in orb:
                a = pair(char, v)
                if a < 0:
                    a = -a
                if a > best:
                    best = a
            if best == 0:
                continue
            least = 0
            for cr in positives:
                a = pair(char, cr)
                if a < 0:
                    a = -a
                if a and (least == 0 or a < least):
                    least = a
            if least and best > bound * least:
                return False
        return True

    # ------------------------------------------------------------------
    # Littlewood--Richardson
    # ------------------------------------------------------------------
    def LR(self, mu, nu):
        """``{eta : c^eta_{mu,nu} != 0}``, the constituents of ``nabla(mu) (x) nabla(nu)``.

        By Mathieu's theorem these are, in every characteristic, exactly the
        constituent highest weights of a good filtration of the tensor product,
        so only the support of the product is wanted -- never the multiplicities,
        and never the characters themselves.

        Computed by the Littlewood--Richardson rule on Schur functions.  Two
        reductions make this cheap.  A determinant twist factors out of the
        tensor product: writing ``mu = mu_0 + a.det`` and ``nu = nu_0 + b.det``
        with ``mu_0, nu_0`` partitions, the constituents are those of
        ``mu_0 (x) nu_0`` shifted by ``(a+b).det``, so the cache is keyed on the
        untwisted pair and hits across every twist.  And ``s_eta`` vanishes in
        ``g`` variables once ``eta`` has more than ``g`` parts, which is exactly
        the truncation from ``GL_infinity`` down to ``GL_g``.
        """
        a, b = mu[-1], nu[-1]
        key = (tuple(x - a for x in mu), tuple(x - b for x in nu))
        hit = self._LR_cache.get(key)
        if hit is None:
            from sage.all import Partition
            prod = (self._schur(Partition([x for x in key[0] if x > 0]))
                    * self._schur(Partition([x for x in key[1] if x > 0])))
            hit = tuple(
                tuple(Integer(x) for x in list(eta) + [0] * (self.g - len(eta)))
                for eta in prod.support() if len(eta) <= self.g
            )
            self._LR_cache[key] = hit
        shift = a + b
        if shift == 0:
            return hit
        return tuple(tuple(x + shift for x in eta) for eta in hit)

    # ------------------------------------------------------------------
    def buildWeights(self, kmin, kmax):
        """Dominant (decreasing) weights with all entries in ``[kmin, kmax]``."""
        key = (kmin, kmax)
        hit = self._weights_cache.get(key)
        if hit is not None:
            return hit

        # A decreasing tuple is a weakly increasing one read backwards, so
        # enumerate weakly increasing sequences and reverse.
        res = []

        def rec(i, lo, acc):
            if i == self.g:
                res.append(tuple(Integer(x) for x in reversed(acc)))
                return
            for v in range(lo, kmax + 1):
                acc.append(v)
                rec(i + 1, v, acc)
                acc.pop()

        rec(0, kmin, [])
        self._weights_cache[key] = res
        return res

    def subsets_of_size(self, s, n):
        """Subsets of ``s`` of cardinality exactly ``n``."""
        from itertools import combinations
        return [list(c) for c in combinations(s, n)]

    # ------------------------------------------------------------------
    # The good filtration of Lambda^j (u_0^-)^v
    # ------------------------------------------------------------------
    def block_ranges(self, T):
        """The blocks of ``{0,...,g-1}`` cut out by ``I_0``."""
        out = []
        start = 0
        for j in range(self.g - 1):
            if T[j] == 0:
                out.append(list(range(start, j + 1)))
                start = j + 1
        out.append(list(range(start, self.g)))
        return out

    def _gl_weight_multiplicities(self, n, lam):
        """Weight multiplicities of the ``GL_n`` irreducible of highest weight ``lam``.

        ``lam`` is a decreasing ``n``-tuple.  Shifting it to a partition and
        expanding the Schur function in the monomial basis gives the Kostka
        numbers, which are exactly these multiplicities.
        """
        key = (n, tuple(lam))
        hit = self._gl_cache.get(key)
        if hit is not None:
            return hit
        from sage.all import Partition, SymmetricFunctions, ZZ
        from itertools import permutations
        c = lam[-1]
        lam0 = [x - c for x in lam]
        m = SymmetricFunctions(ZZ).monomial()
        expansion = m(self._schur(Partition([x for x in lam0 if x > 0])))
        out = {}
        for part, coeff in expansion.monomial_coefficients().items():
            base = list(part) + [0] * (n - len(part))
            if len(base) != n:
                continue
            for perm in set(permutations(base)):
                out[tuple(Integer(x + c) for x in perm)] = Integer(coeff)
        self._gl_cache[key] = out
        return out

    def nabla_L0_character(self, T, nu):
        """Weight multiplicities of ``nabla_{L_0}(nu)``, as a dict.

        ``L_0`` is a product of general linear groups, one per block, acting on
        disjoint sets of coordinates, so the character is the product of the
        block characters.
        """
        blocks = self.block_ranges(T)
        combos = {tuple([0] * self.g): Integer(1)}
        for blk in blocks:
            lam = [nu[i] for i in blk]
            per_block = self._gl_weight_multiplicities(len(blk), lam)
            new = {}
            for base, c0 in combos.items():
                for w, c1 in per_block.items():
                    cur = list(base)
                    for pos, i in enumerate(blk):
                        cur[i] = w[pos]
                    k = tuple(cur)
                    new[k] = new.get(k, 0) + c0 * c1
            combos = new
        return combos

    def good_filtration_hw(self, T, j):
        """``W_{I_0}(j)``: highest weights of the good filtration of ``Lambda^j(u_0^-)^v``.

        The character is ``sum_{|M| = j} e^{-s_M}`` over subsets ``M`` of ``S``.
        Constituents are peeled off greedily from the top: a maximal
        ``L_0``-dominant weight still present must be a highest weight, and its
        whole ``nabla_{L_0}`` character is then removed.  The loop terminates with
        an empty character, which is checked -- a non-empty remainder would mean
        the module had no good filtration.

        This is finer than "every ``L_0``-dominant ``-s_M``": when
        ``<-s_M, alpha^v> >= 2`` the weight ``-s_M - alpha`` is also
        ``L_0``-dominant but already lies inside ``nabla_{L_0}(-s_M)``, and is not
        a constituent of its own.  The two agree for ``g <= 3``.
        """
        key = (tuple(T), j)
        hit = self._gf_cache.get(key)
        if hit is not None:
            return hit

        phi_I, S = self.buildRoots(T)
        zero = tuple(Integer(0) for _ in range(self.g))
        char = {}
        for M in self.subsets_of_size(S, j):
            s_M = zero
            for x in M:
                s_M = add(s_M, x)
            w = tuple(-c for c in s_M)
            char[w] = char.get(w, 0) + 1

        out = []
        while char:
            dominant = [w for w in char if self.L0dominant(phi_I, w)]
            if not dominant:
                raise ArithmeticError(
                    "no L_0-dominant weight left in the character of "
                    "Lambda^%d(u_0^-)^v at I_0 = %s" % (j, tuple(T)))
            # Maximise a functional strictly positive on the positive roots of
            # L_0, so the chosen weight is maximal for the dominance order.
            nu = max(dominant,
                     key=lambda w: sum(-(i + 1) * w[i] for i in range(self.g)))
            mult = char[nu]
            out.extend([nu] * int(mult))
            for w, c in self.nabla_L0_character(T, nu).items():
                left = char.get(w, 0) - mult * c
                if left:
                    char[w] = left
                else:
                    char.pop(w, None)
                if left < 0:
                    raise ArithmeticError(
                        "negative multiplicity while decomposing "
                        "Lambda^%d(u_0^-)^v at I_0 = %s" % (j, tuple(T)))

        self._gf_cache[key] = out
        return out
