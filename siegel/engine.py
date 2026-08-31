"""Iterating the operator of Theorem 6.16' to a fixpoint.

``Cvan[k]`` is the set of weights ``lambda`` for which
``H^i(X, nabla(lambda)(-D)) = 0`` is known for every ``i > k`` -- equivalently,
cohomology concentrated in degrees ``[0:k]``.  The sets increase with ``k``.

Applying Theorem 6.16' at one ``(I_0, e)`` can enlarge ``Cvan[e]``, which may in
turn let another ``(I_0, e)`` succeed, so the operator is swept over all
parabolics and degrees and the sweep repeated until nothing new appears.

The conditions of Definition 6.15' split into a part independent of what is
already known and a part that is a membership test.  :meth:`saturate` computes
the first once per ``(I_0, e)`` and then loops only over the second, so repeating
a sweep costs a subset test per surviving weight rather than a fresh evaluation
of the operator.
"""

from sage.all import Integer

from .theorem616p import Theorem616p


class VanishingEngine:
    """Vanishing results for ``X``, obtained by iterating Theorem 6.16'."""

    def __init__(self, data):
        self.data = data
        self.operator = Theorem616p(data)
        self.Cvan = [set() for _ in range(data.d)]
        self._compiled = {}

    # ------------------------------------------------------------------
    def vanishes(self, k, char):
        """Is ``H^i(X, nabla(char)(-D)) = 0`` known for every ``i > k``?"""
        return tuple(Integer(x) for x in char) in self.Cvan[k]

    def concentration(self, char):
        """Least ``k`` with cohomology known to be concentrated in ``[0:k]``.

        ``None`` if nothing is known for this weight.
        """
        w = tuple(Integer(x) for x in char)
        for k in range(self.data.d):
            if w in self.Cvan[k]:
                return k
        return None

    # ------------------------------------------------------------------
    def _compile(self, T, e, kmin, kmax, verbose):
        key = (tuple(T), e, kmin, kmax)
        if key not in self._compiled:
            weights = self.data.buildWeights(kmin, kmax)
            self._compiled[key] = self.operator.compile(T, e, weights)
            if verbose:
                print("  compiled I_0 = %s, e = %d: %d candidates of %d weights"
                      % (tuple(T), e, len(self._compiled[key]), len(weights)))
        return self._compiled[key]

    def sweep(self, kmin, kmax, verbose=False):
        """One pass over every parabolic and every degree.  ``True`` if it grew."""
        found = False
        for T in self.data.parabolics():
            for e in range(self.data.d):
                if self.apply(T, e, kmin, kmax, verbose=verbose):
                    found = True
        return found

    def apply(self, T, e, kmin, kmax, verbose=False):
        """Apply Theorem 6.16' at one parabolic and degree."""
        compiled = self._compile(T, e, kmin, kmax, verbose)
        d = self.data
        known = self.Cvan[e + 1] if e + 1 < d.d else None
        here = self.Cvan[e]
        gained = []
        for lam, reqs in compiled:
            if lam in here:
                continue
            if known is None:
                gained.append(lam)
                continue
            for req in reqs:
                if req <= known:
                    gained.append(lam)
                    break
        for lam in gained:
            for k in range(e, d.d):
                self.Cvan[k].add(lam)
        return bool(gained)

    def saturate(self, kmin, kmax, max_rounds=64, verbose=False):
        """Sweep until no new weight appears.  Returns the number of sweeps."""
        rounds = 0
        while rounds < max_rounds:
            rounds += 1
            if verbose:
                print("sweep %d" % rounds)
            if not self.sweep(kmin, kmax, verbose=verbose and rounds == 1):
                break
        return rounds

    # ------------------------------------------------------------------
    def convert(self):
        """Weights bucketed by the least degree at which each was obtained."""
        res = []
        seen = set()
        for c in self.Cvan:
            res.append([list(x) for x in sorted(c) if x not in seen])
            seen |= c
        return res

    def statistics(self, show=True):
        res = self.convert()
        counts = [len(x) for x in res]
        if show:
            total = 0
            for i in range(self.data.d):
                total += counts[i]
                print("H^* concentrated in degrees [0:%d] for %5d characters "
                      "(%5d new)" % (i, len(self.Cvan[i]), counts[i]))
        return counts

    # ------------------------------------------------------------------
    def _path(self, k, directory):
        return "%s/g%dp%d_%d.txt" % (directory, self.data.g, self.data.p, k)

    def save(self, directory="save", verbose=True):
        res = self.convert()
        for k in range(self.data.d):
            path = self._path(k, directory)
            with open(path, "w") as f:
                for item in res[k]:
                    f.write("".join("%s " % x for x in item) + "\n")
            if verbose:
                print("Results saved in " + path)

    def load(self, directory="save", verbose=True):
        for k in range(self.data.d):
            path = self._path(k, directory)
            chars = []
            with open(path, "r") as f:
                for line in f.read().split("\n"):
                    y = line.split(" ")
                    y.pop()
                    char = tuple(Integer(z) for z in y)
                    if len(char) == self.data.g:
                        chars.append(char)
            for i in range(k, self.data.d):
                self.Cvan[i].update(chars)
            if verbose:
                print("Results loaded from " + path)
