# Vanishing results over the Siegel variety

Computes vanishing results for the coherent cohomology of automorphic vector bundles
over the Siegel variety in positive characteristic.

Let `p` be a prime, `N >= 3` an integer with `p` not dividing `N`, and let `X` be a
smooth projective toroidal compactification of the Siegel variety of dimension
`d = g(g+1)/2` over `F_p`, with boundary `D`. For a character `λ` of `GL_g`, write
`∇(λ)` for the costandard automorphic vector bundle of highest weight `λ`. This
code computes characters `λ` for which

```
H^i(X, ∇(λ)(-D)) = 0     for all i > k
```

by iterating the operator `g_{I_0,e}` of Definition 6.15', as licensed by
**Theorem 6.16'**, which assumes `p > g^2`.

Accompanies [arXiv:2202.06691](https://arxiv.org/abs/2202.06691) —
T. Alexandre, *Vanishing results for the coherent cohomology of automorphic vector
bundles over the Siegel variety in positive characteristic*, Algebra & Number
Theory **19** (2025), no. 1, 143–193.

> **Which statement is implemented.** The primed numbers — Definition 6.6',
> 6.10', 6.15' and Theorem 6.16' — are those of the **corrigendum**, not of the
> version published in Algebra & Number Theory. Three steps of Section 6 of that
> version do not hold as stated, and the repair adds a threshold `δ_{I_0,e}` to
> the dominance condition (D). This code implements the corrected operator only;
> it does not implement the Theorem 6.16 of the published version, and its output
> is therefore not that paper's. The threshold vanishes identically for `g ≤ 3`,
> so every `g = 2` and `g = 3` vanishing statement of the published version is
> unaffected; at `g ≥ 5` and `e ≥ 3` the corrected operator loses weights the
> published one admitted. Unprimed references below (Theorem 5.12, Remark 6.17)
> are to the paper and are untouched by the corrigendum.

**Where to read the corrected statements.** The corrigendum is the version of
record for the primed numbers used here, and delimits precisely what is affected.
The paper itself is at [arXiv:2202.06691](https://arxiv.org/abs/2202.06691). Paper,
corrigendum and code are all in one convention — that of the article, including
the labelling of its Definition 3.10 — so a weight may be carried between them
unchanged (see [Conventions](#conventions) below).

## Requirements

SageMath (developed against 10.5) and `matplotlib`. Root systems, the Weyl group of
`C_g`, the plethysm `Λ^n Sym² std` and the Littlewood–Richardson rule all come from
Sage.

## Quick start

```python
import sys; sys.path.insert(0, '.')
from siegel import SiegelData, VanishingEngine

data = SiegelData(g=2, p=7)      # enforces p > g^2
X = VanishingEngine(data)
X.saturate(-70, 0)               # iterate to a fixpoint over [-70,0]^g
X.statistics()

X.vanishes(1, (-4, -6))          # is H^i = 0 known for all i > 1?
X.concentration((-4, -6))        # least degree of concentration
X.save()
```

`main.ipynb` is a guided tour: the ingredients of Definition 6.15', certificates for
individual weights, which parabolics carry a computation, plots at `g = 2` and
`g = 3`, behaviour in `p`, a genus 4 run, and — §11 — the regeneration of the tables
behind Figures 5–9 of the paper.

## What is computed

For a standard parabolic `I_0` and a degree `e`, write `S` for the positive roots of
`L` outside `L_0`, `r_0 = |S|` and `2ρ_{I_0} = Σ S`. For a constituent `μ^{d-e}_top`
of Definition 6.6' put `χ = λ - μ^{d-e}_top` and `λ' = χ + 2ρ_{I_0}`. Then `λ` lies in
`g^{μ_top}_{I_0,e}(C)` exactly when

- **(D)** `χ ∈ X*(P_0)` and `⟨χ, α^∨⟩ ≥ δ_{I_0,e}(α)` for every simple `α` outside `I_0`;
- **(A)** `λ'` lies in the ample cone at `I_0`;
- **(C1)** for `1 ≤ k ≤ min(e, r_0)`, every `-s_M ∈ W_{I_0}(r_0-k)` with `λ' - s_M`
  dominant, and every `j`: `LR(μ^{d-e+k}_j, λ' - s_M) ⊆ C`;
- **(C2)** for every `j`: `LR(μ^{d-e}_j, χ) \ {λ} ⊆ C`,

and `g_{I_0,e}(C)` is the union over the **maximal** constituents `μ_top`. Theorem 6.16': if
`p > g^2` and `C ⊆ C^{e+1}_van`, then `g_{I_0,e}(C) ⊆ C^e_van`.

Two points where the implementation follows the statement closely rather than
loosely:

- `μ^n_j` are the **constituent highest weights** of `Λ^n Sym² std` (Definition 6.6'),
  not all of its `T`-weights. They are pairwise incomparable, so there is in general
  no single highest weight and the operator unions over the maximal ones — at `g = 3`,
  `n = 3` there are already two.
- `W_{I_0}(j)` is the multiset of highest weights of the good filtration of
  `Λ^j(u_0^-)^∨`, obtained by decomposing its character, not the larger set of all
  `L_0`-dominant `-s_M` of size `j`. The two agree for `g ≤ 3` but not beyond: if
  `⟨-s_M, α^∨⟩ ≥ 2` then `-s_M - α` is `L_0`-dominant yet already lies inside
  `∇_{L_0}(-s_M)`.

The ampleness test uses the sufficient subset of the `D`-ample cone furnished by
Theorem 5.12 and Remark 6.17 — orbitally `p`-close and `Z_0`-ample. The `D`-ample cone
itself is not reconstructible from the paper, so **every count produced here is a
lower bound** on the true output.

## Conventions

Everything here is in the convention of the paper, including the labelling of its
Definition 3.10: an automorphic bundle is labelled by the character one induces
from, so that `π_*L_λ ≃ ∇(λ)` is label-preserving, and the underlying module has
highest weight the image of the label under the twist

```
τ(a_1,...,a_g) = (-a_g,...,-a_1)
```

(reverse and negate — the action of `-w_0` for `w_0` the longest element of
`GL_g`; the same map is written `w_0 w_{0,L}` when `w_0` is taken in `Sp_2g`).
`SiegelData.changeConvention` is `τ`, and it is applied in exactly the two places
where the package meets representation theory stated intrinsically: the plethysm
constituents of `Λ^n Sym² std` in `_build_form_weights`, and the `p`-small weights.
Everything else — roots, `S`, `2ρ_{I_0}`, `χ`, the threshold, the ampleness test —
is a character of a parabolic or a weight in the root lattice, and carries no twist.

The anchor to check against: `Sym² std` has intrinsic highest weight `(2,0,...,0)`
and label `(0,...,0,-2)`, which is what `data.mu[0]` holds. At `g = 3, n = 3` the
two constituents are `(4,1,1), (3,3,0)` intrinsically and `(-1,-1,-4), (0,-3,-3)`
as labels, which is what `data.mu[2]` holds.

Both labellings share one notion of dominance, so "**dominant means decreasing**"
holds either way and does not identify which is in force: `(-2,-8)` is dominant,
`(-8,-2)` is not. The files in `save/` are in the labelling above, so their weights
are directly comparable with those printed in the paper and in the corrigendum, with
no twist to apply.

`Cvan[k]` is the set of weights whose cohomology is known to be concentrated in
degrees `[0:k]`; these sets increase with `k`.

## Performance

Conditions (D), (A) and the Littlewood–Richardson sets of (C1)–(C2) do not depend on
what is already known — only the final memberships do. `Theorem616p.compile` evaluates
the first group once per `(I_0, e)` and returns, for each weight, the set whose
membership would suffice; a sweep then costs one subset test per candidate, so
iterating to a fixpoint is cheap. Two further points matter:

- only the **support** of `∇(μ) ⊗ ∇(ν)` is ever needed, so it is computed by the
  Littlewood–Richardson rule on Schur functions rather than by building characters;
  a determinant twist factors out, and caching on the untwisted pair hits across every
  twist;
- `max_w |⟨χ, w(α^∨)⟩|` depends only on the Weyl orbit of `α^∨`, so the ampleness test
  scans orbits instead of pairs `(w, α)` — `2g²` pairings per weight rather than
  `|W| · |Φ⁺|`.

Genus 3 over `[-25,0]³` saturates in about a second, and genus 4 is within reach.

## Layout

```
siegel/rootdata.py     root systems, forms, thresholds, ampleness, LR, good filtrations
siegel/theorem616p.py   the operator of Definition 6.15'
siegel/engine.py       iteration to a fixpoint, save/load
siegel/analysis.py     certificates, per-parabolic contributions, threshold tables
siegel/plotting.py     scatter plots, slices, comparison across primes
siegel/psmall.py       p-small weights for Sp_2g
main.ipynb             guided tour
tests.py               self-checks
save/                  computed results
corrigendum/           reproduction of the corrigendum's Section 5 (see below)
```

## Tests

```
PYTHONPATH=. sage -python tests.py
```

Checks the Littlewood–Richardson rule against Sage's Weyl character ring, that each
good filtration reproduces the character it decomposes, that the threshold matches
Definition 6.10' at `g = 2,...,5`, and that a saturated computation is a fixpoint,
deterministic, and independent of the order in which parabolics are swept. It also
pins the labelling: that the twist is an involution fixing `2ρ` and preserving
dominance, and that `mu[0]` is the anchor `(0,...,0,-2)` — dominance alone cannot
tell the two labellings apart, so the anchor is what fixes which one is in force.

## Data

`save/g{g}p{p}_{k}.txt` — weights first obtained at degree `k`, one per line.
`save/g{g}p{p}_psmall.txt` — `p`-small weights.

One window per genus: `[-70,0]` at `g = 2` and `[-25,0]` at `g = 3`, so that the tables
of a given genus are comparable across primes and one table set per `(g, p)` serves every
figure of the paper.

| g | p | window | counts per degree |
|---|---|---|---|
| 2 | 5 | [-70,0] | 899, 1086, 1139 |
| 2 | 7 | [-70,0] | 1354, 1537, 1596 |
| 2 | 11 | [-70,0] | 1744, 1919, 1983 |
| 2 | 31 | [-70,0] | 2138, 2273, 2341 |
| 3 | 11 | [-25,0] | 986, 1184, 1212, 1226, 1269, 1709 |
| 3 | 691 | [-25,0] | 1351, 1627, 1655, 1687, 1752, 2231 |

Every one of these tables is written by a cell of `main.ipynb`, and `save/` holds nothing
else: run the notebook top to bottom on an empty `save/` and it comes back exactly as
listed. The output grows with `p` and then stops — at `g = 3` it has stabilised by
`p = 31`, which agrees with `p = 691` entry for entry over this window; §8 of the
notebook computes `p = 31` for that comparison without saving it.

### The tables behind the figures

Figures 5–9 of the paper are `pgfplots` scatter plots, one table per degree of
concentration. §11 of `main.ipynb` recomputes all of those tables from Theorem 6.16',
under the names the figures read:

| figure | g | p | tables |
|---|---|---|---|
| 5 (and 1) | 2 | 5 | `save/g2p5_{0,1,2}.txt`, `save/g2p5_psmall.txt` |
| 6 | 2 | 11 | `save/g2p11_{0,1,2}.txt`, `save/g2p11_psmall.txt` |
| 7 | 2 | 31 | `save/g2p31_{0,1,2}.txt`, `save/g2p31_psmall.txt` |
| 8 | 3 | 11 | `save/g3p11_{0,...,5}.txt` |
| 9 | 3 | 691 | `save/g3p691_{0,...,5}.txt` |

The window is the one the genus uses throughout, `[-70,0]` at `g = 2` and `[-25,0]` at
`g = 3` — the wider of the windows the published figures were drawn on, Figure 7 reaching
`-70` where Figures 5 and 6 stop at `-50`. Each figure clips to its own axis range.

The corrigendum reprints Figures 8 and 9, from these files. The genus 2 pictures are
correct as printed and are not reproduced there; their tables are regenerated here all
the same, so that nothing in §7 of the paper is beyond the deposit's reach.

## Reproducing the corrigendum's numerical results

Section 5 of the corrigendum compares the published Theorem 6.16 with Theorem 6.16',
degree by degree, and exhibits a weight the correction loses at `g = 5`. Both are
recomputed by one script:

```
cd corrigendum
PYTHONPATH=.. sage -python reproduce.py          # every case, about two minutes
PYTHONPATH=.. sage -python reproduce.py 2 3 5    # seconds; skips the genus 4 run
```

It prints the count tables of §5.2 and checks every assertion of Proposition C.18,
comparing each figure against the value the corrigendum carries; the run ends in
`ALL FIGURES REPRODUCED` or lists what differs. Nothing is written to disk.

What it compares is the printed **statement** against the corrected one. The
published Definition 6.15 is implemented in `corrigendum/definition615.py`, with the
coding defect of the corrigendum's §2.4 already repaired — the skipped subset
`M = ∅` — so that any difference in output is a
difference of mathematics rather than of code. The paper's own implementation is not
reproduced in any form.

> `corrigendum/definition615.py` is **not part of the deposit**. The statement it
> implements does not hold as printed, and no vanishing result may be derived from
> it; it is kept out of the `siegel` package and out of `main.ipynb`, and exists only
> so that the comparison can be checked from source. The deposited statement is
> Theorem 6.16', in `siegel/theorem616p.py`.

Both operators are driven by the same fixpoint loop, and each is saturated at two
windows and reported only on the sub-box where its two runs agree — every count over
a box being a lower bound, and the two operators being unequally sensitive to the
box. `corrigendum/README.md` explains the control and lists the windows.

## Licence

MIT. See `LICENSE`.
