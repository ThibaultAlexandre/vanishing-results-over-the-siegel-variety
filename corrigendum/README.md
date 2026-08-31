# Reproducing the corrigendum's numerical results

One script, `reproduce.py`, recomputes every number Section 5 of the corrigendum
prints: the count tables of §5.2 and the lost weight of Proposition C.18. Each
figure is compared against the value the corrigendum carries, and the run ends in
`ALL FIGURES REPRODUCED` or lists what differs.

```
cd corrigendum
PYTHONPATH=.. sage -python reproduce.py          # every case, about two minutes
PYTHONPATH=.. sage -python reproduce.py 2 3 5    # seconds; skips the genus 4 run
```

Nothing is written to disk, and nothing here is imported by the `siegel` package
or by `main.ipynb`.

## What is compared

| | operator | where it lives |
|---|---|---|
| `Thm. 6.16` | Definition 6.15 **as printed**, with the coding defect of §2.4 repaired | `definition615.py`, this folder |
| `Thm. 6.16'` | the corrected statement | `siegel/theorem616p.py`, the package |

> **`definition615.py` is not part of the deposit.** The statement it implements
> does not hold as printed — that is what the corrigendum is about — and no
> vanishing result may be derived from it. It is here so that the comparison can
> be recomputed from source instead of taken on trust, and it is deliberately
> kept out of the `siegel` package. The published paper's *own* implementation is
> not reproduced here in any form; what is compared is the printed **statement**,
> with the coding slip (`M = ∅` skipped) already fixed, so
> that any difference in output is a difference of mathematics.

Both operators are driven by the same fixpoint loop, in `reproduce.py`, so a
difference in the output is a difference between the operators and not between
their drivers.

## Reading the counts

Every count comes from saturating an operator over a finite box `[-W,0]^g`. The
operators ask whether certain weights already lie in `C^{e+1}_van`; a weight
outside the box never does, so a candidate near the face is rejected for a reason
that has nothing to do with the statement being tested. **Every count is a lower
bound.**

Worse, the bias differs between the two operators: Definition 6.15 quantifies
over all `binom(d,n)` translated `T`-weights, Definition 6.15' only over the
Littlewood–Richardson constituents of the good filtrations, so the two consult
different weights and reach different depths. A difference at a single window is
therefore not yet a difference of operators.

The script controls for this. Each operator is saturated at two windows
`W1 < W2`, and the report is restricted to the largest sub-box `[-B,0]^g` on
which its two runs agree — `B` is measured, not assumed, and is checked against
the corrigendum's value. The windows:

| g | p | windows | B |
|---|---|---|---|
| 2 | 7 | `[-60,0]`, `[-80,0]` | 59 |
| 3 | 11 | `[-26,0]`, `[-34,0]` | 23 |
| 4 | 31 | `[-18,0]`, `[-24,0]` | 13 |

## The genus 5 witness

Proposition C.18 exhibits `λ = (-7,-11,-13,-16,-18)` at `p = 29`, `e = 3`: a
weight the published statement admits and the corrected one loses. The script
checks each assertion separately — the two constituents of `Λ¹² Sym² std₅`, which
of them makes `χ = λ - μ` dominant, that `λ' = χ + 2ρ` is `Z₀`-ample and orbitally
`p`-close, that `δ_{∅,3} = 2` while `⟨χ, α₁^∨⟩ = 1`, and that no parabolic and no
choice of base point admits `λ` under Theorem 6.16'.

This needs no saturation, which is the point: the witness lies at depth 18, out
of reach of any `g = 5` fixpoint the tables can afford, and the loss it exhibits
at `e = 3` is invisible in a saturation run at that genus.
