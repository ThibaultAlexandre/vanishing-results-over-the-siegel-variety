# Vanishing results over the Siegel variety

SageMath code computing vanishing results for the coherent cohomology of automorphic
vector bundles over the Siegel variety in positive characteristic. It accompanies

> T. Alexandre, *Vanishing results for the coherent cohomology of automorphic vector
> bundles over the Siegel variety in positive characteristic*, Algebra & Number Theory
> **19** (2025), no. 1, 143–193; [arXiv:2202.06691](https://arxiv.org/abs/2202.06691),

and its corrigendum.

**Which theorem is implemented.** The code iterates the corrected vanishing theorem:
Theorem 6.16′ of the corrigendum, which is Theorem 6.16 of arXiv:2202.06691v3. It does
not implement Theorem 6.16 as published. The algorithm of the published version is kept
unchanged at the tag `v1-as-published`.

## Requirements

SageMath and `matplotlib`.

## Quick start

```python
import sys; sys.path.insert(0, '.')
from siegel import SiegelData, VanishingEngine

data = SiegelData(g=2, p=7)      # the code requires p > g^2
X = VanishingEngine(data)
X.saturate(-70, 0)               # iterate to a fixpoint over [-70,0]^g
X.vanishes(1, (-4, -6))          # True if H^i = 0 is known for all i > 1
```

Weights are written in the labelling of the paper (its Definition 3.10); for instance
`Sym² std` has label `(0,...,0,-2)`. Ampleness is tested by a sufficient criterion, so
the computed sets are lower bounds for what the theorem gives.

`main.ipynb` is a guided tour, and writes the tables in `save/`.

## Contents

```
siegel/        the package
main.ipynb     guided tour
tests.py       self-checks:  PYTHONPATH=. sage -python tests.py
save/          computed results: save/g{g}p{p}_{k}.txt lists the weights first
               obtained in degree k
corrigendum/   reproduction of the corrigendum's computations (see its README)
```

## Licence

MIT. See `LICENSE`.
