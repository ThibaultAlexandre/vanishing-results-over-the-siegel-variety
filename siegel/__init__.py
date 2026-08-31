"""Vanishing results for the coherent cohomology of automorphic vector bundles
over the Siegel variety in positive characteristic.

Implements the operator of Definition 6.15' and Theorem 6.16' of the corrigendum
to

    T. Alexandre, *Vanishing results for the coherent cohomology of automorphic
    vector bundles over the Siegel variety in positive characteristic*,
    Algebra & Number Theory 19 (2025), no. 1, 143-193; arXiv:2202.06691.

The primed statements are the *corrected* ones.  They are not in the published
paper: the degeneration apparatus feeding its Theorem 6.16 does not hold as
stated, and the repair adds the threshold ``delta_{I_0,e}`` to the dominance
condition.  Nothing here implements the published Theorem 6.16.

Typical use::

    from siegel import SiegelData, VanishingEngine

    data = SiegelData(g=2, p=7)
    X = VanishingEngine(data)
    X.saturate(-70, 0)
    X.statistics()
"""

from .rootdata import SiegelData
from .engine import VanishingEngine
from .theorem616p import Theorem616p

__all__ = ["SiegelData", "VanishingEngine", "Theorem616p"]
