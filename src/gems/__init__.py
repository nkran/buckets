"""gems — mechanistic flux-generation module for CHO iFB soft-sensor features.

Pipeline: DESeq2-normalized counts + bioprocess metadata are integrated onto the
iCHO3K CHO-S genome-scale model via GIMME, constrained by measured exchange
(q-)rates, and solved with FBA/pFBA to yield an intracellular flux matrix.
"""

__version__ = "0.1.0"
