"""Central configuration: paths, model identifiers, ID maps, and the flux panel.

Everything here is intentionally data-driven and tunable. The values that must be
revised once the real AMBR/DESeq2 data lands are flagged with ``# TODO(real-data)``.
"""

from __future__ import annotations

from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
RESULTS = ROOT / "results"

MODEL_JSON = DATA / "models" / "iCHO3K_cho_generic_unblocked.json"
BASE_MODEL_JSON = DATA / "models" / "base_cho_s.json"        # curated cell-line model (cached)
CROSSWALK_TSV = DATA / "maps" / "ensembl_picrh_to_entrez.tsv"

SYNTH_DIR = DATA / "synthetic"
FLUX_DIR = RESULTS / "fluxes"
FIG_DIR = RESULTS / "figures"

# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------
# iCHO3K ships three biomass reactions; CHO-S is our sampled cell line.
BIOMASS_RXN = "biomass_cho_s"
CHO_S_GROWTH_TOL = 1e-3  # a "healthy" base model must exceed this on the defined medium

# Measured extracellular metabolites -> iCHO3K exchange reaction id.
# Sign convention (cobra): uptake is a NEGATIVE flux, secretion POSITIVE.
METABOLITE_EXCHANGE = {
    "Gluc": "EX_glc_e",     # glucose
    "Gln": "EX_gln_L_e",    # L-glutamine
    "Glu": "EX_glu_L_e",    # L-glutamate
    "Lac": "EX_lac_L_e",    # L-lactate
    "NH3": "EX_nh4_e",      # ammonium
}
# The product (IgG) has no metabolite in a generic GSMM; qP is reported but not
# imposed as an exchange constraint in this build.
PRODUCT_EXCHANGE = None

# ---------------------------------------------------------------------------
# Defined base medium (CHO-S)
# ---------------------------------------------------------------------------
# Exchange reactions opened for uptake when building the base cell-line model.
# Placeholder chemically-defined-medium proxy: glucose + glutamine as the primary
# C/N sources, all 20 amino acids, gases/salts, and common vitamin/cofactor
# precursors. Bounds are generous (mmol/gDW/h); real per-run bounds are imposed
# later per sample from q-rates.
# TODO(real-data): replace with the actual AMBR basal medium + feed composition.
_AMINO_ACIDS = [
    "EX_ala_L_e", "EX_arg_L_e", "EX_asn_L_e", "EX_asp_L_e", "EX_cys_L_e",
    "EX_gln_L_e", "EX_glu_L_e", "EX_gly_e", "EX_his_L_e", "EX_ile_L_e",
    "EX_leu_L_e", "EX_lys_L_e", "EX_met_L_e", "EX_phe_L_e", "EX_pro_L_e",
    "EX_ser_L_e", "EX_thr_L_e", "EX_trp_L_e", "EX_tyr_L_e", "EX_val_L_e",
]
_VITAMINS_COFACTORS = [
    "EX_chol_e", "EX_inost_e", "EX_pydx_e", "EX_pydxn_e", "EX_pydam_e",
    "EX_thm_e", "EX_ribflv_e", "EX_fol_e", "EX_ncam_e", "EX_nac_e",
    "EX_pnto_R_e", "EX_btn_e", "EX_ptrc_e", "EX_adn_e", "EX_ade_e",
    "EX_gua_e", "EX_ura_e", "EX_cytd_e", "EX_thymd_e", "EX_hxan_e",
    "EX_gthrd_e", "EX_crn_e",
]
# uptake bound (mmol/gDW/h) applied to defined-medium components
MEDIUM_UPTAKE = 10.0
# essentially-unlimited components (water, protons, gases, bulk ions)
MEDIUM_FREE = ["EX_h2o_e", "EX_h_e", "EX_o2_e", "EX_pi_e", "EX_so4_e", "EX_hco3_e"]

BASE_MEDIUM_UPTAKES = ["EX_glc_e"] + _AMINO_ACIDS + _VITAMINS_COFACTORS

# ---------------------------------------------------------------------------
# GIMME parameters
# ---------------------------------------------------------------------------
# Reactions whose expression score falls below GIMME_QUANTILE are penalised.
GIMME_QUANTILE = 0.25          # 25th percentile of reaction expression scores
# Required Metabolic Functionality: biomass held at this fraction of its optimum.
GIMME_BIOMASS_FRACTION = 0.5
# Genes/reactions with no expression data get this neutral score (not penalised).
UNMAPPED_SCORE = float("nan")

# ---------------------------------------------------------------------------
# Flux panel — subsystems whose reactions are reported in the flux matrix.
# ---------------------------------------------------------------------------
PANEL_SUBSYSTEMS = [
    "GLYCOLYSIS/GLUCONEOGENESIS",
    "CITRIC ACID CYCLE",
    "PENTOSE PHOSPHATE PATHWAY",
    "PYRUVATE METABOLISM",
    "OXIDATIVE PHOSPHORYLATION",
    "GLUTAMATE METABOLISM",
    "FATTY ACID METABOLISM",
    "GLYCEROLIPID METABOLISM",
    "GLYCEROPHOSPHOLIPID METABOLISM",
    "SQUALENE AND CHOLESTEROL SYNTHESIS",
]
# Always include these individually, regardless of subsystem.
PANEL_EXTRA_RXNS = [BIOMASS_RXN] + list(METABOLITE_EXCHANGE.values())

# ---------------------------------------------------------------------------
# Experimental design (used by the synthetic generator + validation)
# ---------------------------------------------------------------------------
N_PROJECTS = 5
N_CLONES = 12
CONDITIONS = ["FB", "iFB"]
TIMEPOINTS = [3, 5, 7, 9, 11, 13]   # culture days (6 longitudinal samples)

# Metadata / count-matrix schema
SAMPLE_KEYS = ["project", "clone", "condition", "timepoint"]
BIOPROCESS_COLS = ["VCD", "Viability", "Titer", "Gluc", "Gln", "Glu", "Lac", "NH3"]
# VCD in 1e6 cells/mL, metabolites in mM, Titer in mg/L.
