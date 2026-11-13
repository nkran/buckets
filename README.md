# gems — CHO iFB mechanistic flux-generation module

`gems` turns clone-level **transcriptomics + bioprocess data** into **predicted
intracellular metabolic fluxes** by integrating them onto the **iCHO3K (CHO-S)**
genome-scale metabolic model with **GIMME** and Flux Balance Analysis (FBA).

The resulting flux matrix is the **mechanistic feature layer** for a downstream
"Soft Sensor" that predicts whether a CHO clone will thrive under Intensified
Fed-Batch (iFB) — i.e. classify clones as **Responder** vs **Non-Responder**.

> **Scope of this module: fluxes only.**
> It *consumes* normalized counts + bioprocess metadata and *emits* a flux table.
> The ML classifier itself is a separate, downstream project and is **not** built here.
> q-rate (specific exchange rate) calculation *is* included, because it is what
> constrains the model's exchange reactions.

---

## 1. Big picture — how it works

```
 DESeq2 normalized counts (Ensembl PICRH gene ids)   bioprocess metadata
 + sample metadata                                    (VCD, Titer, Gluc, Lac, …)
        │                                                     │
        ▼                                                     ▼
 ┌─────────────┐                                      ┌───────────────┐
 │ gene_mapping│  Ensembl PICRH → Entrez → reaction   │  qrate_calc   │  ΔS / IVCD
 │             │  scores via GPR (AND=min, OR=max)    │               │  → exchange
 └──────┬──────┘                                      └──────┬────────┘    bounds
        │              ┌───────────────┐                     │
        │              │  model_setup  │  iCHO3K JSON        │
        │              │  base_cho_s   │  → CHO-S biomass    │
        │              └──────┬────────┘  + defined medium   │
        └──────────────┬──────┴─────────────────────────────┘
                       ▼
                ┌──────────────┐   GIMME as ONE weighted-parsimonious LP:
                │  gimme_fba   │   minimize Σ(1 + penalty_i)·|v_i|
                │              │   s.t. steady state, exchange bounds,
                └──────┬───────┘        biomass ≥ fraction · biomass_max (RMF)
                       ▼
                ┌──────────────┐   flux_table.csv / .parquet  (long format)
                │   export     │   manifest.json (coverage, params, infeasibles)
                └──────┬───────┘
                       ▼
                ┌──────────────┐   subsystem heatmap · key-reaction trajectories
                │     viz      │   differential flux: iFB-vs-FB, Responder-vs-NonResp
                └──────────────┘   (volcano · top-divergent · subsystem heatmap)
```

**GIMME in one line:** reactions whose genes are lowly expressed get a penalty;
the LP finds the lowest-total-flux steady state that still hits a required growth
level while avoiding those penalized reactions. Measured q-rates fix the exchange
fluxes. The output flux distribution is therefore both *expression-consistent* and
*measurement-consistent*.

---

## 2. Repository map — what is where

```
gems/
├── README.md                     ← this file
├── pyproject.toml                ← package metadata, deps, pytest config
│
├── src/gems/                     ← the Python package (installed as `gems`)
│   ├── config.py                 ← ALL knobs: paths, model ids, ID maps, medium,
│   │                               GIMME params, flux panel, experiment design
│   ├── data_io.py                ← load real data / validate schema; and the
│   │                               synthetic dataset generator (+ crosswalk)
│   ├── model_setup.py            ← load iCHO3K, build & cache the CHO-S base model
│   ├── qrate_calc.py             ← specific exchange q-rates → per-sample bounds
│   ├── gene_mapping.py           ← Ensembl→Entrez crosswalk + GPR reaction scoring
│   ├── gimme_fba.py              ← native GIMME LP + flux-panel extraction
│   ├── export.py                 ← assemble flux table + build/write manifest
│   ├── viz.py                    ← absolute + differential-flux figures
│   └── pipeline.py               ← orchestrates all stages; CLI entry point
│
├── tests/
│   ├── test_qrate_calc.py        ← q-rate sign/magnitude + bounds bracket estimate
│   ├── test_gene_mapping.py      ← GPR AND=min/OR=max + Ensembl→Entrez collapse
│   └── test_pipeline_synthetic.py← end-to-end smoke test (marked `slow`)
│
├── notebooks/
│   └── 01_explore_fluxes.ipynb   ← QC: distributions, FB-vs-iFB, Responder deltas
│
├── data/
│   ├── models/
│   │   ├── iCHO3K_cho_generic_unblocked.json   ← downloaded iCHO3K (source model)
│   │   └── base_cho_s.json                     ← cached curated CHO-S base model
│   ├── maps/
│   │   └── ensembl_picrh_to_entrez.tsv         ← gene crosswalk (synthetic run auto-writes)
│   ├── synthetic/                              ← generated placeholder counts + metadata
│   └── raw/                                    ← ⬅ put REAL counts.csv + metadata.csv here
│
└── results/
    ├── fluxes/
    │   ├── flux_table.csv / .parquet           ← THE DELIVERABLE (long-format fluxes)
    │   └── manifest.json                       ← provenance: model, solver, coverage, params
    └── figures/                                ← all PNGs + differential-flux CSVs
```

### Module responsibilities in detail

| Module | Key functions | Input → Output |
|---|---|---|
| `config.py` | (constants only) | Single source of truth. Everything tunable lives here; real-data items marked `# TODO(real-data)`. |
| `model_setup.py` | `get_base_model`, `build_base_model`, `apply_defined_medium` | iCHO3K JSON → CHO-S base model (objective `biomass_cho_s`, defined medium, verified growth), cached to `base_cho_s.json`. |
| `data_io.py` | `load_dataset`, `validate`, `make_synthetic_dataset` | `counts.csv`+`metadata.csv` → validated `Dataset`; or generates a schema-correct synthetic one. |
| `qrate_calc.py` | `compute_qrates`, `qrate_bounds_for_sample` | Bioprocess metadata → per-interval specific rates (via IVCD) → `{exchange_rxn: (lb, ub)}`. |
| `gene_mapping.py` | `counts_to_entrez`, `reaction_scores`, `map_expression` | Ensembl-indexed counts + crosswalk → per-reaction expression scores + coverage report. |
| `gimme_fba.py` | `run_gimme`, `extract_panel` | Reaction scores + exchange bounds → full flux vector → flux panel (curated subsystems). |
| `export.py` | `assemble_flux_table`, `build_manifest`, `write_outputs` | Per-sample panels → tidy long table + manifest, written to `results/fluxes/`. |
| `viz.py` | `make_all_figures`, `differential_flux`, `volcano`, … | Flux table → PNGs + differential CSVs in `results/figures/`. |
| `pipeline.py` | `run`, `main` | Wires everything; loops samples; CLI (`--synthetic` / `--raw` / `--limit`). |

---

## 3. The model (iCHO3K)

- Source: [LewisLabUCSD/iCHO3K](https://github.com/LewisLabUCSD/iCHO3K),
  file `iCHO3K/Model/iCHO3K_cho_generic_unblocked.json` (8,368 reactions /
  4,695 metabolites / 2,929 genes — the "unblocked" build).
- Genes are **NCBI Entrez** ids. Three biomass reactions ship in the model;
  we use **`biomass_cho_s`** (CHO-S WT), our sampled cell line.
- **Fetch note:** the raw GitHub CDN rate-limits (HTTP 429); the model was pulled
  via the GitHub **git-blob API** (base64). If you need to re-fetch, use that route.
- **License:** the iCHO3K repo did not show an explicit license — **confirm before
  redistributing** `data/models/*.json`.

---

## 4. Setup

```bash
cd /Users/nkranjc/work/gems
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"          # cobra, pandas, numpy, scipy, matplotlib, seaborn, pyarrow, pytest
```

Solver is open-source **GLPK** (via optlang). HiGHS would be faster but its optlang
interface isn't available on the installed optlang 1.9.1 / Python 3.9.

---

## 5. Run

```bash
# Full synthetic run (auto-generates placeholder data + crosswalk).
python -m gems.pipeline --synthetic

# Quick subset while iterating (recommended for a first look).
python -m gems.pipeline --synthetic --limit 24

# One project (144 samples).
python -m gems.pipeline --synthetic --limit 144

# Rebuild the cached base model (after changing the medium in config.py).
python -m gems.pipeline --synthetic --rebuild-model

# Real data (see §7).
python -m gems.pipeline --raw data/raw
```

**Where output goes.** Real-data runs write `results/fluxes/` and `results/figures/`
(what the explorer and notebooks read). `--synthetic` runs write to `results/synthetic/`
and keep their placeholder crosswalk in `data/synthetic/`, so they never overwrite the
real results or the BioMart crosswalk; tests write to a temp directory. Override with `--out DIR`.
The committed `results/` are the GSE140669 run (§7.1). `data/` and `escher_maps/` are not
in git — rebuild them with `scripts/build_gene_id_master.py`, `scripts/prepare_gse140669.py`
and `scripts/prepare_escher_maps.py`.

**Runtime:** about 11 s per sample on GLPK. Samples run in parallel on all CPU cores by
default (`--jobs 0`; `--jobs 1` for serial), so on a 10-core laptop the full
5×12×2×6 = 720-sample design takes roughly 15 minutes. Results are identical whatever
the number of workers: the solver is reset before every sample, because the LP has many
equally optimal flux distributions and would otherwise return one that depends on the
previously solved sample.

---

## 6. Outputs

### `results/fluxes/flux_table.csv` (and `.parquet`) — the deliverable
Long format, one row per (sample × reaction):

| column | meaning |
|---|---|
| `project, clone, condition, timepoint` | sample keys |
| `responder` | Responder / NonResponder label (if present in metadata) |
| `reaction_id` | iCHO3K reaction id |
| `subsystem` | pathway the reaction belongs to |
| `flux` | predicted flux (mmol/gDW/h scale; sign = direction) |

Only reactions in the **flux panel** are exported (10 central-carbon / lipid
subsystems + biomass + the 5 measured exchanges — see `config.PANEL_SUBSYSTEMS`).
This is the matrix you pivot to `samples × reactions` to feed the ML soft sensor.

### `results/fluxes/manifest.json` — provenance
Model id/size, solver, sample count, **mapping coverage** (Ensembl→Entrez ≈ 89%,
model-gene ≈ 85%), infeasible-sample log, and all GIMME/q-rate parameters used.

### `results/figures/`
- `flux_subsystem_heatmap.png` — mean |flux| per subsystem × condition/day.
- `flux_trajectories.png` — key exchange/biomass fluxes vs day, FB vs iFB.
- `diff_volcano_*`, `diff_top_*`, `diff_subsystem_*` — differential flux for
  **iFB vs FB** and **Responder vs NonResponder**, plus the underlying `diff_*.csv`.

---

## 7. Swapping in real data

1. Put two files in `data/raw/`:
   - **`counts.csv`** — DESeq2 **normalized** counts. Rows = Ensembl PICRH gene ids
     (`ENSCGRG…`, `Cricetulus_griseus_picr` assembly). Columns = sample ids.
   - **`metadata.csv`** — one row per sample id, columns:
     `project, clone, condition, timepoint` + bioprocess measurements
     `VCD, Viability, Titer, Gluc, Gln, Glu, Lac, NH3`; optional `responder`.
     (VCD in 1e6 cells/mL, metabolites in mM, Titer in mg/L.)
2. The **Ensembl PICRH → Entrez crosswalk** is `data/maps/ensembl_picrh_to_entrez.tsv`
   (columns `ensembl_gene_id, entrez_id, source`), the default the pipeline uses. It is
   built from Ensembl BioMart (dataset `cgpicr_gene_ensembl`, `CriGri-PICRH-1.0`,
   Ensembl 116): 21,909 direct BioMart links plus 83 genes with no BioMart Entrez link
   mapped through an unambiguous gene symbol (`source = via_symbol`). It reaches
   **2,800 of 2,929 iCHO3K model genes (95.6%)**. Versioned ids (`ENSCGRG00015007592.1`)
   are accepted; the version suffix is stripped before mapping. See §8 step 2.
3. Run `python -m gems.pipeline --raw data/raw`. For a dataset keyed by a different
   namespace (e.g. GSE140669's NCBI symbol/LOC), pass its own crosswalk with
   `--crosswalk data/maps/<file>.tsv`.

The `mapping_coverage` block in `results/fluxes/manifest.json` reports the achieved
gene coverage; `notebooks/03_map_preparation_report.ipynb` is a transparency report on
the crosswalk + map construction with worked examples.

The schema is enforced by `data_io.validate`, so mismatches fail fast with a clear error.

> **Note — bioprocess data is optional.** If `metadata.csv` lacks the bioprocess
> columns (a transcriptomics-only dataset), the pipeline detects this and runs GIMME
> **expression-only** (no exchange bounds). The manifest records
> `exchange_constraints_applied: false`.

### 7.1 Worked example — a real public dataset (GEO GSE140669)

A complete real-data run is wired up as an example. Dataset:
[GSE140669](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140669) — *"Defining
lncRNAs Associated with CHO Cell Growth and IgG Productivity by RNA-Seq"*: **4 CHO
clones × 3 timepoints (ST/D4/D7) × 3 replicates**, raw gene counts.

```bash
# 1. download the counts matrix (already fetched to data/raw/GSE140669_counts.csv)
#    from https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140669/suppl/
# 2. build pipeline inputs (collapses triplicates, CPM-normalizes, builds crosswalk)
python scripts/prepare_gse140669.py
# 3. run — expression-only, with late-vs-early fed-batch contrasts
python - <<'PY'
from gems import pipeline
pipeline.run(synthetic=False, raw_dir='data/raw',
             crosswalk_path='data/maps/gse140669_symbol_to_entrez.tsv',
             contrasts=[('timepoint',7,4),('timepoint',7,0)])
PY
```

Gene mapping uses **iCHO3K's own `gene_id_mapping.csv`** (symbol/LOC → Entrez),
downloaded to `data/maps/icho3k_gene_id_mapping.csv` — this is the real-data
equivalent of the BioMart crosswalk and gives **~66% model-gene coverage** here.

Result: 12 samples, 0 infeasible, 22k flux rows. The **D7-vs-D4** contrast shows a
biologically sensible late-culture shift — **oxidative phosphorylation and TCA flux
down, glycolysis / pentose-phosphate up** (see `results/figures/diff_subsystem_7_vs_4.png`).
Per-reaction p-values are weak (only 4 samples/timepoint, real clone heterogeneity),
so read the **subsystem-level** trends, not individual reactions.

`scripts/prepare_gse140669.py` is a good template for adapting any GEO counts matrix:
parse sample names → schema, collapse replicates, normalize, and point the crosswalk
at whatever gene-id namespace the dataset uses.

### 7.0 Start here: the ZeLa fed-batch tutorial

The best way to learn what the model does is `notebooks/04_zela_fedbatch_tutorial.ipynb`.
It runs a real fed-batch dataset end to end: wild-type CHO vs ZeLa (a line engineered to
secrete no lactate), 8 bioreactors, 43 samples over days 2-14, with **both** expression
and measured uptake/secretion rates (published with iCHO3K).

```bash
python scripts/prepare_zela.py          # downloads + formats the data into data/zela/
python -m gems.pipeline --raw data/zela --crosswalk data/zela/crosswalk.tsv --out results/zela
```

What it shows: 92% of model genes have data, all 43 samples solve, the growth the measured
uptake can support tracks the measured growth rate (correlation 0.74, about 1.5x too high),
and lactate dehydrogenase flux is switched off in ZeLa.

**Measured rates as input.** If a data folder contains `exchange_rates.csv` (samples x
exchange reactions, mmol/gDW/h, negative = uptake) the pipeline uses those instead of
computing rates from concentrations. Because raw measurements are rarely exactly consistent
with the network, `qrate_calc.reconcile_rates` first finds the closest feasible set of rates
(optionally growing at least as fast as `measured_growth.csv`) and writes every adjustment
to `exchange_reconciliation.csv`. Per-sample growth capacity and totals go to
`sample_summary.csv`. In the explorer app, pick the run with the **Dataset** selector.

### 7.2 Master gene-ID bridge (one table linking every namespace)

CHO data arrives keyed by many different identifiers. `scripts/build_gene_id_master.py`
builds **`data/maps/gene_id_master.tsv`** — one row per Chinese-hamster gene
(CriGri-PICRH-1.0, ~34.8k genes) linking them all:

```
ensembl_gene_id · entrez_id · symbol · refseq_mrna · refseq_protein ·
ensembl_protein · uniprot · gene_biotype · chromosome · in_icho3k · icho3k_n_reactions
```

Built from Ensembl BioMart (`cgpicr_gene_ensembl`) joined to the iCHO3K model.
**All 2,929 model genes are present.** How many model genes each namespace can reach:

| Source namespace | → iCHO3K model-gene coverage |
|---|---|
| `entrez_id` | 100% |
| `symbol` | 97.1% |
| `ensembl_gene_id` (PICRH) | 95.4% |
| `refseq_mrna` / `refseq_protein` | 95.2% |
| `uniprot` | 74.5% |

`gems.gene_ids` is the accessor:

```python
from gems import gene_ids as gid
gid.load_master()                                  # the whole table
gid.translate(['Has2'], 'symbol', 'refseq_mrna')   # cross-translate any namespaces
gid.detect_namespace(counts.index)                 # guess how a matrix is keyed
gid.make_best_crosswalk(counts.index, 'x.tsv')     # auto-build the entrez crosswalk
```

This bridge also **improved the GSE140669 run from 65.8% → 88.7%** model-gene coverage
(it's symbol-keyed, and the master's symbol→Entrez mapping is far more complete than
iCHO3K's shipped map). See `notebooks/03_map_preparation_report.ipynb` §2b.

---

## 8. Adapting the pipeline to your data (checklist)

Work through these in order. Everything you are expected to edit is in two files:
**`src/gems/config.py`** and **`src/gems/qrate_calc.py`**. Lines tagged
`# TODO(real-data)` are the ones still holding placeholders.

### Step 1 — Put your data in one folder (e.g. `data/raw/`)

| File | Required | Content |
|---|---|---|
| `counts.csv` | yes | DESeq2-normalised counts. Rows = Ensembl PICRH gene ids (`ENSCGRG…`, versions OK), columns = sample ids |
| `metadata.csv` | yes | index = sample id; columns `project, clone, condition, timepoint` (timepoint in **days**) |
| bioprocess columns in `metadata.csv` | optional | `VCD` (1e6 cells/mL), `Gluc, Gln, Glu, Lac, NH3` (mM), `Titer` (mg/L), `Viability` (%) — rates are then computed from consecutive timepoints |
| `exchange_rates.csv` | optional, preferred if you have it | samples × exchange reaction ids (`EX_glc_e`, `EX_lac_L_e`, …), **mmol/gDW/h**, negative = uptake. Overrides the concentration route |
| `measured_growth.csv` | optional | sample id, specific growth rate (1/h); used as a growth floor and for model-vs-measured comparison |

`condition` labels are free text; if there are exactly two (e.g. `FB`, `iFB`) the
differential-flux figures contrast them automatically. `responder` (optional column)
adds a Responder-vs-NonResponder contrast.

### Step 2 — Check the gene mapping (Ensembl PICRH)

Nothing to edit: the default crosswalk `data/maps/ensembl_picrh_to_entrez.tsv` is
built for PICRH. After a run, read `mapping_coverage.model_gene_coverage` in
`results/fluxes/manifest.json` — on full PICRH data it should be about 0.95. A much
lower value means the ids are from another assembly (CHOK1GS or CriGri_1.0 use
different `ENSCGRG…` ids) or from a different Ensembl release; rebuild with
`python scripts/build_gene_id_master.py`. Quick check before a full run:

```python
import pandas as pd
from gems import gene_mapping, model_setup
m = model_setup.get_base_model()
c = pd.read_csv('data/raw/counts.csv', index_col=0)
x = gene_mapping.load_crosswalk()
print(gene_mapping.map_expression(m, c, x, c.columns[0])[1].as_dict())
```

### Step 3 — Medium composition  ⟵ **edit `src/gems/config.py`**

This is where the AMBR basal medium + feeds go. The base model only lets the cell take
up what is listed here; everything else is closed.

| Setting (config.py) | What to put there |
|---|---|
| `_AMINO_ACIDS`, `_VITAMINS_COFACTORS` | exchange ids of the components actually in your basal medium + feed. Remove anything not in the medium |
| `BASE_MEDIUM_UPTAKES` | the full uptake list (glucose + the two lists above + anything else, e.g. `EX_pyr_e`, `EX_lac_L_e` if fed) |
| `MEDIUM_UPTAKE` | maximum uptake per component, mmol/gDW/h. Default 10 is effectively unlimited, which is why growth is currently the same in every sample. Realistic CHO amino-acid uptakes are 0.005–0.3; set this per component from your spent-media data if you can |
| `MEDIUM_FREE` | water, protons, O2, phosphate, sulfate, bicarbonate — normally leave as is |

To find the exchange id for a medium component, search the model:
```python
from gems import model_setup
m = model_setup.load_raw_model()
[(r.id, r.reactants[0].name) for r in m.exchanges if 'serine' in r.reactants[0].name.lower()]
```
After editing, rebuild the cached base model (it must still grow):
`python -m gems.pipeline --raw data/raw --rebuild-model --limit 1`.

### Step 4 — Which measured metabolites constrain the model  ⟵ `config.py`

`METABOLITE_EXCHANGE` maps your `metadata.csv` column names to exchange reactions.
Add entries for anything else you measure (e.g. `"Asn": "EX_asn_L_e"`, `"Ala": "EX_ala_L_e"`)
and add the same names to `BIOPROCESS_COLS`. Every extra measured metabolite tightens
the model; amino acids matter most because they otherwise stay capped only by Step 3.

### Step 5 — Units of the computed rates  ⟵ **edit `src/gems/qrate_calc.py`**

Only relevant on the concentration route (not with `exchange_rates.csv`).
`q = ΔS / IVCD × FLUX_SCALE`. With VCD in 1e6 cells/mL, ΔS in mM and time in days,
`ΔS / IVCD` is mmol per 10⁹ cells per day. To get mmol/gDW/h:

`FLUX_SCALE = 1 / (cell dry weight per 10⁹ cells in g × 24)`

For a typical CHO dry weight of 350 pg/cell (0.35 g per 10⁹ cells) that is
**0.119**. The current `FLUX_SCALE = 30.0` is a placeholder from the synthetic scaffold
— replace it, ideally with your own measured dry weight. `QRATE_TOL` (default 0.2) is
the ± band allowed around each rate.

### Step 6 — Biomass reaction  ⟵ `config.py`

`BIOMASS_RXN = "biomass_cho_s"` (CHO-S, non-producing). iCHO3K also has
`biomass_cho_prod` (producing line) and `biomass_cho` (generic). For IgG-producing
clones `biomass_cho_prod` is the closer composition. Rebuild the base model after changing.

### Step 7 — GIMME parameters (optional)  ⟵ `config.py`

`GIMME_QUANTILE` (0.25): reactions with expression below this quantile are penalised.
`GIMME_BIOMASS_FRACTION` (0.5): growth must stay at least this fraction of the maximum.
Defaults are reasonable; revisit once you see your expression distribution.

### Step 8 — Run and check

```bash
python -m gems.pipeline --raw data/raw --limit 6      # quick check
python -m gems.pipeline --raw data/raw                # full run -> results/fluxes, results/figures
```
Check in `results/fluxes/`: `manifest.json` (gene coverage, infeasible samples),
`sample_summary.csv` (growth capacity vs measured growth per sample) and, if you used
`exchange_rates.csv`, `exchange_reconciliation.csv` (how far each measured rate had to
move to be feasible — large values point to a missing medium component or a
measurement problem). `notebooks/04_zela_fedbatch_tutorial.ipynb` shows how to read
all of these on a real dataset.

---

## 8b. Interactive Escher explorer (flux + expression on pathway maps)

Browse iCHO3K pathways as **connected metabolic cascades** and overlay, per
clone/timepoint, predicted **flux** or gene **expression** on the reaction arrows
(colour + thickness), or compare them side by side.

```bash
pip install -e ".[explorer]"                 # escher + streamlit
python scripts/prepare_escher_maps.py        # download curated maps + escher.min.js
python -m gems.pipeline --raw data/raw       # emits flux_matrix + rxn_expression
streamlit run app/streamlit_app.py           # open the explorer
```

**Maps (~55 pathways, assembled from multiple sources).**
`scripts/prepare_escher_maps.py` builds the catalogue and writes
`escher_maps/manifest.json` with a **reliability tier** for every map
(**5 High · 40 Medium · 10 Low**):

- **🟢 High (5 maps)** — curated **RECON1 human** Escher maps (hand-drawn, connected,
  BiGG ids matching iCHO3K directly): Glycolysis/TCA/PPP (88%), Carbohydrate (86%),
  Inositol/retinol (82%), Amino acid (79%), Tryptophan (73%).
- **🟡 Medium / 🔴 Low (50 maps)** — **KEGG-derived**, using **human (hsa)** KEGG
  pathways so only human-relevant reactions appear (iCHO3K is Recon3D-derived).
  Layout + connectivity come from KEGG KGML; reaction identities are bridged to iCHO3K
  through a **union crosswalk**: BiGG's direct KEGG links **plus a MetaNetX (MNX)
  bridge** that recovers reactions BiGG doesn't link to KEGG. Metabolites and reactions
  carry **human-readable names** (e.g. "Acetyl-CoA", "ATP-Citrate lyase"), and currency
  metabolites (ATP/H2O/NAD…) are drawn as small secondary nodes. Coverage is now high
  for most: Butanoate 100%, Pantothenate/CoA 94%, Pyrimidine 90%, Arginine 87%,
  Gly/Ser/Thr 79%, Tyrosine 77%, plus TCA, purine, fatty-acid, folate, glycans, …

**How to interpret reliability** (shown in the app per map):
- `coverage` = mapped / total pathway reactions. Reactions not in iCHO3K render grey.
  For the KEGG maps this is now the fraction of the *human* pathway present in iCHO3K.
- RECON1 maps have no KEGG round-trip and are the most reliable (High).
- `ambiguous_frac` = reactions where one KEGG reaction matched several iCHO3K reactions
  (mostly compartment variants); one was chosen.

Sources: RECON1 from escher.github.io; KEGG hsa KGML + compound names from rest.kegg.jp;
BiGG namespace + MetaNetX `reac_xref` from bigg.ucsd.edu / metanetx.org (all cached in
`data/maps/`). The `gems.escher_maps.generate_subsystem_map` model-native generator
remains as a last-resort fallback (disconnected reaction-cards).

**How it renders.** `gems.escher_view.build_escher_html` uses Escher's `Builder` to
embed the map + data, then inlines the `escher.min.js` bundle as a `data:` URL so the
result is fully **offline** (no cobra/solver needed to browse — data comes from the
exported matrices). `app/streamlit_app.py` is the platform; `notebooks/02_escher_explore.ipynb`
is a scripted companion.

Encoding note: Escher ties a reaction's colour and thickness to one data vector, so
flux and expression are shown as **separate layers** (toggle) plus a **Compare** view —
not both channels on one arrow at once.

## 9. Testing

```bash
pytest -q                 # fast unit tests (q-rates, GPR, crosswalk)
pytest -q -m slow         # + end-to-end synthetic pipeline (~35 s)
```

---

## 10. What to do next

**A. Immediate — make it run on real biology**
1. Drop the real DESeq2 `counts.csv` + `metadata.csv` into `data/raw/` (§7).
2. Build the real `ensembl_picrh_to_entrez.tsv` from BioMart and check the
   `mapping_coverage` in the manifest — aim for model-gene coverage > 0.8.
3. Encode the **actual AMBR basal medium + feed** in `config.BASE_MEDIUM_UPTAKES`,
   then `--rebuild-model` and confirm positive growth.

**B. Calibrate the exchange constraints**
4. Convert q-rates to true `mmol/gDW/h` (cell dry weight ≈ 300 pg/cell is a common
   CHO assumption) and set `FLUX_SCALE` accordingly; tighten the medium so glucose
   is the dominant carbon source. Re-run and confirm biomass now varies per clone.
5. Add a q-rate sanity report (carbon/nitrogen balance) before feeding bounds.

**C. Scale + robustness**
6. Run the full 720-sample design; parallelize the per-sample loop (`joblib` /
   `multiprocessing`) — each sample is independent, so this scales linearly.
7. Optionally upgrade optlang to get the **HiGHS** solver for a large speedup.

**D. Toward the soft sensor (separate module, not built here)**
8. Pivot `flux_table` to `samples × reactions`, join with **early-timepoint**
   (Day 3/5) bioprocess features, and label by final iFB outcome.
9. Train/evaluate the Responder classifier with **leave-one-project-out** CV (5
   projects) to prove the markers are *universal*, not product-specific.
10. Feed the ranked differential-flux reactions (e.g. TCA / glutamate / lipid
    hubs from `results/figures/diff_*.csv`) as candidate interpretable features.

**E. Optional enhancements**
11. Add **Escher** pathway maps colored by Δflux (dependency already declared under
    `[project.optional-dependencies] viz`).
12. Add flux **sampling** or **FVA** to quantify confidence intervals per reaction,
    not just point estimates.
