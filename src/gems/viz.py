"""Flux and differential-flux visualization from the long-format flux table.

All figures are driven by the same tidy table produced by :mod:`export`.
Absolute:      subsystem heatmap, key-reaction trajectories.
Differential:  volcano, top-divergent barplot, subsystem-aggregated heatmap for a
               contrast (default iFB vs FB; Responder vs NonResponder when labelled).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import config as C

_EPS = 1e-9


def _sample_col(df: pd.DataFrame) -> pd.Series:
    return (df["project"].astype(str) + "|C" + df["clone"].astype(str)
            + "|" + df["condition"].astype(str) + "|D" + df["timepoint"].astype(str))


# ---------------------------------------------------------------------------
# Absolute flux
# ---------------------------------------------------------------------------
def subsystem_heatmap(flux_table: pd.DataFrame, out: Path) -> Path:
    """Mean |flux| per subsystem (rows) across conditions×timepoints (cols)."""
    d = flux_table.copy()
    d["absflux"] = d["flux"].abs()
    piv = (d.groupby(["subsystem", "condition", "timepoint"])["absflux"].mean()
           .reset_index())
    piv["col"] = piv["condition"] + " D" + piv["timepoint"].astype(str)
    mat = piv.pivot(index="subsystem", columns="col", values="absflux").fillna(0.0)
    # log scale for readability
    mat_log = np.log10(mat + 1.0)

    fig, ax = plt.subplots(figsize=(max(6, 0.5 * mat.shape[1]), max(4, 0.4 * mat.shape[0])))
    im = ax.imshow(mat_log.values, aspect="auto", cmap="viridis")
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(mat.columns, rotation=90, fontsize=7)
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels(mat.index, fontsize=7)
    ax.set_title("Mean |flux| by subsystem  (log10 scale)")
    fig.colorbar(im, ax=ax, label="log10(mean |flux| + 1)")
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out


def key_reaction_trajectories(flux_table: pd.DataFrame, out: Path,
                              reactions: list[str] | None = None) -> Path:
    reactions = reactions or list(C.METABOLITE_EXCHANGE.values()) + [C.BIOMASS_RXN]
    d = flux_table[flux_table["reaction_id"].isin(reactions)]
    grp = (d.groupby(["reaction_id", "condition", "timepoint"])["flux"]
           .mean().reset_index())
    fig, axes = plt.subplots(1, len(reactions), figsize=(3.2 * len(reactions), 3.2),
                             squeeze=False)
    for ax, rid in zip(axes[0], reactions):
        sub = grp[grp["reaction_id"] == rid]
        for cond, s in sub.groupby("condition"):
            s = s.sort_values("timepoint")
            ax.plot(s["timepoint"], s["flux"], marker="o", label=cond)
        ax.set_title(rid, fontsize=8)
        ax.set_xlabel("day")
        ax.axhline(0, color="grey", lw=0.6)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=7)
    axes[0][0].set_ylabel("mean flux")
    fig.suptitle("Key reaction flux trajectories (FB vs iFB)")
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Differential flux
# ---------------------------------------------------------------------------
def differential_flux(flux_table: pd.DataFrame, group_col: str,
                      group_a: str, group_b: str) -> pd.DataFrame:
    """Per-reaction contrast group_a vs group_b: log2 fold-change of mean flux
    magnitude + Welch t-test p-value across samples."""
    from scipy import stats  # local import; scipy optional

    d = flux_table.copy()
    key = _sample_col(d)
    d = d.assign(_sample=key)
    a = d[d[group_col] == group_a]
    b = d[d[group_col] == group_b]
    rows = []
    for rid, ga in a.groupby("reaction_id"):
        gb = b[b["reaction_id"] == rid]
        va, vb = ga["flux"].to_numpy(), gb["flux"].to_numpy()
        if len(va) < 2 or len(vb) < 2:
            continue
        ma, mb = np.nanmean(np.abs(va)), np.nanmean(np.abs(vb))
        l2fc = np.log2((ma + _EPS) / (mb + _EPS))
        try:
            _, p = stats.ttest_ind(va, vb, equal_var=False, nan_policy="omit")
        except Exception:
            p = np.nan
        rows.append({"reaction_id": rid, "subsystem": ga["subsystem"].iloc[0],
                     "mean_a": ma, "mean_b": mb, "log2fc": l2fc, "pvalue": p})
    return pd.DataFrame(rows)


def _fallback_diff(flux_table, group_col, group_a, group_b) -> pd.DataFrame:
    """log2FC only (no p-values) if scipy is unavailable."""
    d = flux_table.copy()
    d["absflux"] = d["flux"].abs()
    piv = d.groupby(["reaction_id", "subsystem", group_col])["absflux"].mean().reset_index()
    wide = piv.pivot_table(index=["reaction_id", "subsystem"], columns=group_col,
                           values="absflux").reset_index()
    if group_a not in wide or group_b not in wide:
        return pd.DataFrame(columns=["reaction_id", "subsystem", "log2fc", "pvalue"])
    wide["log2fc"] = np.log2((wide[group_a] + _EPS) / (wide[group_b] + _EPS))
    wide["pvalue"] = np.nan
    return wide.rename(columns={group_a: "mean_a", group_b: "mean_b"})


def volcano(diff: pd.DataFrame, out: Path, title: str) -> Path:
    fig, ax = plt.subplots(figsize=(6, 5))
    x = diff["log2fc"].to_numpy()
    p = diff["pvalue"].to_numpy()
    y = -np.log10(np.where(np.isnan(p), 1.0, p) + _EPS)
    sig = (np.abs(x) > 1) & (y > -np.log10(0.05))
    ax.scatter(x[~sig], y[~sig], s=8, c="lightgrey")
    ax.scatter(x[sig], y[sig], s=12, c="crimson")
    ax.axvline(1, ls="--", c="grey", lw=0.6); ax.axvline(-1, ls="--", c="grey", lw=0.6)
    ax.axhline(-np.log10(0.05), ls="--", c="grey", lw=0.6)
    ax.set_xlabel("log2 flux fold-change"); ax.set_ylabel("-log10 p")
    ax.set_title(title)
    fig.tight_layout(); fig.savefig(out, dpi=140); plt.close(fig)
    return out


def top_divergent_barplot(diff: pd.DataFrame, out: Path, title: str, n: int = 20) -> Path:
    d = diff.reindex(diff["log2fc"].abs().sort_values(ascending=False).index).head(n)
    d = d.iloc[::-1]
    fig, ax = plt.subplots(figsize=(6, max(3, 0.35 * len(d))))
    colors = ["crimson" if v > 0 else "steelblue" for v in d["log2fc"]]
    ax.barh(d["reaction_id"], d["log2fc"], color=colors)
    ax.axvline(0, color="black", lw=0.6)
    ax.set_xlabel("log2 flux fold-change"); ax.set_title(title)
    ax.tick_params(axis="y", labelsize=7)
    fig.tight_layout(); fig.savefig(out, dpi=140); plt.close(fig)
    return out


def subsystem_diff_heatmap(diff: pd.DataFrame, out: Path, title: str) -> Path:
    agg = diff.groupby("subsystem")["log2fc"].mean().sort_values()
    fig, ax = plt.subplots(figsize=(6, max(3, 0.35 * len(agg))))
    colors = ["crimson" if v > 0 else "steelblue" for v in agg]
    ax.barh(agg.index, agg.values, color=colors)
    ax.axvline(0, color="black", lw=0.6)
    ax.set_xlabel("mean log2 flux FC"); ax.set_title(title)
    ax.tick_params(axis="y", labelsize=7)
    fig.tight_layout(); fig.savefig(out, dpi=140); plt.close(fig)
    return out


def make_all_figures(flux_table: pd.DataFrame, out_dir: Path = C.FIG_DIR,
                     contrasts=None) -> list[Path]:
    """Generate all figures. ``contrasts`` is an optional list of
    ``(group_col, group_a, group_b)`` tuples; if omitted, sensible defaults are
    inferred (iFB-vs-FB and, when labelled, Responder-vs-NonResponder)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    figs = []
    figs.append(subsystem_heatmap(flux_table, out_dir / "flux_subsystem_heatmap.png"))
    figs.append(key_reaction_trajectories(flux_table, out_dir / "flux_trajectories.png"))

    if contrasts is None:
        conds = sorted(map(str, flux_table["condition"].dropna().unique()))
        if {"iFB", "FB"} <= set(conds):
            contrasts = [("condition", "iFB", "FB")]
        elif len(conds) == 2:                      # any two-condition design, e.g. ZeLa vs WT
            contrasts = [("condition", conds[1], conds[0])]
        else:
            contrasts = []
        if "responder" in flux_table.columns and flux_table["responder"].nunique() > 1:
            contrasts.append(("responder", "Responder", "NonResponder"))

    for col, a, b in contrasts:
        try:
            diff = differential_flux(flux_table, col, a, b)
        except ImportError:
            diff = _fallback_diff(flux_table, col, a, b)
        if diff.empty:
            continue
        tag = f"{a}_vs_{b}"
        title = f"{a} vs {b}"
        figs.append(volcano(diff, out_dir / f"diff_volcano_{tag}.png", title))
        figs.append(top_divergent_barplot(diff, out_dir / f"diff_top_{tag}.png", title))
        figs.append(subsystem_diff_heatmap(diff, out_dir / f"diff_subsystem_{tag}.png", title))
        diff.to_csv(out_dir / f"diff_{tag}.csv", index=False)
    return figs
