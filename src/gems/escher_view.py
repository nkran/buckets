"""Render Escher maps to self-contained, offline HTML with flux / expression overlays.

Uses escher-python's ``Builder`` to generate the map HTML (it base64-embeds the map +
data), then swaps the CDN ``escher.min.js`` reference for the bundle shipped inside the
installed ``escher`` package — so the result runs with no network access. No cobra or
solver is needed at render time; reaction data comes from the exported matrices.
"""

from __future__ import annotations

import functools
import json
import re
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# reaction_scale presets (Escher maps one data vector to BOTH color and size)
# ---------------------------------------------------------------------------
# Flux: diverging, centred at 0 -> thickness encodes |flux|, colour encodes direction.
FLUX_SCALE = [
    {"type": "min", "color": "#2166ac", "size": 34},
    {"type": "value", "value": 0, "color": "#e8e8e8", "size": 6},
    {"type": "max", "color": "#b2182b", "size": 34},
]
# Expression: sequential, non-negative -> thickness + colour both encode level.
EXPRESSION_SCALE = [
    {"type": "min", "color": "#f7fbff", "size": 8},
    {"type": "median", "color": "#6baed6", "size": 22},
    {"type": "max", "color": "#08519c", "size": 40},
]

_CDN_RE = re.compile(r'<script[^>]+src="[^"]*escher(?:\.min)?\.js"[^>]*>\s*</script>')


# Prefer the vendored UMD `dist` build (exposes a global `window.escher`, which the
# save_html template calls). The escher package's `static/escher.min.js` is a different
# (non-global) build and does NOT work when inlined standalone.
_VENDORED_JS = Path(__file__).resolve().parents[2] / "escher_maps" / "lib" / "escher.min.js"


@functools.lru_cache(maxsize=1)
def _escher_script_tag() -> str:
    """Return a `<script src="data:...">` tag embedding the escher.min.js bundle.

    The bundle must be loaded via a ``src`` (not inlined) because Escher's webpack
    runtime auto-detects its publicPath from ``document.currentScript.src`` and throws
    "Automatic publicPath is not supported" for an inline script. A base64 ``data:``
    URL gives it a src while staying fully offline.
    """
    if not _VENDORED_JS.exists():
        raise FileNotFoundError(
            f"Vendored Escher JS not found at {_VENDORED_JS}. Run "
            f"`python scripts/prepare_escher_maps.py` (it downloads dist/escher.min.js)."
        )
    import base64
    b64 = base64.b64encode(_VENDORED_JS.read_bytes()).decode("ascii")
    return f'<script src="data:text/javascript;base64,{b64}"></script>'


def reaction_data_for_sample(matrix: pd.DataFrame, sample_id: str) -> dict[str, float]:
    """Column ``sample_id`` of a reaction × sample matrix -> {reaction_id: value},
    dropping NaNs (Escher then styles those reactions as 'no data')."""
    if sample_id not in matrix.columns:
        raise KeyError(f"sample {sample_id!r} not in matrix columns")
    s = matrix[sample_id].dropna()
    return {str(k): float(v) for k, v in s.items()}


def coverage_for_map(escher_map: list, data: dict[str, float]) -> float:
    """Fraction of the map's reactions that have a value in ``data``."""
    rids = [r["bigg_id"] for r in escher_map[1]["reactions"].values()]
    if not rids:
        return 0.0
    return sum(1 for rid in rids if rid in data) / len(rids)


def build_escher_html(map_path, reaction_data: dict[str, float], *, mode: str = "flux",
                      scale: list | None = None, height: int = 640,
                      title: str | None = None) -> str:
    """Return self-contained offline HTML for a map coloured/sized by ``reaction_data``.

    ``mode`` selects the default scale ('flux' or 'expression'); ``scale`` overrides it.
    """
    from escher import Builder

    if scale is None:
        scale = FLUX_SCALE if mode == "flux" else EXPRESSION_SCALE
    map_json = Path(map_path).read_text()

    builder = Builder(
        map_json=map_json,
        reaction_data=reaction_data or None,
        reaction_scale=scale,
        reaction_styles=["color", "size", "text"],
        reaction_no_data_color="#dcdcdc",
        reaction_no_data_size=6,
        hide_secondary_metabolites=False,
        never_ask_before_quit=True,
        menu="zoom",
        enable_editing=False,
    )
    with tempfile.NamedTemporaryFile("w+", suffix=".html", delete=False) as fh:
        builder.save_html(fh.name)
        html = Path(fh.name).read_text()
    Path(fh.name).unlink(missing_ok=True)

    # make offline: swap the CDN <script> for the vendored bundle (data: URL).
    # Function replacement avoids regex escape-template issues.
    tag = _escher_script_tag()
    html, n = _CDN_RE.subn(lambda _m: tag, html)
    if n == 0:  # unexpected template — fall back to injecting before </head>
        html = html.replace("</head>", tag + "</head>", 1)

    # Escher's #map-container is height:100% and needs html/body to have height;
    # inside an embedding iframe they collapse to 0, so pin an explicit height.
    style = ("<style>html,body{height:100%;margin:0;padding:0}"
             f"#map-container{{height:{height}px;width:100%}}</style>")
    if "</head>" in html:
        html = html.replace("</head>", style + "</head>", 1)
    else:
        html = style + html

    # Escher builds asynchronously and opens at 100% (often off-content); its
    # "Zoom to canvas" tool misbehaves for our generated canvas, so fit the zoom
    # group to the real content bounding box directly once it has rendered.
    # Fit the view to the CONTENT (union of zoom-group children, excluding the big
    # background canvas rect) once Escher has drawn it. Fitting the canvas instead
    # would shrink content to sub-pixel on the large RECON1 maps.
    autofit = (
        "<script>(function(){var n=0,last=-1,t=setInterval(function(){"
        "var zg=document.querySelector('.zoom-g');"
        "var mc=document.getElementById('map-container');"
        "if(!(zg&&mc&&zg.querySelector('.reaction'))){if(++n>80)clearInterval(t);return;}"
        "var cnt=zg.querySelectorAll('.reaction').length;"
        "if(cnt!==last){last=cnt;if(++n>80)clearInterval(t);return;}"  # wait until stable
        "var x0=1e9,y0=1e9,x1=-1e9,y1=-1e9,ok=false,k=zg.children;"
        "for(var i=0;i<k.length;i++){var c=k[i];"
        "if((c.getAttribute('class')||'').indexOf('canvas')>=0)continue;"
        "var b;try{b=c.getBBox()}catch(e){continue}"
        "if(!b||(b.width===0&&b.height===0))continue;"
        "ok=true;x0=Math.min(x0,b.x);y0=Math.min(y0,b.y);"
        "x1=Math.max(x1,b.x+b.width);y1=Math.max(y1,b.y+b.height);}"
        "if(ok){var bw=x1-x0,bh=y1-y0,cw=mc.clientWidth,ch=mc.clientHeight;"
        "var s=0.92*Math.min(cw/bw,ch/bh);"
        "var tx=cw/2-s*(x0+bw/2),ty=ch/2-s*(y0+bh/2);"
        "zg.setAttribute('transform','translate('+tx+','+ty+') scale('+s+')');"
        "clearInterval(t);}if(++n>80)clearInterval(t);},150);})();</script>"
    )
    html = html.replace("</body>", autofit + "</body>", 1) if "</body>" in html else html + autofit
    return html


def show_map(map_path, reaction_data: dict[str, float], *, mode: str = "flux",
            scale: list | None = None, height: int = 650, title: str | None = None):
    """Display a map inline in Jupyter as a self-contained offline iframe.

    Prefer this over ``escher.Builder`` in the notebook: Escher's ipywidget frontend
    only registers as a classic-notebook **nbextension**, which JupyterLab / Notebook 7
    do not load — the widget then hangs forever on "Loading widget...". This renders
    the same HTML produced by :func:`build_escher_html` (the one used by the Streamlit
    app) inside an ``<iframe srcdoc=...>``, which works in any Jupyter frontend because
    the iframe gets its own document context where the embedded script actually runs.
    """
    import html as _html

    from IPython.display import HTML

    page = build_escher_html(map_path, reaction_data, mode=mode, scale=scale, height=height,
                             title=title)
    srcdoc = _html.escape(page, quote=True)
    return HTML(f'<iframe srcdoc="{srcdoc}" style="width:100%;height:{height}px;'
               f'border:1px solid #ddd;border-radius:4px;"></iframe>')


def scale_bounds(reaction_data: dict[str, float], mode: str = "flux") -> tuple[float, float]:
    """Convenience min/max (abs for flux) for legends in the app."""
    vals = np.array(list(reaction_data.values()), dtype=float)
    if vals.size == 0:
        return (0.0, 0.0)
    if mode == "flux":
        return (float(-np.nanmax(np.abs(vals))), float(np.nanmax(np.abs(vals))))
    return (float(np.nanmin(vals)), float(np.nanmax(vals)))
