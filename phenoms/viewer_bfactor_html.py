"""
Export a standalone HTML 3D viewer that colors a static-frame PDB by B-factor.

R&D module: uses NGL from a CDN.
"""

from __future__ import annotations

from pathlib import Path
import base64

from phenoms import _html_embed


def bfactor_min_max_from_pdb(pdb_path: str | Path) -> tuple[float, float]:
    pdb_path = Path(pdb_path)
    vals = []
    with pdb_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 66:
                try:
                    vals.append(float(line[60:66]))
                except ValueError:
                    continue
    if not vals:
        return (0.0, 0.0)
    return (min(vals), max(vals))


def export_bfactor_colored_pdb_viewer_html(
    pdb_path: str | Path,
    output_html_path: str | Path,
    *,
    title: str = "PDB colored by B-factor",
    background_color: str = "white",
) -> None:
    """
    Export a standalone HTML that loads `pdb_path` (embedded) and colors by B-factor.
    """
    pdb_path = Path(pdb_path)
    output_html_path = Path(output_html_path)

    pdb_bytes = pdb_path.read_bytes()
    pdb_b64 = base64.b64encode(pdb_bytes).decode("ascii")
    bmin, bmax = bfactor_min_max_from_pdb(pdb_path)

    # Embed NGL JS directly so the viewer does not depend on external networks/CSP.
    # Use a stable NGL version that exists on unpkg.
    ngl_url = "https://unpkg.com/ngl@2.4.0/dist/ngl.js"
    ngl_script_tag = _html_embed.embed_js_from_url(ngl_url)

    html = f"""<!doctype html>
<html>
  <head>
    <meta charset="utf-8" />
    <title>{title}</title>
    <style>
      body {{ margin: 0; padding: 0; background: {background_color}; font-family: sans-serif; }}
      #viewport {{ width: 100vw; height: 92vh; }}
      #legend {{ padding: 8px 12px; font-size: 14px; }}
      #status {{ padding: 0 12px 10px 12px; color: #333; font-size: 13px; }}
      #error {{ padding: 0 12px 10px 12px; color: #b00020; font-size: 13px; white-space: pre-wrap; }}
    </style>
    {ngl_script_tag}
  </head>
  <body>
    <div id="legend">
      <b>{title}</b><br/>
      B-factor range: {bmin:.3f} .. {bmax:.3f}
    </div>
    <div id="status">Loading NGL + PDB...</div>
    <div id="error"></div>
    <div id="viewport"></div>
    <script>
      (function() {{
        const statusEl = document.getElementById("status");
        const errorEl = document.getElementById("error");

        try {{
          if (typeof NGL === "undefined") {{
            statusEl.textContent = "Failed to load NGL (CDN blocked).";
            errorEl.textContent = "window.NGL is undefined.";
            return;
          }}

          const pdbBase64 = "{pdb_b64}";
          const pdbText = atob(pdbBase64);

          const stage = new NGL.Stage("viewport", {{ backgroundColor: "{background_color}" }});
          const blob = new Blob([pdbText], {{ type: "text/plain" }});
          const url = URL.createObjectURL(blob);

          stage.loadFile(url, {{ ext: "pdb" }})
            .then(function(component) {{
              statusEl.textContent = "PDB loaded. Rendering...";

              // Explicit blue (-1) -> white (0) -> red (+1); NGL built-in "bwr" can read yellow-ish.
              var bwrId = NGL.ColormakerRegistry.addScheme(function () {{
                this.atomColor = function (atom) {{
                  var v = atom.bfactor;
                  if (v === undefined || v === null || isNaN(v)) v = 0;
                  v = Math.max(-1, Math.min(1, v));
                  var t = (v + 1) / 2;
                  var r, g, b;
                  if (t <= 0.5) {{
                    var s = t * 2;
                    r = Math.round(255 * s);
                    g = Math.round(255 * s);
                    b = 255;
                  }} else {{
                    var s2 = (t - 0.5) * 2;
                    r = 255;
                    g = Math.round(255 * (1 - s2));
                    b = Math.round(255 * (1 - s2));
                  }}
                  return (r << 16) | (g << 8) | b;
                }};
              }});
              component.addRepresentation("cartoon", {{ color: bwrId }});
              component.addRepresentation("licorice", {{ sele: "not hydrogen", opacity: 0.08 }});

              stage.autoView();
              statusEl.textContent = "Done.";
            }})
            .catch(function(err) {{
              errorEl.textContent = "Failed to load/render PDB via NGL:\\n" + String(err);
              statusEl.textContent = "Render failed.";
            }});
        }} catch (err) {{
          errorEl.textContent = "Unexpected viewer error:\\n" + String(err);
          statusEl.textContent = "Viewer failed.";
        }}
      }})();
    </script>
  </body>
</html>"""

    output_html_path.write_text(html, encoding="utf-8")

