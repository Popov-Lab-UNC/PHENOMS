"""
Connectivity (graph) analysis for backbone N–O hydrogen bonds.

R&D module: builds a residue-residue network from bond-level occupancy/delta values
derived from SimulationSet pivot tables.

Exports an interactive HTML graph (no Jupyter required).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json
import re
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import networkx as nx

from phenoms.analysis import extract_residue_numbers
from phenoms import _html_embed


def bond_label_to_donor_acceptor_residue_numbers(bond_label: str) -> Tuple[int, int]:
    """
    'ALA1 -- GLY2' -> (1, 2)
    Using the convention: first residue is the donor.
    """

    return extract_residue_numbers(bond_label)


def _bond_label_to_donor_acceptor_tokens(bond_label: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Parse donor/acceptor residue tokens from a bond label like 'ALA1 -- GLY2'.
    """
    parts = [p.strip() for p in str(bond_label).split("--")]
    if len(parts) != 2:
        return None, None
    return parts[0], parts[1]


def _three_to_one(resname3: str) -> str:
    mapping = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }
    return mapping.get(str(resname3).upper(), "X")


def _parse_residue_token(token: Optional[str]) -> Tuple[Optional[str], Optional[int]]:
    """
    Parse residue token like 'ALA123' into ('ALA', 123).
    """
    if not token:
        return None, None
    m = re.match(r"^\s*([A-Za-z]{3})(-?\d+)\s*$", token)
    if not m:
        return None, None
    return m.group(1).upper(), int(m.group(2))


def bond_occupancy_from_pivot_tables(
    pivot_tables: Sequence[pd.DataFrame],
    bond_labels: Sequence[str],
) -> pd.Series:
    """
    Compute mean bond presence fraction over frames, averaged across replicates.

    Parameters
    ----------
    pivot_tables
        List of pivot tables: index=bond label, columns=Frame, values=0/1 (or counts).
    bond_labels
        Ordered list of bond labels to align to.

    Returns
    -------
    pd.Series
        Index=bond_labels, values=occupancy fraction in [0, 1] (approx).
    """
    if not pivot_tables:
        return pd.Series(index=bond_labels, data=np.zeros(len(bond_labels), dtype=float))

    occs = []
    for pt in pivot_tables:
        pt_aligned = pt.reindex(bond_labels, fill_value=0)
        binary = (pt_aligned > 0).astype(float)
        # Presence fraction per bond for this replicate.
        occ = binary.mean(axis=1)
        occs.append(occ)

    occ_mean = pd.concat(occs, axis=1).mean(axis=1)
    # Ensure correct index/order.
    return occ_mean.reindex(bond_labels)


def bond_delta_from_pivot_tables(
    pivot_tables_a: Sequence[pd.DataFrame],
    pivot_tables_b: Sequence[pd.DataFrame],
    bond_labels: Sequence[str],
    *,
    flip_difference: bool,
    impute_threshold: Optional[float] = None,
) -> pd.Series:
    """
    Compute bond-level delta values between two sets.

    The sign convention matches ComparisonSet.compare():
    - if flip_difference=False: delta = average_a - average_b
    - if flip_difference=True:  delta = average_b - average_a
    """
    occ_a = bond_occupancy_from_pivot_tables(pivot_tables_a, bond_labels)
    occ_b = bond_occupancy_from_pivot_tables(pivot_tables_b, bond_labels)

    delta = occ_a - occ_b
    if flip_difference:
        delta = -delta

    if impute_threshold is not None:
        low_both = (occ_a <= impute_threshold) & (occ_b <= impute_threshold)
        delta = delta.mask(low_both, 0.0)

    return delta.reindex(bond_labels)


@dataclass(frozen=True)
class ResidueGraph:
    nodes: List[dict]
    edges: List[dict]


@dataclass(frozen=True)
class CommunityGraphResult:
    graph: ResidueGraph
    community_by_node: Dict[str, int]
    modularity: float


def _color_from_signed_value(value: float, vmin: float, vmax: float) -> str:
    """
    Simple blue-white-red mapping without adding JS deps.
    For the HTML graph, we return a single color string.
    """
    # Normalize around 0.
    if value == 0 or vmin == vmax:
        return "#888888"
    if value > 0:
        return "#d62728"  # red
    return "#1f77b4"  # blue


def build_residue_graph_from_bond_deltas(
    bond_deltas: pd.Series,
    *,
    top_k_edges: int = 30,
    min_abs_edge_delta: Optional[float] = None,
    aggregate_edges_by_residue_pair: bool = True,
    directed: bool = False,
) -> ResidueGraph:
    """
    Build residue-residue graph from bond-level delta values.

    Nodes are residues (by donor/acceptor residue number).
    Edges correspond to donor->acceptor backbone H-bonds.
    """
    if bond_deltas.empty:
        return ResidueGraph(nodes=[], edges=[])

    # Compute edge rows per bond label.
    rows = []
    residue_token_map: Dict[int, str] = {}
    for bond_label, delta in bond_deltas.items():
        d_res, a_res = bond_label_to_donor_acceptor_residue_numbers(bond_label)
        delta_f = float(delta)
        if not np.isfinite(delta_f):
            continue
        donor_token, acceptor_token = _bond_label_to_donor_acceptor_tokens(str(bond_label))
        donor_name3, donor_num = _parse_residue_token(donor_token)
        acceptor_name3, acceptor_num = _parse_residue_token(acceptor_token)
        if donor_name3 is not None and donor_num is not None:
            residue_token_map[int(donor_num)] = donor_name3
        if acceptor_name3 is not None and acceptor_num is not None:
            residue_token_map[int(acceptor_num)] = acceptor_name3
        rows.append((d_res, a_res, delta_f, bond_label))
    edges_df = pd.DataFrame(rows, columns=["donor_res", "acceptor_res", "delta", "bond_label"])

    edges_df["abs_delta"] = edges_df["delta"].abs()
    edges_df = edges_df.replace([np.inf, -np.inf], np.nan).dropna(subset=["delta", "abs_delta"])

    if min_abs_edge_delta is not None:
        edges_df = edges_df[edges_df["abs_delta"] >= float(min_abs_edge_delta)]

    if edges_df.empty:
        return ResidueGraph(nodes=[], edges=[])

    if aggregate_edges_by_residue_pair:
        # Aggregate multiple bond labels (if any) into a single residue pair edge.
        grouped = edges_df.groupby(["donor_res", "acceptor_res"], as_index=False).agg(
            delta_sum=("delta", "sum"),
            abs_delta_sum=("abs_delta", "sum"),
            bond_labels=("bond_label", lambda xs: ", ".join(list(xs)[:6])),
        )
        edges_df = grouped.rename(columns={"delta_sum": "delta", "abs_delta_sum": "abs_delta"})
    else:
        edges_df = edges_df.rename(columns={"delta": "delta", "abs_delta": "abs_delta"})

    # Select top edges by abs(delta)
    edges_df = edges_df.sort_values("abs_delta", ascending=False)
    if top_k_edges is not None:
        edges_df = edges_df.head(int(top_k_edges))

    if edges_df.empty:
        return ResidueGraph(nodes=[], edges=[])

    max_abs = float(edges_df["abs_delta"].max()) if len(edges_df) else 1.0
    max_abs = max(max_abs, 1e-12)

    # Node score: outgoing sum(delta) using the (donor) role.
    node_scores = edges_df.groupby("donor_res")["delta"].sum().to_dict()

    residues = set(edges_df["donor_res"].tolist()) | set(edges_df["acceptor_res"].tolist())
    residues = sorted(residues)

    nodes = []
    for r in residues:
        score = float(node_scores.get(r, 0.0))
        if not np.isfinite(score):
            continue
        size = 10.0 + 30.0 * (abs(score) / (abs(max(node_scores.values(), default=0.0)) + 1e-12))
        color = _color_from_signed_value(score, -1.0, 1.0)
        name3 = residue_token_map.get(int(r))
        one = _three_to_one(name3) if name3 else "X"
        hover = f"{one}{int(r)}"
        if name3:
            hover = f"{hover} ({name3}{int(r)})"
        nodes.append(
            {
                "id": str(r),
                "label": str(r),
                "size": size,
                "color": color,
                "title": f"Residue: {hover}<br/>Node score: {score:.4f}",
                "residue_name3": name3 if name3 else "",
                "residue_one_letter": one,
                "residue_number": int(r),
            }
        )

    edges = []
    for _, row in edges_df.iterrows():
        donor_res = int(row["donor_res"])
        acceptor_res = int(row["acceptor_res"])
        delta = float(row["delta"])
        abs_delta = float(row["abs_delta"])
        if not np.isfinite(delta) or not np.isfinite(abs_delta):
            continue
        width = 1.0 + 8.0 * (abs_delta / max_abs)
        color = _color_from_signed_value(delta, -1.0, 1.0)

        label = f"{donor_res} -> {acceptor_res} (Δ={delta:.3f})"
        title = f"Donor {donor_res}, Acceptor {acceptor_res}<br/>Δ occupancy = {delta:.4f}<br/>Bonds: {row.get('bond_labels','')}"

        edge = {
            "from": str(donor_res),
            "to": str(acceptor_res),
            "width": width,
            "color": color,
            "weight": delta,
            "abs_weight": abs_delta,
            "title": title,
        }
        if directed:
            edge["arrows"] = "to"

        edges.append(edge)

    return ResidueGraph(nodes=nodes, edges=edges)


def _build_networkx_graph(graph: ResidueGraph, *, directed: bool = False) -> nx.Graph:
    """
    Convert ResidueGraph payload to a NetworkX graph.

    Parameters
    ----------
    graph
        Graph payload as used for vis-network export.
    directed
        If True, build a DiGraph. Community detection uses an undirected copy.
    """
    g = nx.DiGraph() if directed else nx.Graph()
    for node in graph.nodes:
        node_id = str(node.get("id"))
        g.add_node(node_id, **node)
    for edge in graph.edges:
        src = str(edge.get("from"))
        dst = str(edge.get("to"))
        weight = float(edge.get("weight", edge.get("width", 1.0)))
        abs_weight = float(edge.get("abs_weight", abs(weight)))
        attrs = dict(edge)
        attrs["weight"] = weight
        attrs["abs_weight"] = abs_weight
        g.add_edge(src, dst, **attrs)
    return g


def assign_greedy_modularity_communities(
    graph: ResidueGraph,
    *,
    directed: bool = False,
    weight_attr: str = "abs_weight",
    resolution: float = 1.0,
) -> CommunityGraphResult:
    """
    Compute residue communities from a ResidueGraph using greedy modularity.

    Returns per-node community IDs and modularity score.
    """
    nx_graph = _build_networkx_graph(graph, directed=directed)
    if nx_graph.number_of_nodes() == 0:
        return CommunityGraphResult(graph=graph, community_by_node={}, modularity=0.0)

    undirected = nx_graph.to_undirected()
    communities = list(
        nx.algorithms.community.greedy_modularity_communities(
            undirected,
            weight=weight_attr,
            resolution=resolution,
        )
    )
    community_by_node: Dict[str, int] = {}
    for cid, nodes in enumerate(sorted(communities, key=lambda s: (-len(s), sorted(s)))):
        for node_id in sorted(nodes):
            community_by_node[str(node_id)] = int(cid)

    modularity = float(
        nx.algorithms.community.modularity(
            undirected,
            communities,
            weight=weight_attr,
            resolution=resolution,
        )
    ) if communities else 0.0
    return CommunityGraphResult(
        graph=graph,
        community_by_node=community_by_node,
        modularity=modularity,
    )


def _community_palette() -> List[str]:
    # Distinct but lightweight palette with no extra deps.
    return [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    ]


def _color_for_community(community_id: int) -> str:
    palette = _community_palette()
    return palette[int(community_id) % len(palette)]


def apply_community_annotations(
    graph: ResidueGraph,
    community_by_node: Dict[str, int],
) -> ResidueGraph:
    """
    Return a graph with node/edge annotations derived from community assignments.
    """
    nodes = []
    for node in graph.nodes:
        node_id = str(node.get("id"))
        cid = community_by_node.get(node_id, -1)
        colored = dict(node)
        colored["community_id"] = int(cid)
        if cid >= 0:
            colored["color"] = _color_for_community(cid)
        title_base = colored.get("title", f"Residue {node_id}")
        colored["title"] = f"{title_base}<br/>Community: {cid}"
        nodes.append(colored)

    edges = []
    for edge in graph.edges:
        src = str(edge.get("from"))
        dst = str(edge.get("to"))
        src_c = int(community_by_node.get(src, -1))
        dst_c = int(community_by_node.get(dst, -1))
        annotated = dict(edge)
        annotated["source_community"] = src_c
        annotated["target_community"] = dst_c
        annotated["is_inter_community"] = bool(src_c != dst_c)
        title = annotated.get("title", "")
        annotated["title"] = (
            f"{title}<br/>Communities: {src_c} -> {dst_c}"
            if title
            else f"Communities: {src_c} -> {dst_c}"
        )
        edges.append(annotated)

    return ResidueGraph(nodes=nodes, edges=edges)


def build_residue_graph_with_communities(
    bond_values: pd.Series,
    *,
    top_k_edges: int = 30,
    min_abs_edge_delta: Optional[float] = None,
    aggregate_edges_by_residue_pair: bool = True,
    directed: bool = False,
    resolution: float = 1.0,
) -> CommunityGraphResult:
    """
    Build residue graph and assign community IDs with greedy modularity.
    """
    base_graph = build_residue_graph_from_bond_deltas(
        bond_values,
        top_k_edges=top_k_edges,
        min_abs_edge_delta=min_abs_edge_delta,
        aggregate_edges_by_residue_pair=aggregate_edges_by_residue_pair,
        directed=directed,
    )
    result = assign_greedy_modularity_communities(
        base_graph,
        directed=directed,
        resolution=resolution,
    )
    return CommunityGraphResult(
        graph=apply_community_annotations(base_graph, result.community_by_node),
        community_by_node=result.community_by_node,
        modularity=result.modularity,
    )


def community_tables(
    graph_result: CommunityGraphResult,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build node-level and community-level summary tables.
    """
    graph = graph_result.graph
    community_by_node = graph_result.community_by_node
    nodes_rows = []
    node_size_by_id = {str(n.get("id")): float(n.get("size", 0.0)) for n in graph.nodes}
    for node_id, cid in sorted(community_by_node.items(), key=lambda x: (x[1], int(x[0]))):
        nodes_rows.append(
            {
                "residue": int(node_id),
                "community_id": int(cid),
                "node_size": node_size_by_id.get(node_id, 0.0),
            }
        )
    nodes_df = pd.DataFrame(nodes_rows, columns=["residue", "community_id", "node_size"])

    if not community_by_node:
        summary_df = pd.DataFrame(
            [{"community_id": -1, "num_nodes": 0, "intra_edge_weight_sum": 0.0, "modularity": 0.0}]
        )
        return nodes_df, summary_df

    intra_weight = {}
    for edge in graph.edges:
        src = str(edge.get("from"))
        dst = str(edge.get("to"))
        src_c = community_by_node.get(src, -1)
        dst_c = community_by_node.get(dst, -1)
        if src_c < 0 or dst_c < 0:
            continue
        if src_c == dst_c:
            intra_weight[src_c] = intra_weight.get(src_c, 0.0) + float(
                edge.get("abs_weight", edge.get("width", 1.0))
            )

    counts = {}
    for _, cid in community_by_node.items():
        counts[cid] = counts.get(cid, 0) + 1

    summary_rows = []
    for cid in sorted(counts):
        summary_rows.append(
            {
                "community_id": int(cid),
                "num_nodes": int(counts[cid]),
                "intra_edge_weight_sum": float(intra_weight.get(cid, 0.0)),
                "modularity": float(graph_result.modularity),
            }
        )
    summary_df = pd.DataFrame(
        summary_rows,
        columns=["community_id", "num_nodes", "intra_edge_weight_sum", "modularity"],
    )
    return nodes_df, summary_df


def export_residue_graph_html(
    graph: ResidueGraph,
    output_html_path: str | Path,
    *,
    title: str = "Residue H-bond connectivity graph",
    directed: bool = False,
    summary_stats: Optional[dict] = None,
    show_delta_legend: bool = True,
) -> None:
    """
    Export an interactive HTML graph using vis-network (loaded from a CDN).
    """
    output_html_path = Path(output_html_path)

    nodes_json = json.dumps(graph.nodes, allow_nan=False)
    edges_json = json.dumps(graph.edges, allow_nan=False)

    # Simple legend for sign.
    legend = """
      <div style="margin: 10px 0; font-family: sans-serif; font-size: 14px;">
        <span style="display:inline-block; width:12px; height:12px; background:#d62728; margin-right:6px;"></span>Edge Δ > 0
        <span style="display:inline-block; width:12px; height:12px; background:#1f77b4; margin-left:14px; margin-right:6px;"></span>Edge Δ < 0
      </div>
    """.strip()
    if not show_delta_legend:
        legend = ""
    summary_stats = summary_stats or {}
    summary_lines = [
        f"Nodes: {summary_stats.get('n_nodes', len(graph.nodes))}",
        f"Edges: {summary_stats.get('n_edges', len(graph.edges))}",
    ]
    if "n_communities" in summary_stats:
        summary_lines.append(f"Communities: {summary_stats['n_communities']}")
    if "modularity" in summary_stats:
        summary_lines.append(f"Modularity: {float(summary_stats['modularity']):.3f}")
    summary_html = (
        "<div id=\"summary\" style=\"margin: 8px 0; font-family:sans-serif; font-size:13px;\">"
        + " | ".join(summary_lines)
        + "</div>"
    )

    # Embed vis-network JS directly so the HTML doesn't depend on external networks/CSP.
    vis_url = "https://unpkg.com/vis-network/standalone/umd/vis-network.min.js"
    vis_script_tag = _html_embed.embed_js_from_url(vis_url)

    # Lightweight fallback output (always visible) so results are still usable
    # if vis-network fails to render for any reason.
    fallback_rows = []
    for e in graph.edges[:50]:
        # Our edge title already contains the key info.
        fallback_rows.append(f"<li>{e.get('title','')}</li>")
    fallback_html = (
        "<div id=\"fallback\" style=\"margin: 10px 0; font-family:sans-serif; font-size:13px;\">"
        "<b>Graph fallback (top edges)</b>"
        "<ul style=\"margin-top:6px;\">"
        + "".join(fallback_rows)
        + "</ul></div>"
    )

    html = f"""<!doctype html>
<html>
  <head>
    <meta charset="utf-8" />
    <title>{title}</title>
    <style>
      body {{ margin: 0; padding: 12px; }}
      #mynetwork {{ width: 100%; height: 700px; border: 1px solid #ddd; border-radius: 6px; }}
      #error {{ margin-top: 10px; color: #b00020; font-size: 13px; white-space: pre-wrap; }}
    </style>
    {vis_script_tag}
  </head>
  <body>
    <h3 style="margin:0 0 6px 0; font-family: sans-serif;">{title}</h3>
    {legend}
    {summary_html}
    <div style="margin: 8px 0; font-family:sans-serif; font-size:13px;">
      <label><input id="toggleInter" type="checkbox" /> Show inter-community edges only</label>
    </div>
    <div id="status" style="font-family:sans-serif; font-size:13px; margin-bottom:6px;">Loading vis-network...</div>
    <div id="error"></div>
    {fallback_html}
    <div id="mynetwork"></div>
    <script>
      (function() {{
        const statusEl = document.getElementById("status");
        const errorEl = document.getElementById("error");
        const fallbackEl = document.getElementById("fallback");
        try {{
          if (typeof vis === "undefined") {{
            statusEl.textContent = "vis-network failed to load (CDN blocked).";
            errorEl.textContent = "window.vis is undefined.";
            return;
          }}

          statusEl.textContent = "Rendering graph...";
          const nodes = new vis.DataSet({nodes_json});
          const edgesRaw = {edges_json};
          const edges = new vis.DataSet(edgesRaw);
          const container = document.getElementById('mynetwork');
          const data = {{ nodes, edges }};
          // vis-network v10+: edges.scaling / smooth must be objects if set; booleans throw mergeTarget errors.
          const options = {{
            autoResize: true,
            physics: {{ stabilization: {{ iterations: 200 }} }},
            nodes: {{ shape: "dot", font: {{ size: 14 }} }},
            edges: {{ width: 1 }}
          }};
          new vis.Network(container, data, options);
          const toggleInter = document.getElementById("toggleInter");
          if (toggleInter) {{
            toggleInter.addEventListener("change", function() {{
              const filtered = this.checked
                ? edgesRaw.filter(e => Boolean(e.is_inter_community))
                : edgesRaw;
              edges.clear();
              edges.add(filtered);
            }});
          }}
          statusEl.textContent = "Done.";
          if (fallbackEl) fallbackEl.style.display = "none";
        }} catch (err) {{
          errorEl.textContent = "Graph render error:\\n" + String(err);
          statusEl.textContent = "Failed to render graph.";
        }}
      }})();
    </script>
  </body>
</html>"""

    output_html_path.write_text(html, encoding="utf-8")

