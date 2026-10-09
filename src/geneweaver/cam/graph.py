"""The CAM agent as a LangGraph StateGraph.

    load_genome -> call_genes -> translate -> qc_check --(flagged?)--> mark_for_review -> export
                                                     \\--(clean)----------------------> export

Checkpointing (SqliteSaver) records the state after every node, so a run can be
inspected, resumed or replayed later. Later, CAM becomes a subgraph of the full
GeneWeaver graph and the 'review' branch feeds the PAM agent.
"""
from __future__ import annotations

import sqlite3
from pathlib import Path

from langgraph.graph import END, START, StateGraph

from . import nodes
from .state import CAMState


def build_cam_graph():
    g = StateGraph(CAMState)
    g.add_node("load_genome", nodes.load_genome)
    g.add_node("call_genes", nodes.call_genes)
    g.add_node("translate", nodes.translate)
    g.add_node("qc_check", nodes.qc_check)
    g.add_node("mark_for_review", nodes.mark_for_review)
    g.add_node("export", nodes.export)

    g.add_edge(START, "load_genome")
    g.add_edge("load_genome", "call_genes")
    g.add_edge("call_genes", "translate")
    g.add_edge("translate", "qc_check")
    g.add_conditional_edges(
        "qc_check", nodes.route_after_qc, {"review": "mark_for_review", "clean": "export"}
    )
    g.add_edge("mark_for_review", "export")
    g.add_edge("export", END)
    return g


def compile_cam(checkpoint_db: str | Path | None = None):
    """Compile the graph, optionally with a SQLite checkpointer."""
    graph = build_cam_graph()
    if checkpoint_db is None:
        return graph.compile()
    from langgraph.checkpoint.sqlite import SqliteSaver

    Path(checkpoint_db).parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(checkpoint_db), check_same_thread=False)
    return graph.compile(checkpointer=SqliteSaver(conn))


def mermaid() -> str:
    """Mermaid source for the architecture diagram (README / dissertation figure)."""
    return build_cam_graph().compile().get_graph().draw_mermaid()
