"""
Shared pytest fixtures.
"""

import pytest


@pytest.fixture(autouse=True)
def _isolate_phenoms_output_dir(tmp_path, monkeypatch):
    """
    SimulationSet now writes a standard output bundle by default. Point
    PHENOMS_OUTPUT_DIR at a per-test tmp dir so that default write never
    touches the repo working directory or leaks between tests.
    """
    monkeypatch.setenv("PHENOMS_OUTPUT_DIR", str(tmp_path / "phenom_outputs"))


@pytest.fixture(autouse=True)
def _stub_html_cdn_embeds(monkeypatch):
    """
    The default output bundle now writes standalone HTML viewers (structure
    viewer, connectivity graph) that embed a JS library fetched from a CDN.
    Stub that out in tests so the suite stays fast, deterministic, and
    network-independent -- the embed helper itself is unit-tested separately
    (or trusted as a thin, already-exercised download+fallback wrapper).
    """
    monkeypatch.setattr(
        "phenoms._html_embed.embed_js_from_url",
        lambda url, **kwargs: f'<script src="{url}"></script>',
    )
