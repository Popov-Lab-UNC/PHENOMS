"""
Shared helper for embedding third-party JS libraries into standalone HTML exports
(the NGL structure viewer and the vis-network connectivity graph).
"""

from __future__ import annotations

import ssl
import urllib.request
from functools import lru_cache


@lru_cache(maxsize=None)
def embed_js_from_url(url: str, *, timeout_s: float = 15.0) -> str:
    """
    Download JS and embed it as an inline ``<script>`` tag so the exported HTML
    is standalone. Falls back to a remote ``<script src=...>`` tag if the
    download fails (e.g. no network, CDN blocked) -- the export still
    succeeds, it just needs network access at view time in that case.

    Cached per URL for the life of the process, since the default output
    bundle can call this once per run/compare.
    """
    for ctx in (None, ssl._create_unverified_context()):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "phenoms/1.0"})
            with urllib.request.urlopen(req, timeout=timeout_s, context=ctx) as resp:
                js_text = resp.read().decode("utf-8", errors="ignore")
            if not js_text.strip():
                raise RuntimeError("Downloaded JS is empty")
            return f"<script>\n{js_text}\n</script>"
        except Exception:
            continue
    return f'<script src="{url}"></script>'
