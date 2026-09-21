# -- Project information -----------------------------------------------------
"""Sphinx configuration for PHENOMS documentation."""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

try:
    from phenoms import __version__ as release
except Exception:  # pragma: no cover - docs build without full runtime deps
    release = "0.2.0"

project = "PHENOMS"
author = "PHENOMS contributors"
copyright = f"{datetime.now():%Y}, {author}"
version = release

# -- General configuration ---------------------------------------------------
extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
]

autosummary_generate = True
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_param = True
napoleon_use_rtype = True

# Heavy / optional native deps: mock so autodoc can import the package in CI
# even when mdtraj / rust extension are unavailable in a minimal docs env.
autodoc_mock_imports = [
    "mdtraj",
    "polars",
    "pandas",
    "numpy",
    "matplotlib",
    "seaborn",
    "sklearn",
    "networkx",
    "joblib",
    "tqdm",
    "Bio",
    "phenoms.phenoms_hbond_rs",
]

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
master_doc = "index"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
}

# -- Options for HTML output -------------------------------------------------
html_theme = "shibuya"
html_static_path = ["_static"]
html_title = "PHENOMS"
html_copy_source = False
html_show_sourcelink = False

html_theme_options = {
    "accent_color": "violet",
    "github_url": "https://github.com/Popov-Lab-UNC/PHENOMS",
    "globaltoc_expand_depth": 1,
    "nav_links": [
        {"title": "Installation", "url": "installation"},
        {"title": "Usage", "url": "usage"},
        {"title": "CLI", "url": "cli"},
        {"title": "Outputs", "url": "outputs"},
        {"title": "Data", "url": "data"},
        {"title": "API", "url": "api/index"},
    ],
}

html_context = {
    "source_type": "github",
    "source_user": "Popov-Lab-UNC",
    "source_repo": "PHENOMS",
    "source_version": "main",
    "source_docs_path": "/docs/",
}
