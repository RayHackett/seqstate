# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html


# ----------------------- read _version.py --------------------------------------
import ast
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

version_file = next(Path(sys.path[0]).rglob("_version.py"))
tree = ast.parse(version_file.read_text(encoding="utf-8"))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "seqstate"
copyright = "2026, Raymund Hackett"
author = "Raymund Hackett"
release = ast.literal_eval(tree.body[0].value)

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinxarg.ext",
    "myst_parser",
    "sphinx_rtd_theme",
]

napoleon_include_special_with_doc = True


def skip(app, what, name, obj, would_skip, options):
    if name == "__init__":
        return False
    return would_skip


def setup(app):
    app.connect("autodoc-skip-member", skip)


templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

# -- Options for InterSphinx -------------------------------------------------

default_role = "py:obj"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", (None, "python-inv.txt")),
    "[pyhmmer]": ("https://pyhmmer.readthedocs.io/en/latest/", None),
    "[pyrodigal]": ("https://pyrodigal.readthedocs.io/en/latest/", None),
    "[gb_io]": ("https://gb-io.readthedocs.io/en/latest/", None),
    # "[polars]": ("https://docs.pola.rs/api/python/stable/reference/", None),
}
