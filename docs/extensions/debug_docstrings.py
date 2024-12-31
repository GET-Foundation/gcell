# Just do the following to see the rst of a function:
# rm ./_build/doctrees/api/generated/gcell.<what you want>.doctree; DEBUG=1 make html
import os

import sphinx.ext.napoleon
from sphinx.application import Sphinx

_pd_orig = sphinx.ext.napoleon._process_docstring


def pd_new(app, what, name, obj, options, lines):  # noqa: PLR0917
    _pd_orig(app, what, name, obj, options, lines)
    print(*lines, sep="\n")


def setup(app: Sphinx):
    if os.environ.get("DEBUG") is not None:
        sphinx.ext.napoleon._process_docstring = pd_new
