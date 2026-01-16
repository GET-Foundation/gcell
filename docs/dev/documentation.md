# Documentation

(building-the-docs)=

## Building the docs

To build the docs, run `uv run sphinx-build docs docs/_build`.
Afterwards, open {file}`docs/_build/html/index.html` in your browser.

Your browser and Sphinx cache docs which have been built previously.
Sometimes these caches are not invalidated when you've updated the docs.
If docs are not updating the way you expect, first try "force reloading" your browser page – e.g. reload the page without using the cache.
Next, if problems persist, clear the sphinx cache (`rm -rf docs/_build`) and try building them again.

## Adding to the docs

Once you've added a new function to the documentation, you'll need to make sure there is a link somewhere in the documentation site pointing to it.
This should be added to `docs/api.md` under a relevant heading.

For tutorials and more in depth examples, consider adding a notebook to the [gcell-tutorials][] repository.

The tutorials are tied to this repository via a submodule.
To update the submodule, run `git submodule update --remote` from the root of the repository.
Subsequently, commit and push the changes in a PR.
This should be done before each release to ensure the tutorials are up to date.

[gcell-tutorials]: https://github.com/GET-Foundation/gcell-tutorials/

## docstrings format

We use the numpydoc style for writing docstrings.
We'd primarily suggest looking at existing docstrings for examples, but the [napolean guide to numpy style docstrings][] is also a great source.
If you're unfamiliar with the reStructuredText (rST) markup format, check out the [Sphinx rST primer][].

Some key points:

- We have some custom sphinx extensions activated. When in doubt, try to copy the style of existing docstrings.
- We autopopulate type information in docstrings when possible, so just add the type information to signatures.
- When docs exist in the same file as code, line length restrictions still apply. In files which are just docs, go with a sentence per line (for easier `git diff`s).
- Check that the docs look like what you expect them too! It's easy to forget to add a reference to function, be sure it got added and looks right.

Look at [sc.tl.louvain](https://github.com/GET-Foundation/gcell/blob/a811fee0ef44fcaecbde0cad6336336bce649484/gcell/tools/_louvain.py#L22-L90) as an example for everything mentioned here.

[napolean guide to numpy style docstrings]: https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html#example-numpy
[sphinx rst primer]: https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html

### `Params` section

The `Params` abbreviation is a legit replacement for `Parameters`.

To document parameter types use type annotations on function parameters.
These will automatically populate the docstrings on import, and when the documentation is built.

Use the python standard library types (defined in {mod}`collections.abc` and {mod}`typing` modules) for containers, e.g.
{class}`~collections.abc.Sequence`s (like `list`),
{class}`~collections.abc.Iterable`s (like `set`), and
{class}`~collections.abc.Mapping`s (like `dict`).
Always specify what these contain, e.g. `{'a': (1, 2)}` → `Mapping[str, Tuple[int, int]]`.
If you can’t use one of those, use a concrete class like `AnnData`.
If your parameter only accepts an enumeration of strings, specify them like so: `Literal['elem-1', 'elem-2']`.

### `Returns` section

There are three types of return sections – prose, tuple, and a mix of both.

1. Prose is for simple cases.
2. Tuple return sections are formatted like parameters. Other than in numpydoc, each tuple is first characterized by the identifier and *not* by its type. Provide type annotation in the function header.
3. Mix of prose and tuple is relevant in complicated cases, e.g. when you want to describe that you *added something as annotation to an \`AnnData\` object*.
