# Installation

To use `gcell` from another project, install it using your favourite environment manager:

```console
$ pip install git+https://github.com/GET-Foundation/gcell.git@main
```

(dev-install-instructions)=

## Development Version

To work with the latest version [on GitHub][]: clone the repository and `cd` into its root directory.

```console
$ gh repo clone GET-Foundation/gcell
$ cd gcell
```

::::{tabs}

:::{group-tab} uv (recommended)
```console
$ uv sync --all-extras        # install all dependencies
$ uv run pytest               # run tests
$ uv run sphinx-build docs docs/_build  # build docs
```
:::

:::{group-tab} Pip/PyPI
If you are using `pip>=21.3`, an editable install can be made:

```console
$ python -m venv .venv
$ source .venv/bin/activate
$ pip install -e '.[dev,test]'
```
:::

:::{group-tab} Conda
If you want to let `conda` handle the installations of dependencies, do:

```console
$ pipx install beni
$ beni pyproject.toml > environment.yml
$ conda env create -f environment.yml
$ conda activate gcell
$ pip install -e '.[dev,doc,test]'
```

For instructions on how to work with the code, see the {ref}`contribution guide <contribution-guide>`.
:::

::::

[on github]: https://github.com/GET-Foundation/gcell

## Docker

If you're using [Docker][], you can use e.g. the image [gcfntnu/gcell][] from Docker Hub.

[docker]: https://en.wikipedia.org/wiki/Docker_(software)
[gcfntnu/gcell]: https://hub.docker.com/r/gcfntnu/gcell

## Troubleshooting

If you get a `Permission denied` error, never use `sudo pip`. Instead, use virtual environments or:

```console
$ pip install --user gcell
```
