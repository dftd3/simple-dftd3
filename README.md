# The D3 dispersion model

[![Latest Version](https://img.shields.io/github/v/release/awvwgk/simple-dftd3)](https://github.com/awvwgk/simple-dftd3/releases/latest)
[![LGPL-3.0-or-later](https://img.shields.io/github/license/awvwgk/simple-dftd3)](COPYING)
[![CI](https://github.com/awvwgk/simple-dftd3/workflows/CI/badge.svg)](https://github.com/awvwgk/simple-dftd3/actions)
[![docs](https://github.com/awvwgk/simple-dftd3/workflows/docs/badge.svg)](https://awvwgk.github.io/simple-dftd3/)

A simple drop-in replacement for ``dftd3``.

This program provides a small and easy to use implementation of the DFT-D3
dispersion correction
(see [*JCP* **132**, 154104 (2010)](https://dx.doi.org/10.1063/1.3382344)
and [*JCC* **32**, 1456 (2011)](https://dx.doi.org/10.1002/jcc.21759) for details).

It is mostly based on the [`dftd4`](https://github.com/dftd4/dftd4) program and
borrows one or two ideas from the implementation in [`ased3`](https://github.com/ehermes/ased3).


## Installation


### Conda package

[![Conda Version](https://img.shields.io/conda/vn/conda-forge/simple-dftd3.svg?label=simple-dftd3)](https://anaconda.org/conda-forge/simple-dftd3)

This project is packaged for the *conda* package manager and available on the *conda-forge* channel.
To install the *conda* package manager we recommend the [miniforge](https://github.com/conda-forge/miniforge/releases) installer.
If the *conda-forge* channel is not yet enabled, add it to your channels with

```
conda config --add channels conda-forge
```

Once the *conda-forge* channel has been enabled, this project can be installed with:

```
conda install simple-dftd3
```

It is possible to list all of the versions available on your platform with:

```
conda search simple-dftd3 --channel conda-forge
```

Now you are ready to use ``s-dftd3``.


### Building from Source

To build this project from the source code in this repository you need to have
- a Fortran compiler supporting Fortran 2008
- [meson](https://mesonbuild.com) version 0.53 or newer
- a build-system backend, *i.e.* [ninja](https://ninja-build.org) version 1.7 or newer

Optional dependencies are
- BLAS (enabled with `-Dblas=netlib`)
- asciidoctor to build the manual page
- FORD to build the developer documentation

Setup a build with

```
meson setup _build
```

You can select the Fortran compiler by the `FC` environment variable, this project is currently tested with GCC 9 on Ubuntu, MacOS and Windows.
To compile and run the projects testsuite use

```
meson test -C _build --print-errorlogs
```

If the testsuite passes you can install with

```
meson configure _build --prefix=/path/to/install
meson install -C _build
```

This might require administrator access depending on the chosen install prefix.


## Usage

DFT-D3 calculations can be performed with the ``s-dftd3`` executable.
To calculate the dispersion correction for PBE0-D3(BJ)-ATM run:

```
s-dftd3 --bj pbe0 --atm coord
```

In case you want to access the DFT-D3 results from other programs, dump the results to JSON with
(the ``--noedisp`` flag prevents the ``.EDISP`` file generation):

```
s-dftd3 --bj pbe0 --atm --json --noedisp --grad struct.xyz
```

Dispersion related properties can be calculated as well:

```
s-dftd3 --property geo.gen
```

For an overview over all command line arguments use the ``--help`` argument or checkout the [``s-dftd3(1)``](man/s-dftd3.1.adoc) manpage.


## License

This project is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This project is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the
GNU Lesser General Public License for more details.

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in this project by you, as defined in the
GNU Lesser General Public license, shall be licensed as above, without any
additional terms or conditions.
