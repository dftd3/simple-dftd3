# The D3 dispersion model

A simple drop-in replacement for ``dftd3``.

This program provides a small and easy to use implementation of the DFT-D3
dispersion correction
(see [*JCP* **132**, 154104 (2010)](https://dx.doi.org/10.1063/1.3382344)
and [*JCC* **32**, 1456 (2011)](https://dx.doi.org/10.1002/jcc.21759) for details).

It is mostly based on the [`dftd4`](https://github.com/dftd4/dftd4) program and
borrows one or two ideas from the implementation in [`ased3`](https://github.com/ehermes/ased3).

## Installation

To build this project from the source code in this repository you need to have
- a Fortran compiler supporting Fortran 2008
- [meson](https://mesonbuild.com) version 0.53 or newer
- a build-system backend, *i.e.* [ninja](https://ninja-build.org) version 1.7 or newer

Setup a build with

```
meson setup _build
```

You can select the Fortran compiler by the `FC` environment variable, currently this project supports GCC and Intel compilers.
To compile the project run

```
meson compile -C _build
```

You can run the projects testsuite with

```
meson test -C _build --print-errorlogs
```

If the testsuite passes you can install with

```
meson install -C _build
```

This might require administrator access.
You can alter the install prefix by ``meson configure _build --prefix=/path/to/install``.


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
