# Reimplementation of the DFT-D3 program

This program provides a small and easy to use implementation of the DFT-D3
dispersion correction
(see [*JCP* **132**, 154104 (2010)](https://dx.doi.org/10.1063/1.3382344)
and [*JCC* **32**, 1456 (2011)](https://dx.doi.org/10.1002/jcc.21759) for details).

It is mostly based on the [`dftd4`](https://github.com/dftd4/dftd4) program and
borrows one or two ideas from the implementation in [`ased3`](https://github.com/ehermes/ased3).

## Installing

To build this project we use [`meson`](https://mesonbuild.com/) and require
a fairly new version like 0.49 or newer for a decent Fortran support.
For the default backend [`ninja`](https://ninja-build.org/) version 1.5 or newer
has to be provided.

To perform a build run

```bash
meson setup build_gcc
ninja -C build_gcc test
```

You can install the program locally by running install with `ninja`

```bash
[sudo] ninja -C build_gcc install
```

## License

`s-dftd3` is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

`s-dftd3` is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the
GNU General Public License for more details.
