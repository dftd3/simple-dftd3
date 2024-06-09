# The D3 dispersion model

[![Latest Version](https://img.shields.io/github/v/release/dftd3/simple-dftd3)](https://github.com/dftd3/simple-dftd3/releases/latest)
[![LGPL-3.0-or-later](https://img.shields.io/github/license/dftd3/simple-dftd3)](COPYING)
[![CI](https://github.com/dftd3/simple-dftd3/workflows/CI/badge.svg)](https://github.com/dftd3/simple-dftd3/actions)
[![Documentation](https://readthedocs.org/projects/dftd3/badge/?version=latest)](https://dftd3.readthedocs.io/en/latest/)
[![docs](https://github.com/dftd3/simple-dftd3/actions/workflows/docs.yml/badge.svg)](https://dftd3.github.io/simple-dftd3/)
[![codecov](https://codecov.io/gh/dftd3/simple-dftd3/branch/main/graph/badge.svg)](https://codecov.io/gh/dftd3/simple-dftd3)

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
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/dftd3-python.svg?label=dftd3-python)](https://anaconda.org/conda-forge/dftd3-python)

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

If you want to enable the Python API as well install

```
conda install dftd3-python
```

It is possible to list all of the versions available on your platform with:

```
conda search simple-dftd3 --channel conda-forge
```

Now you are ready to use ``s-dftd3``.


### Building from Source

To build this project from the source code in this repository you need to have
a Fortran compiler supporting Fortran 2008 and one of the supported build systems:
- [meson](https://mesonbuild.com) version 0.55 or newer, with
  a build-system backend, *i.e.* [ninja](https://ninja-build.org) version 1.7 or newer
- [cmake](https://cmake.org) version 3.14 or newer, with
  a build-system backend, *i.e.* [ninja](https://ninja-build.org) version 1.10 or newer
- [fpm](https://github.com/fortran-lang/fpm) version 0.3.0 or newer

Meson is the primary build system and provides feature-complete functionality of this project.
CMake and fpm support is available but the functionality of the project is limited.
This project is currently tested with GCC 9 on Ubuntu, MacOS and Windows as well as Intel Fortran on Ubuntu.


#### Building with meson

Optional dependencies are
- asciidoctor to build the manual page
- FORD to build the developer documentation
- Python 3.6 or newer with the CFFI package installed to build the Python API

Setup a build with

```
meson setup _build
```

You can select the Fortran compiler by the `FC` environment variable.
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

To include ``s-dftd3`` in your project add the following wrap file to your subprojects directory:

```ini
[wrap-git]
directory = s-dftd3
url = https://github.com/dftd3/simple-dftd3
revision = head
```

You can retrieve the dependency from the wrap fallback with

```meson
sdftd3_dep = dependency('s-dftd3', fallback: ['s-dftd3', 'sdftd3_dep'])
```

and add it as dependency to your targets.


#### Building with CMake

Alternatively, this project can be build with CMake (in this case ninja 1.10 or newer is required):

```
cmake -B _build -G Ninja
```

To compile the project with CMake run

```
cmake --build _build
```

You can run the project testsuite with

```
pushd _build && ctest && popd
```

To include ``s-dftd3`` in your CMake project retrieve it using the ``FetchContent`` module:

```cmake
if(NOT TARGET s-dftd3)
  set("s-dftd3-url" "https://github.com/dftd3/simple-dftd3")
  message(STATUS "Retrieving s-dftd3 from ${s-dftd3-url}")
  include(FetchContent)
  FetchContent_Declare(
    "s-dftd3"
    GIT_REPOSITORY "${s-dftd3-url}"
    GIT_TAG "HEAD"
  )
  FetchContent_MakeAvailable("s-dftd3")
endif()
```

And link against the ``"s-dftd3"`` interface library.

```cmake
target_link_libraries("${PROJECT_NAME}-lib" PUBLIC "s-dftd3")
```


#### Building with fpm

Invoke fpm in the project root with

```
fpm build
```

To run the testsuite use

```
fpm test
```

You can access the ``s-dftd3`` program using the run subcommand

```
fpm run -- --help
```

To use ``s-dftd3`` for testing include it as dependency in your package manifest

```toml
[dependencies]
s-dftd3.git = "https://github.com/dftd3/simple-dftd3"
```


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


## API access

This DFT-D3 implementation provides first class API support Fortran, C and Python.
Other programming languages should try to interface via one of those three APIs.
To provide first class API support for a new language the interface specification should be available as meson build files.


### Fortran API

The recommended way to access the Fortran module API is by using ``dftd3`` as a meson subproject.
Alternatively, the project is accessible by the Fortran package manager ([fpm](https://github.com/fortran-lang/fpm)) or as CMake subproject as explained above.

The complete API is available from ``dftd3`` module, the individual modules are available to the user as well but are not part of the public API and therefore not guaranteed to remain stable.
API compatibility is only guaranteed for the same minor version, while ABI compatibility cannot be guaranteed in a pre 1.0 stage.

The communication with the Fortran API uses the ``error_type`` and ``structure_type`` of the modular computation tool chain library (mctc-lib) to handle errors and represent geometries, respectively.


### C API

The C API provides access to the basic Fortran objects and their most important methods to interact with them.
All Fortran objects are available as opaque ``void*`` in C and can only be manipulated with the correct API calls.
To evaluate a dispersion correction in C four objects are available:

1. the error handler:

   Simple error handler to carry runtime exceptions created by the library.
   Exceptions can be handled and/or transfered to the downstream error handling system by this means.

2. the molecular structure data:

   Provides a representation of the molecular structure with immutable number of atoms, atomic species, total charge and boundary conditions.
   The object provides a way to update coordinates and lattice parameters, to update immutable quantities the object has to be recreated.

3. the dispersion model:

   Instantiated for a given molecular structure type, it carries no information on the geometry but relies on the atomic species of the structure object.
   Recreating a structure object requires to recreate the dispersion model as well.

4. the damping parameters:

   Damping parameter object determining the short-range behaviour of the dispersion correction.
   Standard damping parameters like the rational damping are independent of the molecular structure and can easily be reused for several structures or easily exchanged.

The user is responsible for creating and deleting the objects to avoid memory leaks.


### Python API

The Python API is disabled by default and can be built in-tree or out-of-tree.
The in-tree build is mainly meant for end users and packages.
To build the Python API with the normal project set the ``python`` option in the configuration step with

```sh
meson setup _build -Dpython=true -Dpython_version=$(which python3)
```

The Python version can be used to select a different Python version, it defaults to `'python3'`.
Python 2 is not supported with this project, the Python version key is meant to select between several local Python 3 versions.

Proceed with the build as described before and install the projects to make the Python API available in the selected prefix.

For the out-of-tree build see the instructions in the [``python``](./python) directory.


## Contributing

This is a volunteer open source projects and contributions are always welcome. 
Please, take a moment to read the [contributing guidelines](CONTRIBUTING.md).


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
