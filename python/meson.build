# This file is part of s-dftd3.
# SPDX-Identifier: LGPL-3.0-or-later
#
# s-dftd3 is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# s-dftd3 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.

pymod = import('python')
python = pymod.find_installation(
  'python@0@'.format(get_option('python_version')),
  modules: [
    'cffi',
  ],
)
python_dep = python.dependency(required: true)

cc = meson.get_compiler('c')

# Python's CFFI is horrible in working with preprocessor statements,
# therefore, we have to preprocess the header before passing it to the ffibuilder
sdftd3_pp = configure_file(
  command: [cc, '-E', '@INPUT@'],
  input: sdftd3_header[0],
  output: '_libdftd3.h',
  capture: true,
)

# This is the actual out-of-line API processing of the ffibuilder
sdftd3_cffi_srcs = configure_file(
  command: [python, files('ffibuilder.py'), '@INPUT@', '@BASENAME@'],
  input: sdftd3_pp,
  output: '@BASENAME@.c',
)

# Actual generation of the Python extension, since the shared_module does not work
# well with dependency objects, we will trick it by linking a whole static lib
sdftd3_pyext = python.extension_module(
  '_libdftd3',
  link_whole: static_library(
    '_libdftd3',
    sdftd3_cffi_srcs,
    include_directories: [inc_dirs],
    dependencies: [sdftd3_dep, python_dep],
  ),
  dependencies: [sdftd3_dep, python_dep],
  subdir: 'dftd3',
  install: install,
)

if install
  python.install_sources(
    'dftd3'/'__init__.py',
    'dftd3'/'libdftd3.py',
    subdir: 'dftd3'
  )
endif
