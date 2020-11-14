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

project(
  's-dftd3',
  'fortran',
  version: '0.2.0',
  license: 'LGPL-3.0-or-later',
  meson_version: '>=0.53',
  default_options: [
    'buildtype=debugoptimized',
    'default_library=both',
  ],
)
install = not (meson.is_subproject() and get_option('default_library') == 'static')

# General configuration information
inc_dirs = []
lib_deps = []
subdir('config')

# Collect source of the project
srcs = []
subdir('src')

# Library target
sdftd3_lib = library(
  meson.project_name(),
  sources: srcs,
  version: meson.project_version(),
  dependencies: lib_deps,
  include_directories: inc_dirs,
  install: install,
)

# Export dependency for other projects and test suite
sdftd3_inc = sdftd3_lib.private_dir_include()
sdftd3_dep = declare_dependency(
  link_with: sdftd3_lib,
  include_directories: sdftd3_inc,
  dependencies: lib_deps,
)

# Add applications
subdir('app')

# Package the license files
sdftd3_lic = files(
  'COPYING',
  'COPYING.LESSER',
)

if install
  # Distribute the license files in share/licenses/<name>
  install_data(
    sdftd3_lic,
    install_dir: get_option('datadir')/'licenses'/meson.project_name()
  )

  pkg = import('pkgconfig')
  pkg.generate(
    sdftd3_lib,
    description: 'Simple reimplementation of the DFT-D3 dispersion model',
  )
endif

# add the testsuite
subdir('test')
