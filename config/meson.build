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

os = host_machine.system()
fc = meson.get_compiler('fortran')
fc_id = fc.get_id()

if fc_id == 'gcc'
  add_project_arguments(
    '-ffree-line-length-none',
    '-fbacktrace',
    language: 'fortran',
  )
elif fc_id == 'intel'
  add_project_arguments(
    '-traceback',
    language: 'fortran',
  )
elif fc_id == 'pgi'
  add_project_arguments(
    '-Mbackslash',
    '-Mallocatable=03',
    '-traceback',
    language: 'fortran',
  )
endif

if get_option('openmp')
  omp_dep = dependency('OpenMP')
  lib_deps += omp_dep
endif

if get_option('blas') != 'none'
  add_project_arguments('-DWITH_BLAS', language: 'fortran')
endif

if get_option('ilp64')
  ilp64 = true
  add_project_arguments('-DIK=i8', language:'fortran')
else
  ilp64 = false
  add_project_arguments('-DIK=i4', language:'fortran')
endif

if get_option('blas') == 'netlib'
  blas_dep = ilp64 ? dependency('blas64', required: false) : dependency('blas', required: false)
  if not blas_dep.found()
    blas_dep = ilp64 ? fc.find_library('blas64', required: true) : fc.find_library('blas', required: false)
  endif
  lib_deps += blas_dep
elif get_option('blas') == 'custom'
  blas_dep = []
  libs = get_option('blas_libs')
  if libs[0].startswith('-L')
    foreach lib: libs
      if lib != libs[0]
        blas_dep += fc.find_library(lib, dirs: libs[0].substring(2))
      endif
    endforeach
  else
    foreach lib: libs
      blas_dep += fc.find_library(lib)
    endforeach
  endif
  if (not fc.links('external dsytrs; call dgemm(); end', dependencies: [blas_dep,omp_dep]))
    error('Custom BLAS libraries do not link')
  endif
  lib_deps += blas_dep
endif

# Create the tool chain library as subproject
mctc_dep = dependency(
  'mctc-lib',
  version: '>=0.3.0',
  fallback: ['mctc-lib', 'mctc_dep'],
  default_options: ['default_library=static'],
)
lib_deps += mctc_dep

# Create the TOML Fortran library as subproject
tomlf_dep = dependency(
  'toml-f',
  fallback: ['toml-f', 'tomlf_dep'],
  default_options: ['default_library=static'],
)
exe_deps += tomlf_dep
