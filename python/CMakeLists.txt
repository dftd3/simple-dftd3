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

cmake_minimum_required(VERSION 3.18)

set(_sdftd3_python_standalone FALSE)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(_sdftd3_python_standalone TRUE)
  project(
    dftd3-python
    LANGUAGES C
    VERSION 1.6.0
    DESCRIPTION "Python API of the DFT-D3 project"
  )
endif()

include(GNUInstallDirs)

find_package(Python3 3.8 REQUIRED COMPONENTS "Interpreter" "Development.Module")
execute_process(
  COMMAND "${Python3_EXECUTABLE}" -c "import cffi, setuptools"
  RESULT_VARIABLE _sdftd3_cffi_status
  ERROR_QUIET
)
if(NOT _sdftd3_cffi_status EQUAL 0)
  message(FATAL_ERROR
    "The Python interpreter does not provide the cffi and setuptools modules"
  )
endif()

set(_sdftd3_target "")
set(_sdftd3_include_dir "")

# An in-tree build already provides the native target. A standalone Python
# build first tries an installed CMake package and then the source bundled in
# the Python sdist.
if(TARGET s-dftd3::s-dftd3)
  set(_sdftd3_target "s-dftd3::s-dftd3")
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../include/dftd3.h")
    set(_sdftd3_include_dir "${CMAKE_CURRENT_SOURCE_DIR}/../include")
  endif()
elseif(TARGET s-dftd3::s-dftd3-lib)
  set(_sdftd3_target "s-dftd3::s-dftd3-lib")
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../include/dftd3.h")
    set(_sdftd3_include_dir "${CMAKE_CURRENT_SOURCE_DIR}/../include")
  endif()
elseif(_sdftd3_python_standalone)
  enable_language(Fortran)
  find_package(s-dftd3 "${PROJECT_VERSION}" CONFIG QUIET)
  if(TARGET s-dftd3::s-dftd3)
    set(_sdftd3_target "s-dftd3::s-dftd3")
  elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/subprojects/s-dftd3/CMakeLists.txt")
    set(BUILD_SHARED_LIBS FALSE CACHE BOOL "Build shared libraries" FORCE)
    set(BUILD_TESTING FALSE CACHE BOOL "Build tests in bundled dependencies" FORCE)
    set(SDFTD3_WITH_API TRUE CACHE BOOL "Build the DFT-D3 C API" FORCE)
    set(SDFTD3_WITH_PYTHON FALSE CACHE BOOL "Build the DFT-D3 Python API" FORCE)
    set(SDFTD3_WITH_TESTS FALSE CACHE BOOL "Build the DFT-D3 tests" FORCE)
    add_subdirectory("subprojects/s-dftd3" "subprojects/s-dftd3")
    set(_sdftd3_target "s-dftd3::s-dftd3")
    set(_sdftd3_include_dir "${CMAKE_CURRENT_SOURCE_DIR}/subprojects/s-dftd3/include")
  endif()
endif()

if(NOT _sdftd3_target)
  message(FATAL_ERROR
    "Could not find s-dftd3 >= ${PROJECT_VERSION}; install its CMake package "
    "or build from a Python sdist containing subprojects/s-dftd3"
  )
endif()
if(DEFINED SDFTD3_WITH_API AND NOT SDFTD3_WITH_API)
  message(FATAL_ERROR "The dftd3 Python extension requires the s-dftd3 C API")
endif()

if(NOT _sdftd3_include_dir)
  get_target_property(
    _sdftd3_interface_include_dirs
    "${_sdftd3_target}"
    INTERFACE_INCLUDE_DIRECTORIES
  )
  unset(_sdftd3_include_dir)
  unset(_sdftd3_include_dir CACHE)
  find_path(
    _sdftd3_include_dir
    NAMES "dftd3.h"
    HINTS ${_sdftd3_interface_include_dirs}
  )
  if(NOT _sdftd3_include_dir)
    message(FATAL_ERROR "Could not locate dftd3.h from ${_sdftd3_target}")
  endif()
endif()

set(_sdftd3_header "${CMAKE_CURRENT_SOURCE_DIR}/include/_dftd3.h")
set(_sdftd3_preprocessed_header "${CMAKE_CURRENT_BINARY_DIR}/_libdftd3.h")
set(_sdftd3_cffi_source "${CMAKE_CURRENT_BINARY_DIR}/_libdftd3.c")

# CFFI cannot parse all preprocessor constructs in the public C header.
add_custom_command(
  OUTPUT "${_sdftd3_preprocessed_header}"
  COMMAND
    "${CMAKE_C_COMPILER}"
    "-I${_sdftd3_include_dir}"
    "-DSDFTD3_CFFI"
    "-E"
    "${_sdftd3_header}"
    "-o"
    "${_sdftd3_preprocessed_header}"
  DEPENDS
    "${_sdftd3_header}"
    "${_sdftd3_include_dir}/dftd3.h"
  VERBATIM
)

add_custom_command(
  OUTPUT "${_sdftd3_cffi_source}"
  COMMAND
    "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ffibuilder.py"
    "${_sdftd3_preprocessed_header}"
    "_libdftd3"
  DEPENDS
    "${_sdftd3_preprocessed_header}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ffibuilder.py"
  VERBATIM
)

Python3_add_library("_libdftd3" MODULE WITH_SOABI "${_sdftd3_cffi_source}")
target_include_directories("_libdftd3" PRIVATE "${_sdftd3_include_dir}")
target_link_libraries("_libdftd3" PRIVATE "${_sdftd3_target}")
if(CMAKE_Fortran_COMPILER_LOADED)
  set_property(TARGET "_libdftd3" PROPERTY LINKER_LANGUAGE "Fortran")
endif()

execute_process(
  COMMAND
    "${Python3_EXECUTABLE}"
    -c
    "import sysconfig, sys; print(sysconfig.get_path('platlib', vars={'base': '', 'platbase': '', 'platlibdir': sys.argv[1]}).lstrip('/'))"
    "${CMAKE_INSTALL_LIBDIR}"
  OUTPUT_VARIABLE _sdftd3_python_install_root
  OUTPUT_STRIP_TRAILING_WHITESPACE
  RESULT_VARIABLE _sdftd3_sysconfig_status
)
if(NOT _sdftd3_sysconfig_status EQUAL 0)
  message(FATAL_ERROR "Could not determine the Python platform-library directory")
endif()
set(
  SDFTD3_PYTHON_INSTALL_DIR
  "${_sdftd3_python_install_root}/dftd3"
  CACHE PATH
  "Directory for the dftd3 Python package"
)

install(
  TARGETS "_libdftd3"
  LIBRARY DESTINATION "${SDFTD3_PYTHON_INSTALL_DIR}"
  COMPONENT "python"
  RUNTIME DESTINATION "${SDFTD3_PYTHON_INSTALL_DIR}"
  COMPONENT "python"
)

set(
  _sdftd3_python_sources
  "dftd3/__init__.py"
  "dftd3/ase.py"
  "dftd3/interface.py"
  "dftd3/library.py"
  "dftd3/parameters.py"
  "dftd3/pyscf.py"
  "dftd3/qcschema.py"
  "dftd3/test_ase.py"
  "dftd3/test_interface.py"
  "dftd3/test_library.py"
  "dftd3/test_parameters.py"
  "dftd3/test_pyscf.py"
  "dftd3/test_qcschema.py"
  "dftd3/py.typed"
)
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/dftd3/parameters.toml")
  list(APPEND _sdftd3_python_sources "dftd3/parameters.toml")
endif()

install(
  FILES ${_sdftd3_python_sources}
  DESTINATION "${SDFTD3_PYTHON_INSTALL_DIR}"
  COMPONENT "python"
)
