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

# Include the mctc-lib project
message(STATUS "Looking for mctc-lib...")
if(NOT TARGET mctc-lib)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mctc-lib/CMakeLists.txt)
    add_subdirectory("mctc-lib")
  else()
    set("mctc-lib-url" "https://github.com/grimme-lab/mctc-lib")
    message(STATUS "Retrieving mctc-lib from ${mctc-lib-url}")
    include(FetchContent)
    FetchContent_Declare(
      "mctc-lib"
      GIT_REPOSITORY "${mctc-lib-url}"
      GIT_TAG "HEAD"
    )
    FetchContent_MakeAvailable("mctc-lib")
  endif()
else()
  message(STATUS "mctc-lib already present")
endif()

list(
  APPEND lib-deps
  "mctc-lib"
)

# Include the mstore project
message(STATUS "Looking for mstore...")
if(NOT TARGET mstore)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mstore/CMakeLists.txt)
    add_subdirectory("mstore")
  else()
    set("mstore-url" "https://github.com/grimme-lab/mstore")
    message(STATUS "Retrieving mstore from ${mstore-url}")
    include(FetchContent)
    FetchContent_Declare(
      "mstore"
      GIT_REPOSITORY "${mstore-url}"
      GIT_TAG "HEAD"
    )
    FetchContent_MakeAvailable("mstore")
  endif()
endif()

set(lib-deps "${lib-deps}" PARENT_SCOPE)
