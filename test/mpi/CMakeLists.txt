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

add_executable(
  "${PROJECT_NAME}-mpi-tester"
  "main.f90"
)
target_link_libraries(
  "${PROJECT_NAME}-mpi-tester"
  PRIVATE
  "${PROJECT_NAME}-lib"
  "mstore::mstore"
  "MPI::MPI_Fortran"
)

foreach(n RANGE 1 3)
  add_test(
    NAME "${PROJECT_NAME}/mpi-ranks-${n}"
    COMMAND
    "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${n}"
    "$<TARGET_FILE:${PROJECT_NAME}-mpi-tester>"
  )
  # one thread per rank, oversubscribing the machine makes the ranks crawl
  set_tests_properties(
    "${PROJECT_NAME}/mpi-ranks-${n}"
    PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=1"
    TIMEOUT 120
  )
endforeach()
