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

mpiexec = find_program('mpiexec', required: false)
if not mpiexec.found()
  subdir_done()
endif

mpi_tester = executable(
  'mpi-tester',
  sources: files('main.f90'),
  dependencies: [sdftd3_dep, mstore_dep, mpi_dep],
)

foreach n : ['1', '2', '3']
  test(
    'ranks-' + n,
    mpiexec,
    args: ['-n', n, mpi_tester],
    suite: 'mpi',
    # one thread per rank, oversubscribing the machine makes the ranks crawl
    env: {'OMP_NUM_THREADS': '1'},
    timeout: 120,
  )
endforeach
