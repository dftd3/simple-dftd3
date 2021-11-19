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

app_tester = find_program(meson.current_source_dir() / 'tester.py')

test(
  'energy',
  app_tester,
  args: [
    sdftd3_exe,
    files('01-energy.json'),
    files('01-ammonia.tmol'),
    '--noedisp',
    '--zero', 'tpssh',
  ],
  suite: 'app',
)
test(
  'gradient',
  app_tester,
  args: [
    sdftd3_exe,
    files('02-gradient.json'),
    files('02-nitralin.mol'),
    '--noedisp',
    '--bj', 'scan',
    '--grad',
  ],
  suite: 'app',
)
