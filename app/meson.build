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

sdftd3_exe = executable(
  meson.project_name(),
  sources: files('main.f90'),
  dependencies: sdftd3_dep,
  install: install,
)

if os != 'windows'
  test('app-version', sdftd3_exe, args: '--version')
  test('app-help', sdftd3_exe, args: '--help')
  test('app-noargs', sdftd3_exe, should_fail: true)
  test('app-example', sdftd3_exe, args: files('example.xyz'))
  test('app-zero', sdftd3_exe, args: ['--json', '--zero-param', '1.0', '1.868', '1.613', files('example.xyz')])
  test('app-rational', sdftd3_exe, args: ['--bj-param', '1.0', '1.2576', '0.3768', '4.5865', files('example.xyz')])
endif
