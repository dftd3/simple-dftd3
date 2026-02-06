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

api_tester = executable(
  'api-tester',
  sources: files('api-test.c'),
  dependencies: sdftd3_dep,
)

test('api-test', api_tester, suite: 'api')

# C++ compatibility test
has_cpp = add_languages('cpp', required: false)
if has_cpp
  cpp_tester = executable(
    'cpp-tester',
    sources: files('cpp-test.cpp'),
    dependencies: sdftd3_dep,
  )

  test('cpp-test', cpp_tester, suite: 'api')
endif
