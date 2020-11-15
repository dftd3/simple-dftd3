#!/usr/bin/env python3
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

import cffi, sys

if len(sys.argv) != 3:
    raise RuntimeError("Requires two arguments")

header_file = sys.argv[1]
module_name = sys.argv[2]

ffibuilder = cffi.FFI()
ffibuilder.set_source(
    module_name,
    '#include "dftd3.h"',
)

with open(header_file) as f:
    ffibuilder.cdef(f.read())

if __name__ == '__main__':
    ffibuilder.distutils_extension('.')
