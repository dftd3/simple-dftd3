/* This file is part of s-dftd3.
 * SPDX-Identifier: LGPL-3.0-or-later
 *
 * s-dftd3 is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * s-dftd3 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.
**/
#include <cstdio>

#include "dftd3.h"

int
test_cpp_compilation()
{
   printf("Start test: C++ compilation\n");
   
   // Test that the header can be included and basic API is accessible
   dftd3_error error = dftd3_new_error();
   
   // In C++, we cannot use the dftd3_delete macro (it's C-only)
   // Instead, we must explicitly call the delete functions
   dftd3_delete_error(&error);
   
   return 0;
}

int
main()
{
   return test_cpp_compilation();
}
