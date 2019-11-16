# This file is part of s-dftd3.
#
# Copyright (C) 2019 Sebastian Ehlert
#
# s-dftd3 is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# s-dftd3 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.

project('s-dftd3', 'fortran',
        version: '1.0',
        license: 'LGPL3',
        meson_version: '>=0.49')

conf = configuration_data()
conf.set('name', meson.project_name())
conf.set('version', meson.project_version())
conf.set('commit', run_command(find_program('git'),'show','-s','--format=%h').stdout().strip())
conf.set('date', run_command(find_program('date'),'-I').stdout().strip())
conf.set('author', run_command(find_program('whoami')).stdout().strip())
conf.set('origin', run_command(find_program('hostname')).stdout().strip())
configure_file(input : 'include/version.fh', output : 's-dftd3-version.fh',
  configuration : conf)

## source
ftest  = []
fsrcs  = []
fmain  = 'src/program_dftd.f90'

fsrcs += 'src/main_dftd.f90'
fsrcs += 'src/header.f90'
fsrcs += 'src/help.f90'
fsrcs += 'src/dftd3.f90'

fsrcs += 'src/typedef/environment.f90'
fsrcs += 'src/typedef/argparser.f90'
fsrcs += 'src/typedef/molecule.f90'
fsrcs += 'src/typedef/damping_parameters.f90'
fsrcs += 'src/typedef/options.f90'
fsrcs += 'src/typedef/results.f90'

fsrcs += 'src/io/molecule_reader.f90'

fsrcs += 'src/utils/file_utils.f90'
fsrcs += 'src/utils/periodic.f90'
fsrcs += 'src/utils/symbols.f90'

fsrcs += 'src/param/constants.f90'
fsrcs += 'src/param/units.f90'
fsrcs += 'src/param/atomic_units.f90'
fsrcs += 'src/param/pse.f90'
fsrcs += 'src/param/dftd3_parameters.f90'
fsrcs += 'src/param/r4r2_expectation_values.f90'
fsrcs += 'src/param/covalent_radii.f90'

incdir = include_directories('include')

dftd3_lib = library(meson.project_name(),
                    sources: fsrcs,
                    version: meson.project_version(),
                    include_directories: incdir)

dftd3_dep = declare_dependency(link_with: dftd3_lib)

dftd3_exe = executable(meson.project_name(),
                       sources: fmain,
                       link_with: dftd3_lib,
                       include_directories: incdir)

# test executable
test_env = environment()
test_env.prepend('LD_LIBRARY_PATH', meson.current_build_dir())

test('print version', dftd3_exe, env: test_env, args: '--version')
test('print help', dftd3_exe, env: test_env, args: '--help')
test('no arguments', dftd3_exe, env: test_env, should_fail: true)
test('taxol.xyz', dftd3_exe, env: test_env,
     args: meson.current_source_dir()/'assets'/'inputs'/'xyz'/'taxol.xyz')
