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
 */

/* Check that the C API reproduces the serial result when distributed over
 * MPI ranks, this driver is only built with MPI support and run under mpiexec */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "s-dftd3.h"

static const double thr = 1.0e-9;

static const int natoms = 7;
static const int attyp[7] = {6, 6, 6, 1, 1, 1, 1};
static const double coord[21] = {
    0.00000000000000, 0.00000000000000, -1.79755622305860,
    0.00000000000000, 0.00000000000000, 0.95338756106749,
    0.00000000000000, 0.00000000000000, 3.22281255790261,
    -0.96412815539807, -1.66991895015711, -2.53624948351102,
    -0.96412815539807, 1.66991895015711, -2.53624948351102,
    1.92825631079613, 0.00000000000000, -2.53624948351102,
    0.00000000000000, 0.00000000000000, 5.23010455462158,
};

static int check(int cond, const char *label) {
  if (!cond) {
    printf("[Fatal] Distributed %s differs from serial result\n", label);
    return 1;
  }
  return 0;
}

static int run(dftd3_error error, dftd3_structure mol, MPI_Comm comm,
               double *energy, double *gradient, double *sigma) {
  dftd3_model disp = dftd3_new_d3_model(error, mol);
  if (dftd3_check_error(error)) return 1;

  if (comm != MPI_COMM_NULL) {
    dftd3_set_model_mpi_comm(error, disp, MPI_Comm_c2f(comm));
    if (dftd3_check_error(error)) {
      dftd3_delete(disp);
      return 1;
    }
  }

  dftd3_param param = dftd3_load_rational_damping(error, "pbe0", true);
  if (dftd3_check_error(error)) {
    dftd3_delete(disp);
    return 1;
  }

  dftd3_get_dispersion(error, mol, disp, param, energy, gradient, sigma);
  dftd3_delete(param);
  dftd3_delete(disp);
  return dftd3_check_error(error) ? 1 : 0;
}

int main(void) {
  int rank, stat = 0;
  double energy, energy_ref;
  double gradient[3 * 7], gradient_ref[3 * 7];
  double sigma[9], sigma_ref[9];
  int i;

  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (!dftd3_has_feature("mpi")) {
    printf("[Fatal] Library reports MPI support as unavailable\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  dftd3_error error = dftd3_new_error();
  dftd3_structure mol =
      dftd3_new_structure(error, natoms, attyp, coord, NULL, NULL);
  if (dftd3_check_error(error)) MPI_Abort(MPI_COMM_WORLD, 1);

  if (run(error, mol, MPI_COMM_NULL, &energy_ref, gradient_ref, sigma_ref))
    MPI_Abort(MPI_COMM_WORLD, 1);
  if (run(error, mol, MPI_COMM_WORLD, &energy, gradient, sigma))
    MPI_Abort(MPI_COMM_WORLD, 1);

  stat += check(fabs(energy - energy_ref) < thr, "energy");
  for (i = 0; i < 3 * natoms; i++)
    stat += check(fabs(gradient[i] - gradient_ref[i]) < thr, "gradient");
  for (i = 0; i < 9; i++)
    stat += check(fabs(sigma[i] - sigma_ref[i]) < thr, "virial");

  if (stat) MPI_Abort(MPI_COMM_WORLD, 1);
  if (rank == 0) printf("# distributed C API result matches serial result\n");

  dftd3_delete(mol);
  dftd3_delete(error);
  MPI_Finalize();
  return 0;
}
