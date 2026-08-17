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
#include <math.h>
#include <assert.h>
#include <stdio.h>

#include "dftd3.h"

static inline void
show_error(dftd3_error error)
{
   char message[512];
   dftd3_get_error(error, message, NULL);
   printf("[Message] %s\n", message);
}

static inline dftd3_structure
get_test_structure(dftd3_error error)
{
   int const natoms = 7;
   int const attyp[7] = {6,6,6,1,1,1,1};
   double const coord[21] =
      {0.00000000000000, 0.00000000000000,-1.79755622305860,
       0.00000000000000, 0.00000000000000, 0.95338756106749,
       0.00000000000000, 0.00000000000000, 3.22281255790261,
      -0.96412815539807,-1.66991895015711,-2.53624948351102,
      -0.96412815539807, 1.66991895015711,-2.53624948351102,
       1.92825631079613, 0.00000000000000,-2.53624948351102,
       0.00000000000000, 0.00000000000000, 5.23010455462158};
   return dftd3_new_structure(error, natoms, attyp, coord, NULL, NULL);
}

static inline dftd3_structure
get_periodic_test_structure(dftd3_error error)
{
   int const natoms = 2;
   int const attyp[2] = {6,8};
   double const coord[6] =
      {0.40000000000000, 0.80000000000000, 1.20000000000000,
       4.40000000000000, 2.80000000000000, 3.60000000000000};
   double const lattice[9] =
      {8.00000000000000, 0.00000000000000, 0.00000000000000,
       0.00000000000000, 8.00000000000000, 0.00000000000000,
       0.00000000000000, 0.00000000000000, 8.00000000000000};
   bool const periodic[3] = {true, true, true};
   return dftd3_new_structure(error, natoms, attyp, coord, lattice, periodic);
}

int
test_version (void)
{
   printf("Start test: version\n");
   return dftd3_get_version() > 0 ? 0 : 1;
}

int
test_uninitialized_error (void)
{
   printf("Start test: uninitialized error\n");
   dftd3_error error = NULL;
   return dftd3_check_error(error) ? 0 : 1;
}

int
test_uninitialized_structure (void)
{
   printf("Start test: uninitialized structure\n");
   dftd3_error error = NULL;
   dftd3_structure mol = NULL;

   error = dftd3_new_error();

   double xyz[6] = {0.0};
   dftd3_update_structure(error, mol, xyz, NULL);
   if (!dftd3_check_error(error)) goto unexpected;

   show_error(error);

   dftd3_delete(error);
   return 0;

unexpected:
   printf("[Fatal] Unexpected pass for uninitialized-structure test\n");
   dftd3_delete(error);
   return 1;
}

int
test_uninitialized_param (void)
{
   printf("Start test: uninitialized parameters\n");
   dftd3_error error = NULL;
   dftd3_structure mol = NULL;
   dftd3_model disp = NULL;
   dftd3_param param = NULL;
   double energy;

   error = dftd3_new_error();

   dftd3_get_dispersion(error, mol, disp, param, &energy, NULL, NULL);
   if (!dftd3_check_error(error)) goto unexpected;

   show_error(error);

   dftd3_get_pairwise_dispersion(error, mol, disp, param, NULL, NULL);
   if (!dftd3_check_error(error)) goto unexpected;

   show_error(error);

   mol = get_test_structure(error);
   dftd3_get_dispersion(error, mol, disp, param, &energy, NULL, NULL);
   if (!dftd3_check_error(error)) goto unexpected;

   show_error(error);

   dftd3_get_pairwise_dispersion(error, mol, disp, param, NULL, NULL);
   if (!dftd3_check_error(error)) goto unexpected;

   show_error(error);

   disp = dftd3_new_d3_model(error, mol);

   dftd3_get_dispersion(error, mol, disp, param, &energy, NULL, NULL);
   if (!dftd3_check_error(error)) goto unexpected;

   show_error(error);

   dftd3_get_pairwise_dispersion(error, mol, disp, param, NULL, NULL);
   if (!dftd3_check_error(error)) goto unexpected;

   show_error(error);

   dftd3_get_dispersion_hessian(error, mol, disp, param, &energy, NULL);
   if (!dftd3_check_error(error)) goto unexpected;

   show_error(error);

   dftd3_delete(error);
   dftd3_delete(disp);
   dftd3_delete(mol);
   return 0;

unexpected:
   printf("[Fatal] Unexpected pass for uninitialized-parameters test\n");
   dftd3_delete(error);
   dftd3_delete(disp);
   dftd3_delete(mol);
   return 1;
}

int
test_uninitialized_model (void)
{
   printf("Start test: uninitialized model\n");
   dftd3_error error = NULL;
   dftd3_model disp = NULL;

   error = dftd3_new_error();

   dftd3_set_model_realspace_cutoff(error, disp, 0.0, 0.0, 0.0);
   if (!dftd3_check_error(error)) goto unexpected;

   show_error(error);

   dftd3_delete(error);
   dftd3_delete(disp);
   return 0;

unexpected:
   printf("[Fatal] Unexpected pass for uninitialized-model test\n");
   dftd3_delete(error);
   dftd3_delete(disp);
   return 1;
}

int
test_uninitialized_gcp (void)
{
   printf("Start test: uninitialized counter-poise\n");
   dftd3_error error = NULL;
   dftd3_structure mol = NULL;
   dftd3_gcp gcp = NULL;
   double energy;

   error = dftd3_new_error();

   dftd3_set_gcp_realspace_cutoff(error, gcp, 0.0, 0.0);
   if (!dftd3_check_error(error)) goto unexpected;

   show_error(error);

   dftd3_get_counterpoise(error, mol, gcp, &energy, NULL, NULL);
   if (!dftd3_check_error(error)) goto unexpected;

   show_error(error);

   mol = get_test_structure(error);

   dftd3_get_counterpoise(error, mol, gcp, &energy, NULL, NULL);
   if (!dftd3_check_error(error)) goto unexpected;

   show_error(error);

   dftd3_get_counterpoise_hessian(error, mol, gcp, &energy, NULL);
   if (!dftd3_check_error(error)) goto unexpected;

   show_error(error);

   dftd3_delete(error);
   dftd3_delete(mol);
   return 0;

unexpected:
   printf("[Fatal] Unexpected pass for uninitialized-counter-poise test\n");
   dftd3_delete(error);
   dftd3_delete(mol);
   return 1;
}

int
test_invalid_structure (void)
{
   printf("Start test: invalid structure\n");
   dftd3_error error = NULL;
   dftd3_structure mol = NULL;

   int natoms = 2;
   int num[2] = {1, 1};
   double xyz[6] = {0.0};

   error = dftd3_new_error();

   mol = dftd3_new_structure(error, natoms, num, xyz, NULL, NULL);
   if (!dftd3_check_error(error)) goto unexpected;

   show_error(error);

   dftd3_delete(error);
   dftd3_delete(mol);
   return 0;

unexpected:
   printf("[Fatal] Unexpected pass for invalid-structure test\n");
   dftd3_delete(error);
   dftd3_delete(mol);
   return 1;
}

int
test_ghost_atoms(void) {
   double energy = 0.0;
   double gradient[21] = {0.0};
   dftd3_error error = dftd3_new_error();
   dftd3_structure mol = get_test_structure(error);
   dftd3_model disp = dftd3_new_d3_model(error, mol);
   dftd3_param param = dftd3_new_rational_damping(error, 1.0, 0.7875, 0.0, 0.4289, 4.4407, 14.0);
   int const ghost[7] = {0, 1, 2, 3, 4, 5, 6};

   if (dftd3_check_error(error)) return 1;
   dftd3_set_model_ghost_index(error, disp, ghost, 7);
   if (dftd3_check_error(error)) return 1;
   dftd3_get_dispersion(error, mol, disp, param, &energy, gradient, NULL);
   if (dftd3_check_error(error)) return 1;

   if (fabs(energy) > 1.0e-12) return 1;
   for (int i = 0; i < 21; ++i) {
      if (fabs(gradient[i]) > 1.0e-12) return 1;
   }

   dftd3_delete(param);
   dftd3_delete(disp);
   dftd3_delete(mol);
   dftd3_delete(error);
   return 0;
}

int
test_d3 (void) {
   double energy;
   double pair_disp2[49];
   double pair_disp3[49];
   double gradient[21];
   double sigma[9];
   double hessian[441];

   dftd3_error error;
   dftd3_structure mol;
   dftd3_model disp;
   dftd3_param param;

   error = dftd3_new_error();
   assert(!!error);

   mol = get_test_structure(error);
   if (dftd3_check_error(error)) {return 1;};
   assert(!!mol);

   disp = dftd3_new_d3_model(error, mol);
   if (dftd3_check_error(error)) {return 1;}
   assert(!!disp);

   // PBE-D3(BJ)
   param = dftd3_new_rational_damping(error, 1.0, 0.7875, 0.0, 0.4289, 4.4407, 14.0);
   if (dftd3_check_error(error)) {return 1;}
   assert(!!param);
   dftd3_get_dispersion(error, mol, disp, param, &energy, NULL, NULL);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_get_dispersion(error, mol, disp, param, &energy, gradient, sigma);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_get_pairwise_dispersion(error, mol, disp, param, pair_disp2, pair_disp3);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_delete(param);

   // RPBE-D3(0)
   param = dftd3_load_zero_damping(error, "rpbe", false);
   if (dftd3_check_error(error)) {return 1;}
   assert(!!param);
   dftd3_get_dispersion(error, mol, disp, param, &energy, NULL, NULL);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_get_dispersion(error, mol, disp, param, &energy, gradient, sigma);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_get_dispersion_hessian(error, mol, disp, param, &energy, hessian);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_delete(param);

   dftd3_set_model_realspace_cutoff_smooth(error, disp, 50.0, 30.0, 25.0, 5.0, 3.0);
   if (dftd3_check_error(error)) {return 1;}

   // DSD-BLYP-D3(BJ)-ATM
   param = dftd3_load_rational_damping(error, "dsdblyp", true);
   if (dftd3_check_error(error)) {return 1;}
   assert(!!param);
   dftd3_get_dispersion(error, mol, disp, param, &energy, NULL, NULL);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_get_dispersion(error, mol, disp, param, &energy, gradient, sigma);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_delete(param);

   // BLYP-D3(0)-ATM
   param = dftd3_new_zero_damping(error, 1.0, 1.682, 1.0, 1.094, 1.0, 14.0);
   if (dftd3_check_error(error)) {return 1;}
   assert(!!param);
   dftd3_get_dispersion(error, mol, disp, param, &energy, NULL, NULL);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_get_dispersion(error, mol, disp, param, &energy, gradient, sigma);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_delete(param);

   // B3LYP-D3(CSO)
   param = dftd3_new_cso_damping(error, 1.0, 0.0, 0.86, 2.5, 0.0, 6.25, 14.0);
   if (dftd3_check_error(error)) {return 1;}
   assert(!!param);
   dftd3_get_dispersion(error, mol, disp, param, &energy, NULL, NULL);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_get_dispersion(error, mol, disp, param, &energy, gradient, sigma);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_get_pairwise_dispersion(error, mol, disp, param, pair_disp2, pair_disp3);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_delete(param);

   // PBE-D3(Z)
   param = dftd3_new_z_damping(error, 1.0, 1.0, 0.0, 200770.0, 14.0);
   if (dftd3_check_error(error)) {return 1;}
   assert(!!param);
   dftd3_get_dispersion(error, mol, disp, param, &energy, NULL, NULL);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_get_dispersion(error, mol, disp, param, &energy, gradient, sigma);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_get_pairwise_dispersion(error, mol, disp, param, pair_disp2, pair_disp3);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_delete(param);

   // PBE-D3(Z)-ATM
   param = dftd3_new_z_damping(error, 1.0, 1.0, 1.0, 200770.0, 14.0);
   if (dftd3_check_error(error)) {return 1;}
   assert(!!param);
   dftd3_get_dispersion(error, mol, disp, param, &energy, NULL, NULL);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_get_dispersion(error, mol, disp, param, &energy, gradient, sigma);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_delete(param);

   // PBE-D3(CSO)-ATM
   param = dftd3_load_cso_damping(error, "pbe", true);
   if (dftd3_check_error(error)) {return 1;}
   assert(!!param);
   dftd3_get_dispersion(error, mol, disp, param, &energy, NULL, NULL);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_get_dispersion(error, mol, disp, param, &energy, gradient, sigma);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_delete(param);

   dftd3_delete(disp);
   dftd3_delete(mol);
   dftd3_delete(error);

   assert(!param);
   assert(!disp);
   assert(!mol);
   assert(!error);

   return 0;
}

int
test_gcp (void) {
   double energy;
   double pair_disp2[49];
   double pair_disp3[49];
   double gradient[21];
   double sigma[9];

   dftd3_error error;
   dftd3_structure mol;
   dftd3_gcp gcp;

   error = dftd3_new_error();
   assert(!!error);

   mol = get_test_structure(error);
   if (dftd3_check_error(error)) {return 1;};
   assert(!!mol);

   gcp = dftd3_load_gcp_param(error, mol, "pbeh3c", NULL);
   if (dftd3_check_error(error)) {return 1;}
   assert(!!gcp);

   dftd3_get_counterpoise(error, mol, gcp, &energy, NULL, NULL);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_get_counterpoise(error, mol, gcp, &energy, gradient, sigma);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_delete(gcp);

   gcp = dftd3_load_gcp_param(error, mol, "b973c", NULL);
   if (dftd3_check_error(error)) {return 1;}
   assert(!!gcp);

   dftd3_set_gcp_realspace_cutoff(error, gcp, 50.0, 30.0);
   if (dftd3_check_error(error)) {return 1;}

   dftd3_get_counterpoise(error, mol, gcp, &energy, NULL, NULL);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_get_counterpoise(error, mol, gcp, &energy, gradient, sigma);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_delete(gcp);

   dftd3_delete(mol);
   dftd3_delete(error);

   assert(!gcp);
   assert(!mol);
   assert(!error);

   return 0;
}

int
test_hessian (void)
{
   printf("Start test: hessian\n");
   int const natoms = 7;
   int const ndim = 21;
   int const attyp[7] = {6,6,6,1,1,1,1};
   double coord[21] =
      {0.00000000000000, 0.00000000000000,-1.79755622305860,
       0.00000000000000, 0.00000000000000, 0.95338756106749,
       0.00000000000000, 0.00000000000000, 3.22281255790261,
      -0.96412815539807,-1.66991895015711,-2.53624948351102,
      -0.96412815539807, 1.66991895015711,-2.53624948351102,
       1.92825631079613, 0.00000000000000,-2.53624948351102,
       0.00000000000000, 0.00000000000000, 5.23010455462158};
   double ref[21];
   double hessian[441];
   double numhess[441];
   double gradr[21], gradl[21], sigma[9];
   double energy;
   double const step = 1.0e-6;
   int stat = 0;

   dftd3_error error = NULL;
   dftd3_structure mol = NULL;
   dftd3_model disp = NULL;
   dftd3_param param = NULL;

   error = dftd3_new_error();
   mol = dftd3_new_structure(error, natoms, attyp, coord, NULL, NULL);
   if (dftd3_check_error(error)) {stat = 1; goto cleanup;}

   disp = dftd3_new_d3_model(error, mol);
   if (dftd3_check_error(error)) {stat = 1; goto cleanup;}

   // DSD-BLYP-D3(BJ)-ATM
   param = dftd3_load_rational_damping(error, "dsdblyp", true);
   if (dftd3_check_error(error)) {stat = 1; goto cleanup;}

   dftd3_get_dispersion_hessian(error, mol, disp, param, &energy, hessian);
   if (dftd3_check_error(error)) {stat = 1; goto cleanup;}

   for (int i = 0; i < ndim; i++) ref[i] = coord[i];

   for (int i = 0; i < ndim; i++) {
      coord[i] = ref[i] + step;
      dftd3_update_structure(error, mol, coord, NULL);
      dftd3_get_dispersion(error, mol, disp, param, &energy, gradr, sigma);
      if (dftd3_check_error(error)) {stat = 1; goto cleanup;}

      coord[i] = ref[i] - step;
      dftd3_update_structure(error, mol, coord, NULL);
      dftd3_get_dispersion(error, mol, disp, param, &energy, gradl, sigma);
      if (dftd3_check_error(error)) {stat = 1; goto cleanup;}

      coord[i] = ref[i];
      for (int j = 0; j < ndim; j++) {
         numhess[i*ndim + j] = 0.5 * (gradr[j] - gradl[j]) / step;
      }
   }

   for (int i = 0; i < ndim; i++) {
      for (int j = 0; j < ndim; j++) {
         if (fabs(hessian[i*ndim + j] - hessian[j*ndim + i]) > 1.0e-12) {
            printf("[Fatal] Hessian is not symmetric\n");
            stat = 1;
            goto cleanup;
         }
         if (fabs(hessian[i*ndim + j] - numhess[i*ndim + j]) > 1.0e-7) {
            printf("[Fatal] Hessian does not match finite differences\n");
            stat = 1;
            goto cleanup;
         }
      }
   }

cleanup:
   dftd3_delete(param);
   dftd3_delete(disp);
   dftd3_delete(mol);
   dftd3_delete(error);
   return stat;
}

int
test_gcp_hessian (void)
{
   printf("Start test: counter-poise hessian\n");
   int const natoms = 7;
   int const ndim = 21;
   int const attyp[7] = {6,6,6,1,1,1,1};
   double coord[21] =
      {0.00000000000000, 0.00000000000000,-1.79755622305860,
       0.00000000000000, 0.00000000000000, 0.95338756106749,
       0.00000000000000, 0.00000000000000, 3.22281255790261,
      -0.96412815539807,-1.66991895015711,-2.53624948351102,
      -0.96412815539807, 1.66991895015711,-2.53624948351102,
       1.92825631079613, 0.00000000000000,-2.53624948351102,
       0.00000000000000, 0.00000000000000, 5.23010455462158};
   double ref[21];
   double hessian[441];
   double numhess[441];
   double gradr[21], gradl[21], sigma[9];
   double energy;
   double const step = 1.0e-6;
   int stat = 0;

   dftd3_error error = NULL;
   dftd3_structure mol = NULL;
   dftd3_gcp gcp = NULL;

   error = dftd3_new_error();
   mol = dftd3_new_structure(error, natoms, attyp, coord, NULL, NULL);
   if (dftd3_check_error(error)) {stat = 1; goto cleanup;}

   // HF-3c includes the counter-poise, short-range bond and basis corrections
   gcp = dftd3_load_gcp_param(error, mol, "hf3c", NULL);
   if (dftd3_check_error(error)) {stat = 1; goto cleanup;}

   dftd3_get_counterpoise_hessian(error, mol, gcp, &energy, hessian);
   if (dftd3_check_error(error)) {stat = 1; goto cleanup;}

   for (int i = 0; i < ndim; i++) ref[i] = coord[i];

   for (int i = 0; i < ndim; i++) {
      coord[i] = ref[i] + step;
      dftd3_update_structure(error, mol, coord, NULL);
      dftd3_get_counterpoise(error, mol, gcp, &energy, gradr, sigma);
      if (dftd3_check_error(error)) {stat = 1; goto cleanup;}

      coord[i] = ref[i] - step;
      dftd3_update_structure(error, mol, coord, NULL);
      dftd3_get_counterpoise(error, mol, gcp, &energy, gradl, sigma);
      if (dftd3_check_error(error)) {stat = 1; goto cleanup;}

      coord[i] = ref[i];
      for (int j = 0; j < ndim; j++) {
         numhess[i*ndim + j] = 0.5 * (gradr[j] - gradl[j]) / step;
      }
   }

   for (int i = 0; i < ndim; i++) {
      for (int j = 0; j < ndim; j++) {
         if (fabs(hessian[i*ndim + j] - hessian[j*ndim + i]) > 1.0e-12) {
            printf("[Fatal] Counter-poise hessian is not symmetric\n");
            stat = 1;
            goto cleanup;
         }
         if (fabs(hessian[i*ndim + j] - numhess[i*ndim + j]) > 1.0e-7) {
            printf("[Fatal] Counter-poise hessian does not match finite differences\n");
            stat = 1;
            goto cleanup;
         }
      }
   }

cleanup:
   dftd3_delete(gcp);
   dftd3_delete(mol);
   dftd3_delete(error);
   return stat;
}

int
test_ghost_atoms_atm(void) {
   printf("Start test: ghost atoms with ATM\n");
   double energy = 0.0;
   double gradient[21] = {0.0};
   double hessian[441] = {0.0};
   // ghosting a subset while the three-body term is active zeroes single C6
   // coefficients, which must not produce NaN in the derivatives
   int const ghost[2] = {0, 3};
   int stat = 0;

   dftd3_error error = dftd3_new_error();
   dftd3_structure mol = get_test_structure(error);
   dftd3_model disp = dftd3_new_d3_model(error, mol);
   // PBE-D3(BJ)-ATM
   dftd3_param param =
      dftd3_new_rational_damping(error, 1.0, 0.7875, 1.0, 0.4289, 4.4407, 14.0);
   if (dftd3_check_error(error)) {stat = 1; goto cleanup;}

   dftd3_set_model_ghost_index(error, disp, ghost, 2);
   if (dftd3_check_error(error)) {stat = 1; goto cleanup;}

   dftd3_get_dispersion(error, mol, disp, param, &energy, gradient, NULL);
   if (dftd3_check_error(error)) {stat = 1; goto cleanup;}

   if (isnan(energy)) {
      printf("[Fatal] Energy is NaN for ghost atoms with ATM\n");
      stat = 1;
      goto cleanup;
   }
   for (int i = 0; i < 21; ++i) {
      if (isnan(gradient[i])) {
         printf("[Fatal] Gradient is NaN for ghost atoms with ATM\n");
         stat = 1;
         goto cleanup;
      }
   }

   dftd3_get_dispersion_hessian(error, mol, disp, param, &energy, hessian);
   if (dftd3_check_error(error)) {stat = 1; goto cleanup;}
   for (int i = 0; i < 441; ++i) {
      if (isnan(hessian[i])) {
         printf("[Fatal] Hessian is NaN for ghost atoms with ATM\n");
         stat = 1;
         goto cleanup;
      }
   }

cleanup:
   dftd3_delete(param);
   dftd3_delete(disp);
   dftd3_delete(mol);
   dftd3_delete(error);
   return stat;
}

int
test_ewald (void)
{
   printf("Start test: Ewald summation\n");
   double converged, ewald, truncated;
   double gradient[6], sigma[9];

   dftd3_error error = NULL;
   dftd3_structure mol = NULL;
   dftd3_model disp = NULL;
   dftd3_param param = NULL;

   error = dftd3_new_error();
   mol = get_periodic_test_structure(error);
   if (dftd3_check_error(error)) return 1;

   param = dftd3_load_rational_damping(error, "pbe", false);
   if (dftd3_check_error(error)) return 1;

   /* converged real space reference, the coordination number cutoff has to
      match the one used for the reciprocal space evaluation */
   disp = dftd3_new_d3_model(error, mol);
   if (dftd3_check_error(error)) return 1;
   dftd3_set_model_realspace_cutoff(error, disp, 400.0, 40.0, 40.0);
   if (dftd3_check_error(error)) return 1;
   dftd3_get_dispersion(error, mol, disp, param, &converged, NULL, NULL);
   if (dftd3_check_error(error)) return 1;
   dftd3_delete(disp);

   /* truncated real space evaluation with the default cutoff */
   disp = dftd3_new_d3_model(error, mol);
   if (dftd3_check_error(error)) return 1;
   dftd3_get_dispersion(error, mol, disp, param, &truncated, NULL, NULL);
   if (dftd3_check_error(error)) return 1;
   dftd3_delete(disp);

   /* reciprocal space evaluation, the two-body cutoff is ignored */
   disp = dftd3_new_d3_model(error, mol);
   if (dftd3_check_error(error)) return 1;
   dftd3_set_model_ewald(error, disp, 0, 0.0, 10.0);
   if (dftd3_check_error(error)) return 1;
   dftd3_get_dispersion(error, mol, disp, param, &ewald, gradient, sigma);
   if (dftd3_check_error(error)) return 1;

   printf("[Info] converged %.12f  ewald %.12f  truncated %.12f\n",
          converged, ewald, truncated);

   if (fabs(ewald - converged) > 1.0e-7) {
      printf("[Fatal] Ewald summation does not match converged real space\n");
      return 1;
   }

   if (fabs(truncated - converged) < 1.0e-6) {
      printf("[Fatal] Real space truncation error is unexpectedly small\n");
      return 1;
   }

   dftd3_delete(param);
   dftd3_delete(disp);
   dftd3_delete(mol);
   dftd3_delete(error);
   return 0;
}

int
test_ewald_unsupported (void)
{
   printf("Start test: Ewald summation with unsupported damping\n");
   double energy;

   dftd3_error error = NULL;
   dftd3_structure mol = NULL;
   dftd3_model disp = NULL;
   dftd3_param param = NULL;

   error = dftd3_new_error();

   /* setting up the expansion without a model has to be reported */
   dftd3_set_model_ewald(error, disp, 0, 0.0, 0.0);
   if (!dftd3_check_error(error)) goto unexpected;
   show_error(error);

   mol = get_periodic_test_structure(error);
   if (dftd3_check_error(error)) return 1;

   disp = dftd3_new_d3_model(error, mol);
   if (dftd3_check_error(error)) return 1;
   dftd3_set_model_ewald(error, disp, 0, 0.0, 0.0);
   if (dftd3_check_error(error)) return 1;

   /* the modified zero damping has no reciprocal space representation */
   param = dftd3_load_mzero_damping(error, "pbe", false);
   if (dftd3_check_error(error)) return 1;

   dftd3_get_dispersion(error, mol, disp, param, &energy, NULL, NULL);
   if (!dftd3_check_error(error)) goto unexpected;
   show_error(error);

   dftd3_delete(param);
   dftd3_delete(disp);
   dftd3_delete(mol);
   dftd3_delete(error);
   return 0;

unexpected:
   printf("[Fatal] Unexpected pass for unsupported-damping test\n");
   dftd3_delete(param);
   dftd3_delete(disp);
   dftd3_delete(mol);
   dftd3_delete(error);
   return 1;
}

int
test_ewald_realspace_only (void)
{
   printf("Start test: Ewald summation with real space only properties\n");
   double energy;
   double pair_disp2[4], pair_disp3[4];
   double hessian[36];

   dftd3_error error = NULL;
   dftd3_structure mol = NULL;
   dftd3_model disp = NULL;
   dftd3_param param = NULL;

   error = dftd3_new_error();
   mol = get_periodic_test_structure(error);
   if (dftd3_check_error(error)) return 1;

   param = dftd3_load_rational_damping(error, "pbe", false);
   if (dftd3_check_error(error)) return 1;

   disp = dftd3_new_d3_model(error, mol);
   if (dftd3_check_error(error)) return 1;
   dftd3_set_model_ewald(error, disp, 0, 0.0, 10.0);
   if (dftd3_check_error(error)) return 1;

   /* the pairwise decomposition would not add up to the reciprocal space energy */
   dftd3_get_pairwise_dispersion(error, mol, disp, param, pair_disp2, pair_disp3);
   if (!dftd3_check_error(error)) goto unexpected;
   show_error(error);

   /* the second derivatives are only implemented in real space */
   dftd3_get_dispersion_hessian(error, mol, disp, param, &energy, hessian);
   if (!dftd3_check_error(error)) goto unexpected;
   show_error(error);

   /* the energy itself is still available */
   dftd3_get_dispersion(error, mol, disp, param, &energy, NULL, NULL);
   if (dftd3_check_error(error)) return 1;

   dftd3_delete(param);
   dftd3_delete(disp);
   dftd3_delete(mol);
   dftd3_delete(error);
   return 0;

unexpected:
   printf("[Fatal] Unexpected pass for real-space-only property\n");
   dftd3_delete(param);
   dftd3_delete(disp);
   dftd3_delete(mol);
   dftd3_delete(error);
   return 1;
}

int
test_work_partition (void)
{
   printf("Start test: work partition\n");
   static const int nparts = 3;
   static const int nat = 7;
   double energy, part_energy, sum_energy;
   double gradient[21], part_gradient[21], sum_gradient[21];
   double sigma[9], part_sigma[9], sum_sigma[9];
   double hessian[441], part_hessian[441], sum_hessian[441];
   int i, part;

   dftd3_error error = NULL;
   dftd3_structure mol = NULL;
   dftd3_model disp = NULL;
   dftd3_gcp gcp = NULL;
   dftd3_param param = NULL;

   error = dftd3_new_error();
   mol = get_test_structure(error);
   if (dftd3_check_error(error)) return 1;

   param = dftd3_load_rational_damping(error, "pbe", true);
   if (dftd3_check_error(error)) return 1;

   disp = dftd3_new_d3_model(error, mol);
   if (dftd3_check_error(error)) return 1;
   dftd3_get_dispersion(error, mol, disp, param, &energy, gradient, sigma);
   if (dftd3_check_error(error)) return 1;
   dftd3_get_dispersion_hessian(error, mol, disp, param, &part_energy, hessian);
   if (dftd3_check_error(error)) return 1;
   dftd3_delete(disp);

   sum_energy = 0.0;
   for (i = 0; i < 3 * nat; i++) sum_gradient[i] = 0.0;
   for (i = 0; i < 9; i++) sum_sigma[i] = 0.0;
   for (i = 0; i < 9 * nat * nat; i++) sum_hessian[i] = 0.0;

   for (part = 0; part < nparts; part++) {
      disp = dftd3_new_d3_model(error, mol);
      if (dftd3_check_error(error)) return 1;
      dftd3_set_model_work_partition(error, disp, part, nparts);
      if (dftd3_check_error(error)) return 1;

      dftd3_get_dispersion(error, mol, disp, param, &part_energy, part_gradient, part_sigma);
      if (dftd3_check_error(error)) return 1;
      dftd3_get_dispersion_hessian(error, mol, disp, param, &part_energy, part_hessian);
      if (dftd3_check_error(error)) return 1;
      dftd3_delete(disp);

      sum_energy += part_energy;
      for (i = 0; i < 3 * nat; i++) sum_gradient[i] += part_gradient[i];
      for (i = 0; i < 9; i++) sum_sigma[i] += part_sigma[i];
      for (i = 0; i < 9 * nat * nat; i++) sum_hessian[i] += part_hessian[i];
   }

   if (fabs(sum_energy - energy) > 1.0e-12) {
      printf("[Fatal] Partitioned dispersion energy does not match\n");
      return 1;
   }
   for (i = 0; i < 3 * nat; i++) {
      if (fabs(sum_gradient[i] - gradient[i]) > 1.0e-12) {
         printf("[Fatal] Partitioned dispersion gradient does not match\n");
         return 1;
      }
   }
   for (i = 0; i < 9; i++) {
      if (fabs(sum_sigma[i] - sigma[i]) > 1.0e-12) {
         printf("[Fatal] Partitioned dispersion virial does not match\n");
         return 1;
      }
   }
   for (i = 0; i < 9 * nat * nat; i++) {
      if (fabs(sum_hessian[i] - hessian[i]) > 1.0e-12) {
         printf("[Fatal] Partitioned dispersion hessian does not match\n");
         return 1;
      }
   }

   /* the counterpoise correction carries its own partition */
   gcp = dftd3_load_gcp_param(error, mol, "hf3c", NULL);
   if (dftd3_check_error(error)) return 1;
   dftd3_get_counterpoise(error, mol, gcp, &energy, gradient, sigma);
   if (dftd3_check_error(error)) return 1;
   dftd3_delete(gcp);

   sum_energy = 0.0;
   for (i = 0; i < 3 * nat; i++) sum_gradient[i] = 0.0;
   for (part = 0; part < nparts; part++) {
      gcp = dftd3_load_gcp_param(error, mol, "hf3c", NULL);
      if (dftd3_check_error(error)) return 1;
      dftd3_set_gcp_work_partition(error, gcp, part, nparts);
      if (dftd3_check_error(error)) return 1;
      dftd3_get_counterpoise(error, mol, gcp, &part_energy, part_gradient, part_sigma);
      if (dftd3_check_error(error)) return 1;
      dftd3_delete(gcp);

      sum_energy += part_energy;
      for (i = 0; i < 3 * nat; i++) sum_gradient[i] += part_gradient[i];
   }

   if (fabs(sum_energy - energy) > 1.0e-12) {
      printf("[Fatal] Partitioned counterpoise energy does not match\n");
      return 1;
   }
   for (i = 0; i < 3 * nat; i++) {
      if (fabs(sum_gradient[i] - gradient[i]) > 1.0e-12) {
         printf("[Fatal] Partitioned counterpoise gradient does not match\n");
         return 1;
      }
   }

   dftd3_delete(param);
   dftd3_delete(mol);
   dftd3_delete(error);
   return 0;
}

int
test_invalid_work_partition (void)
{
   printf("Start test: invalid work partition\n");

   dftd3_error error = NULL;
   dftd3_structure mol = NULL;
   dftd3_model disp = NULL;
   dftd3_gcp gcp = NULL;
   int i;

   /* part index, number of parts */
   int const invalid[4][2] = {{-1, 3}, {3, 3}, {0, 0}, {1, 1}};

   error = dftd3_new_error();
   mol = get_test_structure(error);
   if (dftd3_check_error(error)) return 1;
   disp = dftd3_new_d3_model(error, mol);
   if (dftd3_check_error(error)) return 1;
   gcp = dftd3_load_gcp_param(error, mol, "hf3c", NULL);
   if (dftd3_check_error(error)) return 1;

   /* the error handle is not cleared by reading it, use a fresh one per case */
   for (i = 0; i < 4; i++) {
      dftd3_delete(error);
      error = dftd3_new_error();
      dftd3_set_model_work_partition(error, disp, invalid[i][0], invalid[i][1]);
      if (!dftd3_check_error(error)) goto unexpected;
      show_error(error);

      dftd3_delete(error);
      error = dftd3_new_error();
      dftd3_set_gcp_work_partition(error, gcp, invalid[i][0], invalid[i][1]);
      if (!dftd3_check_error(error)) goto unexpected;
      show_error(error);
   }

   /* assigning a partition without a handle has to be reported */
   dftd3_delete(error);
   error = dftd3_new_error();
   dftd3_set_model_work_partition(error, NULL, 0, 1);
   if (!dftd3_check_error(error)) goto unexpected;
   show_error(error);

   dftd3_delete(error);
   error = dftd3_new_error();
   dftd3_set_gcp_work_partition(error, NULL, 0, 1);
   if (!dftd3_check_error(error)) goto unexpected;
   show_error(error);

   dftd3_delete(gcp);
   dftd3_delete(disp);
   dftd3_delete(mol);
   dftd3_delete(error);
   return 0;

unexpected:
   printf("[Fatal] Unexpected pass for invalid work partition\n");
   dftd3_delete(gcp);
   dftd3_delete(disp);
   dftd3_delete(mol);
   dftd3_delete(error);
   return 1;
}


int
test_mpi_comm (void)
{
   printf("Start test: mpi communicator\n");

   dftd3_error error = NULL;
   dftd3_structure mol = NULL;
   dftd3_model disp = NULL;
   dftd3_gcp gcp = NULL;

   if (dftd3_has_feature("this-is-not-a-feature")) {
      printf("[Fatal] Unknown feature reported as available\n");
      return 1;
   }

   error = dftd3_new_error();
   mol = get_test_structure(error);
   if (dftd3_check_error(error)) return 1;
   disp = dftd3_new_d3_model(error, mol);
   if (dftd3_check_error(error)) return 1;
   gcp = dftd3_load_gcp_param(error, mol, "hf3c", NULL);
   if (dftd3_check_error(error)) return 1;

   /* MPI is never initialized here, without the feature the missing build is
    * reported instead, either way no communicator may be accepted */
   dftd3_delete(error);
   error = dftd3_new_error();
   dftd3_set_model_mpi_comm(error, disp, 0);
   if (!dftd3_check_error(error)) goto unexpected;
   show_error(error);

   dftd3_delete(error);
   error = dftd3_new_error();
   dftd3_set_gcp_mpi_comm(error, gcp, 0);
   if (!dftd3_check_error(error)) goto unexpected;
   show_error(error);

   dftd3_delete(error);
   error = dftd3_new_error();
   dftd3_set_model_mpi_comm(error, NULL, 0);
   if (!dftd3_check_error(error)) goto unexpected;
   show_error(error);

   dftd3_delete(error);
   error = dftd3_new_error();
   dftd3_set_gcp_mpi_comm(error, NULL, 0);
   if (!dftd3_check_error(error)) goto unexpected;
   show_error(error);

   dftd3_delete(gcp);
   dftd3_delete(disp);
   dftd3_delete(mol);
   dftd3_delete(error);
   return 0;

unexpected:
   printf("[Fatal] Unexpected pass for uninitialized MPI\n");
   dftd3_delete(gcp);
   dftd3_delete(disp);
   dftd3_delete(mol);
   dftd3_delete(error);
   return 1;
}

int
main (void)
{
   int stat = 0;
   stat += test_version();
   stat += test_uninitialized_error();
   stat += test_uninitialized_structure();
   stat += test_uninitialized_model();
   stat += test_uninitialized_param();
   stat += test_uninitialized_gcp();
   stat += test_invalid_structure();
   stat += test_ghost_atoms();
   stat += test_ghost_atoms_atm();
   stat += test_d3();
   stat += test_hessian();
   stat += test_gcp();
   stat += test_gcp_hessian();
   stat += test_ewald();
   stat += test_ewald_unsupported();
   stat += test_ewald_realspace_only();
   stat += test_work_partition();
   stat += test_invalid_work_partition();
   stat += test_mpi_comm();
   return stat;
}
