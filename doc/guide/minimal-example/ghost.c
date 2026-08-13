#include <stdio.h>

#include "dftd3.h"

int main(void)
{
  dftd3_error error = dftd3_new_error();
  dftd3_structure mol = NULL;
  dftd3_model d3 = NULL;
  dftd3_param param = NULL;

  int nat = 6;
  int num[6] = {8, 1, 1, 8, 1, 1};
  double xyz[18] = {
    -2.5142928288858868,  1.8879764891638844,  0.0000000000000000,
    -2.0098096791094657,  3.6303285583904779,  0.0000000000000000,
    -0.9545622395768059,  0.92963588958819932, 0.0000000000000000,
     1.8995798491583649, -1.4196911493711153,  0.0000000000000000,
     1.7895424492078418, -2.5141248888855077, -1.4467274293585810,
     1.7895424492078418, -2.5141248888855077,  1.4467274293585810};
  int ghost[3] = {3, 4, 5};
  double energy;

  mol = dftd3_new_structure(error, nat, num, xyz, NULL, NULL);
  if (dftd3_check_error(error)) goto handle_error;

  param = dftd3_load_rational_damping(error, "PBE0", false);
  if (dftd3_check_error(error)) goto handle_error;

  d3 = dftd3_new_d3_model(error, mol);
  if (dftd3_check_error(error)) goto handle_error;

  dftd3_set_model_ghost_index(error, d3, ghost, 3);
  if (dftd3_check_error(error)) goto handle_error;

  dftd3_get_dispersion(error, mol, d3, param, &energy, NULL, NULL);
  if (dftd3_check_error(error)) goto handle_error;

  printf("Dispersion energy with ghost atoms is %13.10lf Hartree\n", energy);

  dftd3_delete(error);
  dftd3_delete(mol);
  dftd3_delete(d3);
  dftd3_delete(param);
  return 0;

handle_error:
  char msg[512];
  dftd3_get_error(error, msg, NULL);
  printf("Error: %s\n", msg);

  dftd3_delete(error);
  dftd3_delete(mol);
  dftd3_delete(d3);
  dftd3_delete(param);
  return 1;
}
