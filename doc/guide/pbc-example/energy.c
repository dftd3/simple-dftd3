#include <stdio.h>
#include <stdbool.h>

#include "dftd3.h"

int main(void)
{
  dftd3_error error = dftd3_new_error();
  dftd3_structure mol = NULL;
  dftd3_model d3 = NULL;
  dftd3_param param = NULL;

  int nat = 16;
  int num[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 7, 7, 7};
  double xyz[48] = {  // coordinates in Bohr
     3.5141339317522400, 2.5620901186651399, 0.9887044918760900,
     2.5620901186651399, 0.9887044918760900, 3.5141339317522400,
     0.9887044918760900, 3.5141339317522400, 2.5620901186651399,
     8.2027323927494997, 2.1265083423321101, 8.3884924301184203,
     7.2504996070913199, 3.6998939691211699, 5.8628740176711904,
     5.6773029528733403, 1.1742755566739300, 6.8151068033293702,
     6.8151068033293702, 5.6773029528733403, 1.1742755566739300,
     8.3884924301184203, 8.2027323927494997, 2.1265083423321101,
     3.6998939691211699, 5.8628740176711904, 7.2504996070913199,
     1.1742755566739300, 6.8151068033293702, 5.6773029528733403,
     2.1265083423321101, 8.3884924301184203, 8.2027323927494997,
     5.8628740176711904, 7.2504996070913199, 3.6998939691211699,
     1.9543543300807500, 1.9543543300807500, 1.9543543300807500,
     6.6427638185069302, 2.7342441309165002, 7.4228425919137502,
     7.4228425919137502, 6.6427638185069302, 2.7342441309165002,
     2.7342441309165002, 7.4228425919137502, 6.6427638185069302};
   double lattice[9] = {  // lattice vectors in Bohr
     9.3771283249512098, 0.0000000000000000, 0.0000000000000000,
     0.0000000000000000, 9.3771283249512098, 0.0000000000000000,
     0.0000000000000000, 0.0000000000000000, 9.3771283249512098};
  mol = dftd3_new_structure(error, nat, num, xyz, lattice, NULL);
  if (dftd3_check_error(error)) goto handle_error;

  char method[7] = "r2SCAN";
  param = dftd3_load_rational_damping(error, method, true);
  if (dftd3_check_error(error)) goto handle_error;

  d3 = dftd3_new_d3_model(error, mol);
  if (dftd3_check_error(error)) goto handle_error;

  // the defaults truncate the two-body sum at 60 Bohr and the coordination
  // number at 40 Bohr, all cutoffs are given in Bohr
  dftd3_set_model_realspace_cutoff(error, d3, 60.0, 40.0, 40.0);
  if (dftd3_check_error(error)) goto handle_error;

  // to get the gradient and virial, provide arrays for the output
  double energy, gradient[48], virial[9];
  dftd3_get_dispersion(error, mol, d3, param, &energy, gradient, virial);
  if (dftd3_check_error(error)) goto handle_error;

  printf("Dispersion energy for %s-D3(BJ) is %13.10lf Hartree\n", method, energy);

  dftd3_delete(error);
  dftd3_delete(mol);
  dftd3_delete(d3);
  dftd3_delete(param);
  return 0;

handle_error: {
  char msg[512];
  dftd3_get_error(error, msg, NULL);
  printf("Error: %s\n", msg);

  dftd3_delete(error);
  dftd3_delete(mol);
  dftd3_delete(d3);
  dftd3_delete(param);
  return 1;
}
}