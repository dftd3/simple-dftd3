#include <stdio.h>
#include <stdbool.h>

#include <mpi.h>

#include "dftd3.h"

int main(int argc, char **argv)
{
  dftd3_error error = dftd3_new_error();
  dftd3_structure mol = NULL;
  dftd3_model d3 = NULL;
  dftd3_param param = NULL;
  int rank, nranks;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);

  int nat = 5;
  int num[5] = {6, 1, 1, 1, 1};
  double xyz[15] = {  // coordinates in Bohr
     0.00000000, -0.00000000,  0.00000000,
    -1.19220800,  1.19220800,  1.19220800,
     1.19220800, -1.19220800,  1.19220800,
    -1.19220800, -1.19220800, -1.19220800,
     1.19220800,  1.19220800, -1.19220800};
  mol = dftd3_new_structure(error, nat, num, xyz, NULL, NULL);
  if (dftd3_check_error(error)) goto handle_error;

  char method[5] = "PBE0";
  param = dftd3_load_rational_damping(error, method, true);
  if (dftd3_check_error(error)) goto handle_error;

  d3 = dftd3_new_d3_model(error, mol);
  if (dftd3_check_error(error)) goto handle_error;

  // every rank holds the complete structure and evaluates its own share of the
  // pairwise and three-body interactions, parts are zero based
  dftd3_set_model_work_partition(error, d3, rank, nranks);
  if (dftd3_check_error(error)) goto handle_error;

  double energy;
  double gradient[15];
  dftd3_get_dispersion(error, mol, d3, param, &energy, gradient, NULL);
  if (dftd3_check_error(error)) goto handle_error;

  // s-dftd3 performs no communication, the caller reduces the partial results
  MPI_Allreduce(MPI_IN_PLACE, &energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, gradient, 3 * nat, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (rank == 0) {
    printf("Dispersion energy for %s-D3(BJ) is %13.10lf Hartree\n", method, energy);
  }

  dftd3_delete(error);
  dftd3_delete(mol);
  dftd3_delete(d3);
  dftd3_delete(param);
  MPI_Finalize();
  return 0;

handle_error: {
  char msg[512];
  dftd3_get_error(error, msg, NULL);
  printf("Error: %s\n", msg);

  dftd3_delete(error);
  dftd3_delete(mol);
  dftd3_delete(d3);
  dftd3_delete(param);
  MPI_Abort(MPI_COMM_WORLD, 1);
  return 1;
}
}
