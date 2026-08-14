program mpi_partition
   use, intrinsic :: iso_fortran_env, only : r8 => real64
   use dftd3, only: d3_model, d3_param, rational_damping_param, get_rational_damping, &
      & new_rational_damping, new_d3_model, get_dispersion, realspace_cutoff, &
      & work_partition, new_work_partition
   use mctc_env, only: error_type
   use mctc_io, only: structure_type, new
   use mpi_f08, only : MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_SUM, &
      & MPI_Abort, MPI_Finalize, MPI_Init, MPI_Comm_rank, MPI_Comm_size, &
      & MPI_Allreduce
   implicit none

   character(len=:), allocatable :: method
   type(structure_type) :: mol
   type(error_type), allocatable :: error
   integer, allocatable :: num(:)
   real(r8), allocatable :: xyz(:, :)
   real(r8), allocatable :: gradient(:, :)
   real(r8) :: energy

   type(d3_model) :: disp
   type(d3_param) :: inp
   type(rational_damping_param) :: param
   type(work_partition) :: partition
   integer :: rank, nranks

   call MPI_Init()
   call MPI_Comm_rank(MPI_COMM_WORLD, rank)
   call MPI_Comm_size(MPI_COMM_WORLD, nranks)

   method = "PBE0"
   num = [6, 1, 1, 1, 1]
   xyz = reshape([ &  ! coordinates in Bohr
     &  0.0000000_r8, -0.0000000_r8,  0.0000000_r8, &
     & -1.1922080_r8,  1.1922080_r8,  1.1922080_r8, &
     &  1.1922080_r8, -1.1922080_r8,  1.1922080_r8, &
     & -1.1922080_r8, -1.1922080_r8, -1.1922080_r8, &
     &  1.1922080_r8,  1.1922080_r8, -1.1922080_r8],&
     & [3, size(num)])
   call new(mol, num, xyz, charge=0.0_r8, uhf=0)

   call get_rational_damping(inp, method, error, s9=1.0_r8)
   if (allocated(error)) call fatal(error%message)
   call new_rational_damping(param, inp)
   call new_d3_model(disp, mol)

   ! every rank holds the complete structure and evaluates its own share of the
   ! pairwise and three-body interactions, parts are zero based
   call new_work_partition(error, partition, rank, nranks)
   if (allocated(error)) call fatal(error%message)

   allocate(gradient(3, mol%nat))
   call get_dispersion(error, mol, disp, param, realspace_cutoff(), energy, gradient, &
      & partition=partition)
   if (allocated(error)) call fatal(error%message)

   ! s-dftd3 performs no communication, the caller reduces the partial results
   call MPI_Allreduce(MPI_IN_PLACE, energy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD)
   call MPI_Allreduce(MPI_IN_PLACE, gradient, size(gradient), MPI_DOUBLE_PRECISION, &
      & MPI_SUM, MPI_COMM_WORLD)

   if (rank == 0) then
      print "(3a, f13.10, a)", "Dispersion energy for ", method, "-D3(BJ) is ", energy, " Hartree"
   end if

   call MPI_Finalize()

contains

   subroutine fatal(message)
      character(len=*), intent(in) :: message
      print "(2a)", "Error: ", message
      call MPI_Abort(MPI_COMM_WORLD, 1)
   end subroutine fatal

end program mpi_partition
