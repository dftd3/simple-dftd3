!> molecular structure information
! 
!  contains information about the molecular geometry, the atom types
!  the nuclear and total charge, atomic masses and all interatomic distances
!  In periodic calculations the lattice parameters, a Wigner--Seitz cell
!  as well as fractional coordinates are attached to the data type
module d3def_molecule
   use iso_fortran_env, only: wp => real64
   use iso_c_binding
   implicit none
   public :: d3_molecule
   public :: len
   private

   integer, parameter :: small = 30
   integer, parameter :: large = 90

   !> molecular structure information
   type :: d3_molecule
      integer  :: n = 0            !< number of atoms
      real(wp) :: chrg = 0.0_wp    !< total charge
      integer  :: uhf = 0          !< number of unpaired electrons
      logical  :: pbc(3) = .false. !< periodic dimensions
      integer  :: npbc = 0         !< periodicity of system
      character(len=2),allocatable :: sym(:) !< element symbols
      integer, allocatable :: at(:)          !< ordinal numbers
      real(wp),allocatable :: xyz(:,:)       !< cartesian coordinates in bohr
      real(wp),allocatable :: abc(:,:)       !< fractional coordinates
      real(wp),allocatable :: dist(:,:)      !< interatomic distances
      real(wp),allocatable :: atmass(:)      !< atomic masses in amu
      real(wp),allocatable :: z(:)           !< nuclear charges
      real(wp),allocatable :: cn(:)          !< coordination number
      real(wp) :: cellpar(6) = 0.0_wp        !< cell parameters
      real(wp) :: lattice(3,3) = 0.0_wp      !< direct lattice parameters
      real(wp) :: rec_lat(3,3) = 0.0_wp      !< reciprocal lattice parameters
      real(wp) :: volume = 0.0_wp            !< volume of unit cell
      character(len=:), allocatable :: name
      integer :: ftype
   contains
      generic :: new => new_from_natoms
      procedure, private :: new_from_natoms => molecule_new_from_natoms
      procedure :: read => read_molecule_generic
      procedure :: print_info => molecule_print_info
      procedure :: destroy => molecule_destroy
      final :: molecule_finalizer
   end type d3_molecule


   interface
      module subroutine read_molecule_generic(self, unit, format)
         class(d3_molecule), intent(out) :: self
         integer, intent(in) :: unit
         integer, intent(in), optional :: format
      end subroutine read_molecule_generic
   end interface


   interface len
      module procedure :: mol_length
   end interface len

contains


!> constructor for molecular structure
subroutine molecule_new_from_natoms(self, n)
   class(d3_molecule),intent(inout) :: self !< molecular structure information
   integer,intent(in) :: n
   call self%destroy
   self%n = n
   allocate( self%at(n), source = 0 )
   allocate( self%sym(n), source = '  ' )
   allocate( self%xyz(3,n), source = 0.0_wp )
   allocate( self%abc(3,n), source = 0.0_wp )
   allocate( self%dist(n,n), source = 0.0_wp )
   allocate( self%atmass(n), source = 0.0_wp )
   allocate( self%z(n), source = 0.0_wp )
   allocate( self%cn(n), source = 0.0_wp )
end subroutine molecule_new_from_natoms


integer pure elemental function mol_length(self) result(length)
   class(d3_molecule),intent(in) :: self !< molecular structure information
   length = self%n
end function mol_length


subroutine molecule_print_info(self, unit)
   use d3par_atomic_units
   class(d3_molecule), intent(in) :: self
   integer, intent(in) :: unit
   integer :: i
   real(wp) :: conv

   conv = autoaa

   write(unit,'(a)')

   ! atomic coordinates
   write(unit,'(1x,"*",1x,i0,1x,a)') self%n,"atoms in unit cell"
   write(unit,'(a)')
   write(unit,'(5x,"#",3x,"Z",3x,32x,"position/Å",8x,"charge")')
   do i = 1, len(self)
      write(unit,'(i6,1x,i3,1x,a2)',advance='no') i,self%at(i),self%sym(i)
      write(unit,'(3f14.7)',advance='no') self%xyz(:,i)*conv
      write(unit,'(f14.7)') self%z(i)
   enddo
   write(unit,'(a)')

!   ! periodicity
!   write(unit,'(1x,"*",1x,i0,a)') self%npbc,"D periodic system"
!   write(unit,'(a)')
!
!   if (self%npbc > 0) then
!      ! cell parameters
!      write(unit,'(1x,"*",1x,a)') "cell parameter"
!      write(unit,'(a)')
!      write(unit,'(a12,2a15,2x,3a11)') &
!         "|a|/Å", "|b|/Å", "|c|/Å", "α/°", "β/°", "γ/°"
!      write(unit,'(f13.7,2f14.7,1x,3f9.3)') &
!         self%cellpar(:3)*conv,self%cellpar(4:)*180.0_wp/pi
!      write(unit,'(a)')
!
!      ! direct lattice (transformation abc -> xyz)
!      write(unit,'(1x,"*",1x,a)') "direct lattice/Å"
!      write(unit,'(a)')
!      write(unit,'(12x,a,3f14.7)') "a",self%lattice(:,1)*conv
!      write(unit,'(12x,a,3f14.7)') "b",self%lattice(:,2)*conv
!      write(unit,'(12x,a,3f14.7)') "c",self%lattice(:,3)*conv
!      write(unit,'(a)')
!
!      ! reciprocal lattice
!      write(unit,'(1x,"*",1x,a)') "reciprocal lattice/Å⁻¹"
!      write(unit,'(a)')
!      write(unit,'(11x,a,3f14.7)') "a*",self%rec_lat(:,1)/conv
!      write(unit,'(11x,a,3f14.7)') "b*",self%rec_lat(:,2)/conv
!      write(unit,'(11x,a,3f14.7)') "c*",self%rec_lat(:,3)/conv
!      write(unit,'(a)')
!
!      ! geometry in fractional coordinates
!      write(unit,'(1x,"*",1x,a)') "geometry in fractional coordinates"
!      write(unit,'(a)')
!      write(unit,'(5x,"#",3x,"Z",3x,20x,"fractional coordinates")')
!      do i = 1, self%n
!         write(unit,'(i6,1x,i3,1x,a2)',advance='no') i,self%at(i),self%sym(i)
!         write(unit,'(3f14.7)',advance='no') self%abc(:,i)
!         write(unit,'(a)')
!      enddo
!      write(unit,'(a)')
!
!      ! volume of unit cell
!      write(unit,'(1x,"*",1x,a,1x,"=",f14.7)') "volume of direct unit cell/Å³", &
!         self%volume*conv**3
!      write(unit,'(a)')
!   endif
end subroutine molecule_print_info


subroutine molecule_destroy(self)
   class(d3_molecule),intent(inout) :: self !< molecular structure information
   self%n = 0
   self%pbc = .false.
   self%chrg = 0.0_wp
   self%uhf = 0
   self%lattice = 0.0_wp
   if (allocated(self%at)) deallocate(self%at)
   if (allocated(self%sym)) deallocate(self%sym)
   if (allocated(self%xyz)) deallocate(self%xyz)
   if (allocated(self%abc)) deallocate(self%abc)
   if (allocated(self%dist)) deallocate(self%dist)
   if (allocated(self%atmass)) deallocate(self%atmass)
   if (allocated(self%z)) deallocate(self%z)
   if (allocated(self%cn)) deallocate(self%cn)
end subroutine molecule_destroy

subroutine molecule_finalizer(self)
   type(d3_molecule), intent(inout) :: self
   call self%destroy
end subroutine molecule_finalizer


end module d3def_molecule
