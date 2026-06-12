program tester
  use dftd3_version, only : get_dftd3_version
  implicit none
  character(len=:), allocatable :: version
  call get_dftd3_version(string=version)
  print *, version
end program tester
