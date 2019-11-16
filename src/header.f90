module d3mod_header
   implicit none
   public :: sdftd3_header, sdftd3_version, generic_header
   private

contains

subroutine sdftd3_header(unit)
   integer, intent(in) :: unit
   write(unit, '(32("="),">",1x,6(1x,a),2x,"<",32("="))') &
      & ['s','D','F','T','D','3']
   call sdftd3_version(unit)
end subroutine sdftd3_header

subroutine sdftd3_version(unit)
   integer, intent(in) :: unit
   include 's-dftd3-version.fh'
   write(unit, '(1x,"*",*(1x,a))') &
      & name, "version", version, "compiled by", author, "on", date
end subroutine sdftd3_version

subroutine generic_header(unit, string, width, offset)
   integer, intent(in) :: unit
   integer, intent(in) :: offset
   integer, intent(in) :: width
   character(len=*), intent(in) :: string
   character(len=width) :: dum1,dum2
   character(len=2*width) :: outstring
   character(len=width) :: formatstr
   integer :: strlen, ifront, iback
   strlen = len(string)
   ifront = (width - strlen)/2
   iback  = width - ifront - strlen
   write(dum1,*) width
   write(dum2,*) offset
   write(formatstr,'(i0,"x,a,",i0,"x")') ifront, iback
   write(outstring,'("|",'//formatstr//',"|")') string
   write(unit,'('//dum2//'x,1x,'//dum1//'("-"),1x)')
   write(unit,'('//dum2//'x,a)') trim(outstring)
   write(unit,'('//dum2//'x,1x,'//dum1//'("-"),1x)')
end subroutine generic_header

end module d3mod_header
