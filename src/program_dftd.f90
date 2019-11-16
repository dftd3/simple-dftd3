program sdftd3_prog
   use d3mod_header
   use d3mod_main
   use d3def_environment
   use d3def_molecule
   implicit none
   type(d3_environment) :: env
   type(d3_molecule) :: mol
   character(len=:), allocatable :: filename

   call env%new

   call read_command_line_arguments(env, filename)
   call env%checkpoint("reading command line arguments")

   call sdftd3_header(env%unit)

   !call main_run(env, filename)
   !call env%checkpoint("running calculation from '"//filename//"'")

contains

subroutine read_command_line_arguments(env, filename)
   use d3def_argparser
   use d3mod_help
   implicit none
   class(d3_environment), intent(inout) :: env
   character(len=:), allocatable, intent(out) :: filename
   character(len=:), allocatable :: string

   type(d3_argparser) :: args
   integer :: iarg
   call args%new

   if (size(args) == 0) then
      call helpmsg(env%unit)
      call env%set_error("no arguments given, so there is nothing to do")
      return
   endif

   ! check for quick help exit first
   if (args%has_option('help')) then
      call helpmsg(env%unit)
      call env%set_finalize("help message printed")
      return
   endif

   ! check for requested version number
   if (args%has_option('version')) then
      call sdftd3_version(env%unit)
      call env%set_finalize("version number printed")
      return
   endif

   ! look up the number of possible input files, we need at least one
   if (file(args) == 0) then
      call helpmsg(env%unit)
      call env%set_error("no file names given, so there is nothing to do")
      return
   endif

   ! now read the arguments


   ! everything left should be a file now, but we only want one input file
   if (file(args) > 1) then
      call env%set_error("multiple file names given, cannot process more than one")
      return
   endif
   call args%get_file(filename)

   ! at this point all command line arguments must be consumed, otherwise...
   if (len(args) > 0) then
      do iarg = 1, len(args)
         call args%get_argument(string)
         call env%set_error("unknown argument: '"//string//"'")
      enddo
   endif

end subroutine read_command_line_arguments

end program sdftd3_prog
