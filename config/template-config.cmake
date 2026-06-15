@PACKAGE_INIT@

set("SDFTD3_WITH_API" @SDFTD3_WITH_API@)
set("SDFTD3_WITH_OpenMP" @SDFTD3_WITH_OpenMP@)
set("SDFTD3_USE_MCTCLIB" @SDFTD3_USE_MCTCLIB@)

if(NOT TARGET "@PROJECT_NAME@::@PROJECT_NAME@")
  include(CMakeFindDependencyMacro)
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

  if(NOT TARGET "OpenMP::OpenMP_Fortran" AND SDFTD3_WITH_OpenMP)
    find_dependency("OpenMP")
  endif()

  if(NOT TARGET "mctc-lib::mctc-lib" AND SDFTD3_USE_MCTCLIB)
    find_dependency("mctc-lib")
  endif()

  list(REMOVE_AT CMAKE_MODULE_PATH -1)
  include("${CMAKE_CURRENT_LIST_DIR}/@PROJECT_NAME@-targets.cmake")
endif()
