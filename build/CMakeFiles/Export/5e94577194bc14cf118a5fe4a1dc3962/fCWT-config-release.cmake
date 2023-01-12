#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "fCWT" for configuration "Release"
set_property(TARGET fCWT APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(fCWT PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libfCWT.2.0.dylib"
  IMPORTED_SONAME_RELEASE "@rpath/libfCWT.2.0.dylib"
  )

list(APPEND _cmake_import_check_targets fCWT )
list(APPEND _cmake_import_check_files_for_fCWT "${_IMPORT_PREFIX}/lib/libfCWT.2.0.dylib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
