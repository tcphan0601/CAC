#----------------------------------------------------------------
# Generated CMake target import file for configuration "release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "DeePMD::deepmd" for configuration "release"
set_property(TARGET DeePMD::deepmd APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(DeePMD::deepmd PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libdeepmd.so"
  IMPORTED_SONAME_RELEASE "libdeepmd.so"
  )

list(APPEND _cmake_import_check_targets DeePMD::deepmd )
list(APPEND _cmake_import_check_files_for_DeePMD::deepmd "${_IMPORT_PREFIX}/lib/libdeepmd.so" )

# Import target "DeePMD::deepmd_cc" for configuration "release"
set_property(TARGET DeePMD::deepmd_cc APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(DeePMD::deepmd_cc PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELEASE "TensorFlow::tensorflow_cc;TensorFlow::tensorflow_framework"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libdeepmd_cc.so"
  IMPORTED_SONAME_RELEASE "libdeepmd_cc.so"
  )

list(APPEND _cmake_import_check_targets DeePMD::deepmd_cc )
list(APPEND _cmake_import_check_files_for_DeePMD::deepmd_cc "${_IMPORT_PREFIX}/lib/libdeepmd_cc.so" )

# Import target "DeePMD::deepmd_c" for configuration "release"
set_property(TARGET DeePMD::deepmd_c APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(DeePMD::deepmd_c PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELEASE "DeePMD::deepmd_cc"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libdeepmd_c.so"
  IMPORTED_SONAME_RELEASE "libdeepmd_c.so"
  )

list(APPEND _cmake_import_check_targets DeePMD::deepmd_c )
list(APPEND _cmake_import_check_files_for_DeePMD::deepmd_c "${_IMPORT_PREFIX}/lib/libdeepmd_c.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
