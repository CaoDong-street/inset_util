# - Config file for the inset_util package
# It defines the following variables
#  INSET_UTIL_INCLUDE_DIRS - include directories for inset_util
#  INSET_UTIL_LIBRARIES    - libraries to link against

# Compute paths
get_filename_component(INSET_UTIL_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
set(INSET_UTIL_INCLUDE_DIRS "${INSET_UTIL_CMAKE_DIR}/../../../include")
set(INSET_UTIL_LIBRARIES "${INSET_UTIL_CMAKE_DIR}/../../../lib")
