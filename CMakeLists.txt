# CMakeList.txt : CMake project for PG3, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("main")

# Add source to this project's executable.
add_executable (main main.cpp  "dualbase_class.h" "dual1_class.h"    "reverse_mode.cpp")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET main PROPERTY CXX_STANDARD 20)
endif()

# Add the test executable
#add_executable(test_dual_classes tests/test_dual_classes.cpp)

# If your dual*.h files use features from your main library, link it:
# target_link_libraries(test_dual_classes your_main_library)

# Optionally, add the test to CTest (so you can use `ctest`)
#enable_testing()
#add_test(NAME DualClassesTest COMMAND test_dual_classes)

# TODO: Add tests and install targets if needed.
