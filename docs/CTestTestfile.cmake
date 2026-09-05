# CMake generated Testfile for 
# Source directory: /workspace
# Build directory: /workspace/docs
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[unit_tests]=] "/workspace/docs/unit_tests")
set_tests_properties([=[unit_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "/workspace/CMakeLists.txt;154;add_test;/workspace/CMakeLists.txt;0;")
subdirs("lib/progressbar")
