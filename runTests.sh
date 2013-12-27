#!/bin/sh
TEST_LIB="lib/UnitTest++/libUnitTest++.a"
SECAN_SOURCE="lib"
TEST_HEADERS="lib/UnitTest++/src"

if [ ! -f $TEST_LIB ]; then
  echo "UnitTest++ lib not found. Building it..."
  cd lib/UnitTest++
  make all
  cd ../..
fi
g++ -I $SECAN_SOURCE -I $TEST_HEADERS tests/test.cpp \
source/local_alignment.cpp \
source/global_alignment.cpp \
source/local_global_alignment.cpp \
source/uncertain_local_global_alignment.cpp \
source/align_lib.cpp \
$TEST_LIB -o bin/test
echo "---------------\nRunning Tests..."
bin/test
