#!/bin/sh
TEST_LIB="lib/UnitTest++/libUnitTest++.a"
TEST_HEADERS="lib/UnitTest++/src"

if [ ! -f $TEST_LIB ]; then
  echo "UnitTest++ lib not found. Building it..."
  cd lib/UnitTest++
  make all
  cd ../..
fi
g++ -I $TEST_HEADERS tests/test.cpp $TEST_LIB -o bin/test
echo "---------------\nRunning Tests..."
bin/test
