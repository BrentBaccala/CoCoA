#! /bin/bash

# Script to check that the BOOST libraries are linkable to
# a program compiled with the GMP compilation flags.
# We check just libboost_filesystem.a, and if that works OK
# then we assume all the BOOST libs are fine too.
# ASSUMES enviroment variables CXX and CXXFLAGS are set correctly.

if [ $# -ne 2 ]
then
  echo "$0: expecting 2 args (extra CXXFLAGS & BOOST_LIBS for linker)" > /dev/stderr
  exit 1
fi

if [ -z "$CXX" ]
then
  echo "$0: expecting environment variables CXX and CXXFLAGS to be set" > /dev/stderr
fi

CXXFLAGS_EXTRA="$1"
BOOST_LDLIBS="$2"

umask 22

# We create a temp dir and work in there; at the end we
# just delete the whole dir -- that should clean up properly.

TMP_DIR=configuration/boost-check-arch-$UID@$HOSTNAME-$$
/bin/rm -rf $TMP_DIR  &&  mkdir $TMP_DIR 2>/dev/null
if [ $? -ne 0 ]
then
  echo "ERROR: $0: Unable to create temporary working directory ($TMP_DIR)"
  exit 2
fi

pushd $TMP_DIR  >/dev/null

# Here is the simple source code we shall use to test the compiler:
cat > TestProg.C <<EOF
#include "boost/filesystem.hpp"
using namespace boost::filesystem;

int main()
{
  path file = "TestProg.C";
  if (!exists(file) || !is_regular_file(file)) exit(2);
}
EOF

echo "$CXX $CXXFLAGS $CXXFLAGS_EXTRA -I../ExternalLibs/include TestProg.C -o TestProg -L../ExternalLibs/lib $BOOST_LDLIBS" > LogFile
$CXX $CXXFLAGS $CXXFLAGS_EXTRA -I../ExternalLibs/include TestProg.C -o TestProg -L../ExternalLibs/lib $BOOST_LDLIBS  2> /dev/null
if [ $? -ne 0 ]
then
  echo "ERROR: $0: compilation failed"
  echo "$CXX $CXXFLAGS $CXXFLAGS_EXTRA TestProg.C -o TestProg -L../ExternalLibs/lib $BOOST_LDLIBS"
#  popd; /bin/rm -rf $TMP_DIR
#  echo "Compilation of CoCoALib's BOOST test prog failed"
  exit 3
fi

./TestProg  2> /dev/null
if [ $? -ne 0 ]
then
  echo "ERROR: $0: TestProg gave run-time error"
#  popd; /bin/rm -rf $TMP_DIR
  exit 4
fi

popd  >/dev/null
/bin/rm -rf $TMP_DIR
