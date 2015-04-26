#!/bin/bash

# This script checks that GMP is at least version 4.2 (needed for random number generator).
# If all checks pass, print GMP version number and exit with code 0.
# If GMP is too old, or some other problem arises, print error mesg and return with code 1.

# This script expects two args:
#   first is the full path of the GMP library,
#   second is the full path of the GMP header file.

if [ "$#" != "2" ]
then
  echo "$0: ERROR: requires 2 args (full path of libgmp.a & full path of gmp.h)"
  exit 1
fi

if [ -z "$CXX" ]
then
  echo "$0: ERROR: shell variable CXX must be set to a C++ compiler compatible with GMP"
  exit 1
fi

GMP_LIB="$1"
GMP_HDR="$2"

umask 22

# Check for existence and readability of library and header file.
# Not very readable because "test" has peculiar syntax.
if [ ! \( -f "$GMP_LIB" -a -r "$GMP_LIB" -a -f "$GMP_HDR" -a -r "$GMP_HDR" \) ]
then
  echo "ERROR: GMP header file or library does not exist or is unreadable." > /dev/stderr
  echo "ERROR: Sought header file in \"$GMP_HDR\""                          > /dev/stderr
  echo "ERROR: Sought library in \"$GMP_LIB\""                              > /dev/stderr
  exit 1
fi

# Line below assumes that the GMP header file defines the version macros
# in the order major-minor-patch, and that __GNU_MP_VERSION appears only on those lines.
GMP_HDR_VERSION=`fgrep __GNU_MP_VERSION $GMP_HDR | cut -f 3 -d " " | tr "\n" "."`


# Get version number from the library.  A bit convoluted but should work.
# Create tmp directory, put C prog in it, compile, run, get output, delete directory.
# TMP_DIR depends on hostname, userid, and process number to try to avoid unfortunate
# name clashes if several people try to install CoCoALib simultaneously.
TMP_DIR=/tmp/CoCoALib-$UID@$HOSTNAME-$$
/bin/rm -rf $TMP_DIR && mkdir $TMP_DIR
if [ $? -ne 0 ]; then
  echobox "ERROR: $0: failed to create temporary directory \"$TMP_DIR\"";
  exit 1
fi
pushd $TMP_DIR
cat > TestProg.C <<EOF
#include <stdio.h>

extern const char* const __gmp_version;

int main()
{
  printf("%s\n", __gmp_version);
  return 0;
}
EOF


# Use c++ compiler found by the main configure script to compile $PROG.C
# Need CXXFLAGS too in case it contains important compiler flags (e.g. 32/64 bits)
$CXX $CXXFLAGS TestProg.C -o TestProg "$GMP_LIB"  2> /dev/null

# Check whether compilation failed; if so, complain.
if [ $? -ne 0 ]
then
  echo "ERROR: unable to determine version of GMP library $GMP_LIB"         > /dev/stderr
  echo "ERROR: (compilation failed in gmp-version-check.sh)"                > /dev/stderr
  popd
  /bin/rm -rf $TMP_DIR
  exit 1
fi

# Compilation succeeded, so run $PROG which will print out the version.
GMP_LIB_VERSION=`./TestProg`

# Check whether execution failed; if so, complain (probably linker problems).
if [ $? -ne 0 ]
then
  echo "ERROR: unable to determine version of GMP library $GMP_LIB"            > /dev/stderr
  echo "ERROR: (execution of test program failed in gmp-version-check.sh)"     > /dev/stderr
  echo "ERROR: Check that LD_LIBRARY_PATH contains the GMP library directory." > /dev/stderr
  popd
  /bin/rm -rf $TMP_DIR
  exit 1
fi

popd
/bin/rm -rf $TMP_DIR

# NB in the following line the extra dot after GMP_LIB_VERSION is needed.
if [ "$GMP_HDR_VERSION" != "$GMP_LIB_VERSION." -a "$GMP_HDR_VERSION" != "$GMP_LIB_VERSION.0." ]
then
  echo "ERROR: INCOMPATIBLE library and header files"              > /dev/stderr
  echo "ERROR: library $GMP_LIB is from version: $GMP_LIB_VERSION" > /dev/stderr
  echo "ERROR: header $GMP_HDR is from version: $GMP_HDR_VERSION"  > /dev/stderr
  exit 1
fi

# CoCoALib source assumes existence of some fns which appeared only in GMP 4.2.
if [ "$GMP_LIB_VERSION" \< "4.2.0" ]
then
  echo "ERROR: Your version of GMP is too old: $GMP_LIB_VERSION"   > /dev/stderr
  echo "ERROR: CoCoALib requires version 4.2.0 or newer."          > /dev/stderr
  echo "ERROR: Header file is $GMP_HDR"                            > /dev/stderr
  echo "ERROR: Library file is $GMP_LIB"                           > /dev/stderr
  exit 1
fi

# If we get here, all tests have passed, so print version number and exit with code 0.
echo $GMP_LIB_VERSION
