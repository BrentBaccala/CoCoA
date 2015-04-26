#! /bin/bash

# Expects 2 args:
CXX="$1"           # first is name of the compiler
CXXFLAGS="$2"      # second contains compilation flags (or empty)

# This script tries to check that CXX is a working C++ compiler,
# and that CXXFLAGS are suitable flags for it.  The check entails
# compiling a very simple source file (with and without CXXFLAGS).

# Check that CXX is an executable file:
FULLCXX=`which "$CXX" 2>/dev/null`
if [ $? -ne 0 -o \! \( -x "$FULLCXX" -a -r "$FULLCXX" -a -f "$FULLCXX" \) ]
then
  echo "Specified compiler \"$CXX\" is not an executable!"
  exit 1
fi

# Now create a temp dir in which to work:
# I use a temp dir so all junk created by the compiler can be removed
# simply by removing the directory -- also I do not need to worry about
# stomping on existing files.
if [ \! -w . ]
then
  echo "$0: ERROR current directory must be writable."
  exit 1
fi

umask 22
TMP_DIR=verify-compiler-$UID@$HOSTNAME-$$
/bin/rm -rf $TMP_DIR  &&  mkdir $TMP_DIR  2>/dev/null
if [ $? -ne 0 ]
then
  echo "Unable to create working directory ($TMP_DIR)"
  exit 1
fi

cd $TMP_DIR

# Here is the simple source code we shall use to test the compiler:
cat > TestProg.C <<EOF
#include <iostream>
using namespace std;
int main()
{
#ifdef __GNUC__
  cout << "gnu";
#else
  cout << "not gnu";
#endif
}
EOF

# Try plain compiler (without CXXFLAGS):
$CXX TestProg.C -o TestProg  > /dev/null 2>&1
if [ $? -ne 0 -o \! -f TestProg -o \! -x TestProg ]
then
  cd ..; /bin/rm -rf $TMP_DIR
  echo "Are you sure \"$CXX\" is a C++ compiler?"
  exit 1
fi
/bin/rm TestProg  # not necessary, just being tidy :-)

# Try compiler with CXXFLAGS:
$CXX $CXXFLAGS TestProg.C -o TestProg  > /dev/null 2>&1
if [ $? -ne 0 -o \! -f TestProg -o \! -x TestProg ]
then
  cd ..; /bin/rm -rf $TMP_DIR
  echo "Compilation flags \"$CXXFLAGS\" seem to be unsuitable for \"$CXX\""
  exit 1
fi

COMPILER_TYPE=`./TestProg`
cd ..; /bin/rm -rf $TMP_DIR

echo "$COMPILER_TYPE"
