#! /bin/sh

if [ "$#" != "1" ]
then
  echo "ERROR: $0 requires 1 arg (full path of gmp.h)" > /dev/stderr
  exit 1
fi

# Assume the file exists and is readable -- already checked in configure script.
GMP_H="$1"
GMP_INC_DIR=`dirname "$GMP_H"`

# Below we create a small C++ program for printing out the GMP version number.
# We use two temporary files.
TMPFILE=get-gmp-version-$UID@$HOSTNAME-$$
/bin/rm -rf $TMPFILE $TMPFILE.C
cat > $TMPFILE.C <<EOF
#include "gmp.h"
#include <iostream>

int main()
{
  std::cout << __GNU_MP_VERSION << "." << __GNU_MP_VERSION_MINOR << "." << __GNU_MP_VERSION_PATCHLEVEL << std::endl;
}
EOF

$CXX -I "$GMP_INC_DIR" $TMPFILE.C -o $TMPFILE  2> /dev/null
if [ $? -ne 0 ]
then
  /bin/rm -f $TMPFILE.C $TMPFILE
  echo "ERROR:  Failed to find GMP version number; perhaps your GMP is too old?"    > /dev/stderr
  echo "ERROR:  Are you sure that \"$GMP_H\" really is a recent GMP header file?"   > /dev/stderr
  exit 1
fi
GMP_VER=`./$TMPFILE`
/bin/rm -f $TMPFILE.C $TMPFILE

# CoCoALib source assumes existence of some fns which appeared only in GMP 4.2.
if [ "$GMP_VER" \< "4.2.1" ]
then
  echo "ERROR: Your version of GMP is too old: you have $GMP_VER"  > /dev/stderr
  echo "ERROR: but CoCoALib requires version 4.2.1 or newer."      > /dev/stderr
  echo "ERROR: Header file is $GMP_HDR"                            > /dev/stderr
  echo "ERROR: Library file is $GMP_LIB"                           > /dev/stderr
  echo "ERROR: The latest GMP can be found at http://gmplib.org/"  > /dev/stderr
  exit 1
fi

echo $GMP_VER
