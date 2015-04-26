#! /bin/bash

# This script builds the CoCoA-5 GUI; it assumes that CoCoALib
# has already been built, and that the Qt libraries have been installed
# (and that the qmake command is available).

if [ $# -ne 0 ]
then
  echo "ERROR: $0 does not expect any args"           > /dev/stderr
  exit 1
fi

##################################################################
# Check for qmake program

which qmake >/dev/null
if [ $? -ne 0 ]
then
  echo "ERROR: $0 failed because it cannot find qmake" >/dev/stderr
  exit 1
fi

##################################################################
# Check that subdir QCodeEdit/ exists

QCE="QCodeEdit"
if [ \! -d "$QCE" ]
then
  echo "ERROR: $0 cannot find subdirectory $QCE/"
  exit 2
fi


##################################################################
# Set the variable ARCH if we are on a MacOS X system...

SYS_TYPE=`uname`

# On MacOSX we must explicitly state the "architecture".
# The block below should pick the right one (as dictated by libgmp).
if [ "$SYS_TYPE" = "Darwin" ]
then
  ARCH=`arch`  # either "i386" or "ppc"
  if [ "$ARCH" = "i386" ]
  then
    (cd ../../examples; make ex-empty) > /dev/null 2>&1
    arch -x86_64 ../../examples/ex-empty > /dev/null 2>&1
    if [ $? = 0 ]
    then
      ARCH=x86_64
    fi
  fi
  if [ "$ARCH" = "ppc" ]
  then
    (cd ../../examples; make ex-empty) > /dev/null 2>&1
    arch -ppc64 ../../examples/ex-empty > /dev/null 2>&1
    if [ $? = 0 ]
    then
      ARCH=ppc64
    fi
  fi
fi


##################################################################
# Output everything to $QCE/QCodeEdit.pro (deleting old file if nec)

/bin/rm -f $QCE/QCodeEdit.pro

echo "## This file was generated automatically by $0"  > $QCE/QCodeEdit.pro
echo                                                  >> $QCE/QCodeEdit.pro

if [ "$SYS_TYPE" = "Darwin" ]
then
  echo "CONFIG += $ARCH    #MacOSX architecture"      >> $QCE/QCodeEdit.pro
fi

cat $QCE/QCodeEdit.pro.in                             >> $QCE/QCodeEdit.pro


if [ "$SYS_TYPE" = "Darwin" ]
then
  DARWIN="-spec macx-g++"
fi
cd $QCE
qmake $DARWIN QCodeEdit.pro
if [ "$?" -ne 0 ]
then
  echo "ERROR: $0 failed because qmake failed for QCodeEdit"  >/dev/stderr
  exit 4
fi
