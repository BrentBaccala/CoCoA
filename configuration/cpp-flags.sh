#/bin/sh

# Script to determine which CPP flags to set for compiling CoCoALib.
# Calls script cpp-flags-ulong2long.sh

# Move to directory containing this script
cd `dirname "$0"`

if [ -z "$CXX" ]
then
  echo "***ERROR*** [[$0]]  \$CXX variable is not defined"  > /dev/stderr
  exit 1
fi

# Choose the "best" ULong2Long defn
UL2L=`CXX="$CXX"  CXXFLAGS="$CXXFLAGS"  ./cpp-flags-ulong2long.sh`
if [ $? -ne 0 ]
then
  # cpp-flags-ulong2long.sh has already printed error mesg on stderr
  exit 2
fi

echo "$UL2L"
