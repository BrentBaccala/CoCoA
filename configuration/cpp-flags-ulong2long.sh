#/bin/sh

# Script to choose the best ULong2Long defn (in ULong2Long.H).
# The defn is chosen by setting the CPP flag CoCoA_ULONG2LONG.

# !!! ASSUMES .. is CoCoA_ROOT/ !!!

# Try the three possible settings for the flag CoCoA_ULONG2LONG.
# Print out the first which works (and exit with code 0).
# If none works, print an error message on stderr, and exit with code 99.
# If compilation fails, print error message on stderr, & exit with code 1, 2 or 3.

if [ -z "$CXX" ]
then
  echo "$0: ***ERROR*** \$CXX variable is not defined"  > /dev/stderr
  exit 1
fi

# Create temp dir in which to work...
umask 22
TMP_DIR=CheckULong2Long-$HOSTNAME-$UID-$$
/bin/rm -rf "$TMP_DIR" && mkdir "$TMP_DIR"
if [ $? -ne 0 ]
then
  echo "***ERROR*** $0: unable to create temp directory $TMP_DIR"  > /dev/stderr
  exit 2
fi
cd "$TMP_DIR"
cp ../CheckULong2Long.C .


# Now try the various possible settings for CoCoA_ULONG2LONG
EXIT_CODE=99
for defn in 1 2 3
do
  CPPFLAG="-DCoCoA_ULONG2LONG=$defn"
  $CXX $CXXFLAGS $CPPFLAG -I../../include  CheckULong2Long.C -o CheckULong2Long
  if [ $? -ne 0 ]
  then
    echo "***ERROR*** [[$0]]  Compilation of test program failed" >/dev/stderr
    EXIT_CODE=3
    break
  fi
  ./CheckULong2Long
  if [ $? = 0 ]
  then
    echo "$CPPFLAG"
    EXIT_CODE=0
    break
  fi
done

if [ $EXIT_CODE = 99 ]
then
  echo "*** Failed to determine a working defn of ULong2Long ***"  > /dev/stderr
fi

cd ..
/bin/rm -rf "$TMP_DIR"
exit $EXIT_CODE
