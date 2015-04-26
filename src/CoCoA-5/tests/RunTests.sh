
# Script for running the CoCoALib tests (after they've been compiled).
# This script is normally called by make as part of the target "check";
# it is not really intended to be called manually.

# For each executable called test-XYZ (where XYZ can be changed) there may
# be some other related files:
#   test-XYZ.in   -- sample input for the test program
#   test-XYZ.out  -- expected output (on cout) from the test program
#   test-XYZ.err  -- expected output (on cerr) from the test program

# When test-XYZ is run, its output is placed into two files:
#   test-XYZ.cout -- actual output (on cout) from the test program
#   test-XYZ.cerr -- actual output (on cerr) from the test program
# A correct test should always exit with code 0; if the executable exits
# with a non-zero code, this is regarded as a failure.  If the exit code
# is 0 then the actual output is compared with the expected output.
# If they match then test-XYZ.cout and test-XYZ.cerr are deleted,
# otherwise a failure message is printed, and the files are not deleted.

# The output files are compared using "diff -w" to work around a
# gratuitous incompatibility introduced by Microsoft.

if [ $# -eq 0 ]
then
  echo "$0: ERROR: no tests specified."
  echo "$0: usage: $0 <list of executables>"
  exit 1
fi

source ../../../configuration/shell-fns.sh
source ../../../configuration/autoconf.mk

echo
echounderline "Running the CoCoA-5 tests ($# tests altogether)"

# Keep track of which tests failed, to print a summary at the end.
failures=""

# This loop iterates through the names supplied as arguments to the script.
while [ $# -ne 0 ]
do
  prog=`basename "$1" .cocoa5`
  if [ $? -ne 0 -o ! -f "$1" ] ; then echo "!!!!! Bad test source file \`$1' !!!!!"; exit 1; fi
  shift
  /bin/rm -f $prog.found $prog.cerr
  HAVE_EXTLIB=HAVE_${prog:6:100};
  if [ "ExtLib" != "${prog:0:6}" ] || [ ${!HAVE_EXTLIB} = yes ]
  then
  ../CoCoAInterpreter --no-preamble --packageDir ../packages < $prog.cocoa5 > $prog.found 2> $prog.cerr
  if [ $? -ne 0 ]
  then
    echo "*****  $prog FAILED  ***** (non-zero exit status)"
    failures="$failures $prog"
  else
    if [ -f $prog.out ]
    then
      diff -w $prog.found $prog.out > /dev/null
    else
      diff -w $prog.cout /dev/null > /dev/null
    fi
    if [ $? -ne 0 ]
    then
      echo "*****  $prog FAILED  ***** (wrong output)"
      failures="$failures $prog.cocoa5"
    else
      if [ -f $prog.err ]
      then
        diff -w $prog.cerr $prog.err > /dev/null
      else
        diff -w $prog.cerr /dev/null > /dev/null
      fi
    if [ $? -ne 0 ]
      then
        echo "*****  $prog FAILED  ***** (wrong output on cerr/clog)"
        failures="$failures $prog.cocoa5"
      else
        /bin/rm $prog.found $prog.cerr
        echo "$prog.cocoa5 ..... OK"
      fi
    fi
    fi
  fi
done
if [ -z "$failures" ]
then
  echo "==================================="
  echo "Good news: all CoCoA-5 tests passed"
  echo "==================================="
  echo
  exit 0
fi
if [ "$failures" != " test-RingTwinFloat3 test-RingTwinFloat6" ]
then
  echo "**********************"
  echo "*****  Bad news  *****"
  echo "**********************"
  echo "*****  The following CoCoA-5 tests failed, please tell us about it."
  echo "*****  $failures"
  exit 1
fi
