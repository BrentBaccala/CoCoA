
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

source ../../configuration/shell-fns.sh

echo
echounderline "Running the CoCoALib tests ($# tests altogether)"

# Keep track of which tests failed, to print a summary at the end.
failures=""

# This loop iterates through the names supplied as arguments to the script.
while [ $# -ne 0 ]
do
  prog="$1"; shift
  /bin/rm -f $prog.cout $prog.cerr
  if [ -f $prog.in ]
  then
    ./$prog < $prog.in > $prog.cout 2> $prog.cerr
  else
    ./$prog > $prog.cout 2> $prog.cerr
  fi
  if [ $? -ne 0 ]
  then
    echo "*****  $prog FAILED  ***** (non-zero exit status)"
    failures="$failures $prog"
  else
    if [ -f $prog.out ]
    then
      diff -w $prog.cout $prog.out > /dev/null
    else
      diff -w $prog.cout /dev/null > /dev/null
    fi
    if [ $? -ne 0 ]
    then
      echo "*****  $prog FAILED  ***** (wrong output)"
      failures="$failures $prog"
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
        failures="$failures $prog"
      else
        /bin/rm $prog.cout $prog.cerr
        echo "$prog..... OK"
      fi
    fi
  fi
done
if [ -z "$failures" ]
then
  echo "===================================="
  echo "Good news: all CoCoALib tests passed"
  echo "===================================="
  echo
  exit 0
fi
if [ "$failures" != " test-RingTwinFloat3 test-RingTwinFloat6" ]
then
  echo "**********************"
  echo "*****  Bad news  *****"
  echo "**********************"
  echo "*****  The following CoCoALib tests failed, please tell us about it."
  echo "*****  $failures"
  exit 1
fi
echo "*****"
echo "*****"
echo "*****  If the only failed tests are test-RingTwinFloat3 and test-RingTwinFloat6"
echo "*****  you probably have a version of GMP older than 4.1.4; updating is"
echo "*****  recommended, but not vital if you do not want to use RingTwinFloat."
echo "*****"
echo "*****"
