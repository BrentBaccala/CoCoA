#!/bin/bash

# This script looks for a BOOST installation in a standard location.
# If a single BOOST installation is found, it prints out the full path of
# the dir containing subdir "boost" and returns with exit code 0.
# If none is found, it prints out an appropriate message (on /dev/stdout)
# and exits with code 1 (used by "configure" in CoCoA root dir).
# If several are found, it prints out an appropriate message (on stdout)
# and returns with exit code 2 (used by "configure" in CoCoA root dir).

# When a single directory is found, only a few basic sanity checks are performed.

##################################################################
# List of most common directories in which boost header dir is normally found.
STD_BOOST_HDR_DIRS="/usr/local/include:/usr/include:/opt/local/include:/sw/include:/usr/sfw/include:"

IFS=":"
TMPFILE=hppboost-paths
/bin/rm -rf $TMPFILE
for dir in $STD_BOOST_HDR_DIRS
do
  if [ -d "$dir" -a -r "$dir" -a -d "$dir/boost" -a -r "$dir/boost" ]
  then
    echo "$dir" >> $TMPFILE
  fi
done

if [ \! -e $TMPFILE ]
then
  # Did not find any plausible BOOST installation, so return empty handed.
  echo "No BOOST installation found; looked inside $STD_BOOST_HDR_DIRS"
  exit 1
fi

# Slightly odd call to wc is to avoid it printing out the file name.
if [ `wc -l < $TMPFILE` -ne 1 ]
then
  echo "*** Found multiple BOOST libraries ***"
  cat $TMPFILE
  /bin/rm $TMPFILE
  exit 2
fi

# We have found a single suitable directory (it exists & is readable,
# and contains readable subdir called "boost"), so return it.
BOOST_HDR_DIR=`cat $TMPFILE`
/bin/rm $TMPFILE
echo "$BOOST_HDR_DIR"
