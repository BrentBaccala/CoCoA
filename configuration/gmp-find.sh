#!/bin/bash

# This script looks for a GMP library in a standard location.
# If a single GMP library is found, it prints out the full path of
# the (static) library and returns with exit code 0.
# If none is found or several are found, it prints out an error message
# and returns with a non-zero exit code.

##################################################################
# Use find to search through various standard directories.
# NB look through all directories, even if a GMP has already been found.

# List of directories under which libgmp.a and/or libgmp.so is normally found.
STD_GMP_LIBDIRS="/usr/lib:/usr/lib64:/usr/lib32:/usr/local:/opt/local/lib:/sw/lib:/usr/sfw/lib"
# # Some versions of Linux put libgmp.a in an immediate subdirectory of /usr/lib/
# for file in /usr/lib/*
# do
#   if [ -d "$file" ]
#   then
#     STD_GMP_LIBDIRS="$STD_GMP_LIBDIRS:$file"
#   fi
# done
IFS=":"
LIBGMPPATHS=libgmp-paths
/bin/rm -rf $LIBGMPPATHS
for directory in $STD_GMP_LIBDIRS
do
  if [ -d "$directory" ]
  then
#    if [ -f "$directory/libgmp.a" ];  then echo "$directory/libgmp.a"  >> $LIBGMPPATHS; fi
#    if [ -f "$directory/libgmp.so" ]; then echo "$directory/libgmp.so" >> $LIBGMPPATHS; fi
    find "$directory" -name  libgmp.a   -print >> $LIBGMPPATHS  2> /dev/null
    find "$directory" -name  libgmp.so  -print >> $LIBGMPPATHS  2> /dev/null
  fi
done

if [ \! -s $LIBGMPPATHS ]
then
  # Did not find any plausible GMP installation, so return empty handed.
  echo "No GMP installation found; looked inside $STD_GMP_LIBDIRS"
  /bin/rm $LIBGMPPATHS
  exit 1
fi

# Slightly odd call to wc is to avoid it printing out the file name.
if [ `wc -l < $LIBGMPPATHS` -ne 1 ]
then
  echo "*** Found multiple GMP libraries ***"
  cat $LIBGMPPATHS
  /bin/rm $LIBGMPPATHS
  exit 2
fi

# We have found a single file called libgmp.a or libgmp.so; do a couple of quick
# sanity checks before declaring our search successful...
GMP_LIB=`cat $LIBGMPPATHS`
/bin/rm $LIBGMPPATHS
if [ -f "$GMP_LIB" -a -r "$GMP_LIB" ]
then
  echo "$GMP_LIB"
  exit 0
else
  echo "Trouble reading GMP library file $GMP_LIB"
  exit 4
fi
