#!/bin/bash

# Script to generate automatically the header for the doxygen documentation
#  with version number and creation date (copied from include/CoCoA/)

# Check we have one arg.
if [ $# -ne 1 ]
then
  echo                                                                  >/dev/stderr
  echo "ERROR: $0 requires 1 arg (the version ID of the CoCoALibrary)." >/dev/stderr
  exit 1
fi
VERSION="$1"

# Check we're in the CoCoALib/include/CoCoA/ directory.  Give error if not.
CWD=`pwd`
LAST=`basename "$CWD"`
TMP=`dirname "$CWD"`
LASTBUTONE=`basename "$TMP"`
if [ "$LAST" != "doxygen" -o "$LASTBUTONE" != "doc" ]
then
  echo                                                                         >/dev/stderr
  echo "ERROR: $0 should be run only in the directory CoCoALib/doc/doxygen" >/dev/stderr
  exit 1
fi

UNIFIEDHDR=settings/html-header.html

# Move existing unified header into a back-up file (with .old suffix).
if [ -f "$UNIFIEDHDR" ]
then
  /bin/rm -rf "$UNIFIEDHDR.old"
  /bin/mv "$UNIFIEDHDR" "$UNIFIEDHDR.old"
fi


echo "<!-- produced automatically by MakeUnifiedHeader.sh -->" >> "$UNIFIEDHDR"
echo "<html>"                                              >> "$UNIFIEDHDR"
echo "<head>"                                              >> "$UNIFIEDHDR"
echo "  <title>CoCoALib -- C++ Library</title>"            >> "$UNIFIEDHDR"
echo "  <link href=\"doxygen-cocoa.css\" rel=\"stylesheet\" type=\"text/css\">"    >> "$UNIFIEDHDR"
echo "</head>"                                             >> "$UNIFIEDHDR"
echo "<body>"                                              >> "$UNIFIEDHDR"
echo "<table bgcolor=#88bbaa width=100% border=1>"         >> "$UNIFIEDHDR"
echo "<tr><td>"                                            >> "$UNIFIEDHDR"
echo "  <table bgcolor=#ccffff width=100% cellpadding=4>"  >> "$UNIFIEDHDR"
echo "  <tr>"                                              >> "$UNIFIEDHDR"
echo "  <td>  <b>CoCoALib-$VERSION</bb>  </td>"            >> "$UNIFIEDHDR"
echo "  <td align=right>  date: "                          >> "$UNIFIEDHDR"
date "+%d %b %Y"                                           >> "$UNIFIEDHDR"
echo "  </td>"                                             >> "$UNIFIEDHDR"
echo "  </tr>"                                             >> "$UNIFIEDHDR"
echo "  </table>"                                          >> "$UNIFIEDHDR"
echo "</td></tr>"                                          >> "$UNIFIEDHDR"
echo "</table>"                                            >> "$UNIFIEDHDR"
echo "<br><br>"                                            >> "$UNIFIEDHDR"
