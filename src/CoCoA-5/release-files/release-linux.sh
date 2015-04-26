#!/bin/bash

#  release-files/release-linux.sh
#--------------------------------------
VER=5.1
MINORMINOR=1  # <--------------- update
#--------------------------------------

# WEBSITE="www.dima.unige.it:/Volumes/WWW/webcocoa"
# WEBSITE="130.251.60.18:/Volumes/WWW/webcocoa"
WEBSITE="WWW/webcocoa"

#--------------------------------------
# CVS directory with latest version (either arg or default value)
if [ $# = 0 ]
then
  COCOALIBDIRWithCVS=`cd ../../..; pwd`
else
  COCOALIBDIRWithCVS=`cd "$1"; pwd`
fi
#--------------------------------------
REL=$VER.$MINORMINOR
# names for release directories
COCOA_TEXT=cocoa-$VER

RELEASE_DIR=$COCOALIBDIRWithCVS/src/CoCoA-5/release-files

FULLPATH_RELEASE_TEXT_DIR=$RELEASE_DIR/$COCOA_TEXT

copy-packages()
{
  echo "Copying packages in    ..${1:(-40)}/";  # last chars
  # rm -rf $1/packages;
  mkdir -p $1/packages;
  /bin/cp -R  ../packages/*  $1/packages/;
  /bin/rm -rf $1/packages/CVS/;
}

copy-CoCoAManual()
{
  echo "Copying CoCoAManual in ..${1:(-40)}/";  # last chars
  # rm -rf $1/CoCoAManual/;
  mkdir -p $1/CoCoAManual;
  /bin/cp -R  ../CoCoAManual/*  $1/CoCoAManual/;
  /bin/mv $1/CoCoAManual/tex/CoCoAManual.pdf $1/CoCoAManual/.;
  /bin/cp ../../../doc/CoCoATranslationTable.html $1/CoCoAManual/.;
  /bin/rm -rf $1/CoCoAManual/CVS/;
  /bin/rm -rf $1/CoCoAManual/aux-files/;
  /bin/rm -rf $1/CoCoAManual/tex/;
  /bin/rm -rf $1/CoCoAManual/html/CVS;
}

copy-emacs()
{
  echo "Copying emacs in       ..${1:(-40)}/";  # last chars
  mkdir -p $1/emacs;
  /bin/cp ../emacs/cocoa5.emacs  $1/emacs/;
  /bin/cp ../emacs/cocoa5.el     $1/emacs/;
}

copy-emacs-el()
{
  echo "Copying emacs.el in    ..${1:(-40)}/";  # last chars
  mkdir -p $1/emacs/;
  /bin/cp ../emacs/cocoa5.el     $1/emacs/;
}

MakeTGZ()
{
  CURR_DIR=`pwd` # just for printing
  echo "compressing  ..${CURR_DIR:(-40)}/$1";
  tar -czf $1-$2.tgz $1;
}

set-mac-icon()
{
    SetFile -a C $2
    DeRez -only icns $1 > MyIcon.rsrc
    if [ -f $2 ]; then    # Destination is a file
	Rez -append MyIcon.rsrc -o $2
	SetFile -a C $2
    elif [ -d $2 ]; then
    # Create the magical Icon\r file
	touch $2/$'Icon\r'
	Rez -append MyIcon.rsrc -o $2/Icon?
	SetFile -a V $2/Icon?
    fi
    rm MyIcon.rsrc
    osascript -e "tell application \"Finder\" to set label index of alias POSIX file \"$2\" to 1"
#    osascript -e "tell application \"Finder\" to set label index of alias POSIX file \"`cd -P -- "$(dirname -- "$2")" && printf '%s\n' "$(pwd -P)/$(basename -- "$2")"`\" to 1"
}

suggest-scp()
{
  echo "scp $1-$2.tgz $WEBSITE/download/bin/cocoa-$REL-$2.tgz"
}

suggest-sftp()
{
  echo "put $1-$2.tgz $WEBSITE/download/bin/cocoa-$REL-$2.tgz"
}

#------------------------------------------------------------
pushd $RELEASE_DIR/   # <------------ always assume we are here
#------------------------------------------------------------
if [ \! -x "$HOME/bin/CoCoAInterpreter-32" ]
then
  echo "$0: ERROR cannot find CoCoAInterpreter-32"
  exit 1
fi
if [ \! -x "$HOME/bin/CoCoAInterpreter-64" ]
then
  echo "$0: ERROR cannot find CoCoAInterpreter-64"
  exit 1
fi

make texdoc
make htmldoc

#------------------------------------------------------------
DATE=`date "+%Y%m%d-%H%M"`
echo " --vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv--"
echo " GENERATING RELEASE from sources in"
echo "   $COCOALIBDIRWithCVS"
echo " HAVE YOU CHECKED OUT AND COMPILED?  If not, do it with either:"
echo "   cd $COCOALIBDIRWithCVS; ./configure --again; make"
echo "   cd $COCOALIBDIRWithCVS; make"
echo " RELEASE DIR(s) will be:"
echo "   $FULLPATH_RELEASE_TEXT_DIR"
echo " --^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^--"
#make -j3 library; make -j3

#------------------------------------------------------------
if [ -d "$FULLPATH_RELEASE_TEXT_DIR" ]
  then
    echo "  renaming --> $COCOA_TEXT-$DATE"
    /bin/mv $FULLPATH_RELEASE_TEXT_DIR $FULLPATH_RELEASE_TEXT_DIR-$DATE
fi
echo "  rm -rf $FULLPATH_RELEASE_TEXT_DIR-2014*"

#------------------------------------------------------------
echo " --======-- CoCoA for linux 32/64 --======--"
echo " --remember to copy the executables in ~/bin: --vvvvvvvvvvvvvvvvv"
echo "cp $COCOALIBDIRWithCVS/src/CoCoA-5/CoCoAInterpreter ~/bin/CoCoAInterpreter-64"
echo " --remember! --^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
copy-packages    $FULLPATH_RELEASE_TEXT_DIR
copy-CoCoAManual $FULLPATH_RELEASE_TEXT_DIR
copy-emacs       $FULLPATH_RELEASE_TEXT_DIR

mkdir -p $FULLPATH_RELEASE_TEXT_DIR/bin
/bin/cp ~/bin/*  $FULLPATH_RELEASE_TEXT_DIR/bin/.
/bin/cp cocoa5-linux  $FULLPATH_RELEASE_TEXT_DIR/cocoa5
/bin/cp ConfigEmacs-linux.sh  $FULLPATH_RELEASE_TEXT_DIR/emacs/ConfigEmacs.sh
chmod +x  $FULLPATH_RELEASE_TEXT_DIR/cocoa5
chmod +x  $FULLPATH_RELEASE_TEXT_DIR/emacs/ConfigEmacs.sh
MakeTGZ "$COCOA_TEXT" "linux"
echo "...done"

#------------------------------------------------------------
echo " --======-- suggest-sftp --======--"
echo "sftp storage1.dima.unige.it"
suggest-sftp      $FULLPATH_RELEASE_TEXT_DIR "linux"
echo "exit"
echo "touch $WEBSITE/download/download5.shtml"
echo " --======-- end --======--"
