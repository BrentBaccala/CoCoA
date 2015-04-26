#! /bin/bash
#----------------------------------------------------------------------
# to configure Emacs for CoCoA run:
#   ./ConfigEmacs.sh
# NB: you will need to run it again if you move the cocoa-5 directory
#----------------------------------------------------------------------

pushd . > /dev/null
SCRIPT_PATH="${BASH_SOURCE[0]}"
while [ -L "$SCRIPT_PATH" ]
do
  cd "`dirname "$SCRIPT_PATH"`"
  SCRIPT_PATH="$(readlink "`basename "$SCRIPT_PATH"`")"
done
cd "`dirname "$SCRIPT_PATH"`" > /dev/null
SCRIPT_DIR="`pwd`"
popd  > /dev/null

COCOA5_EMACS_DIR=$SCRIPT_DIR"/"

FILE=$HOME/.emacs

echo "[1;34m" # to produce coloured text (in some xterms)
echo "configuring emacs.."
echo ";;------[ CoCoA-5 settings ]----------------------"        >> $FILE
echo "(setq cocoa5-emacs-dir"                                    >> $FILE
echo "    \"$COCOA5_EMACS_DIR\")"                                >> $FILE
echo "(if (file-readable-p (concat cocoa5-emacs-dir \"cocoa5.emacs\"))">> $FILE
echo "    (load-file (concat cocoa5-emacs-dir \"cocoa5.emacs\")))"   >> $FILE
echo                                                             >> $FILE
echo ";;-- *Emacs* settings for non-experts --"                  >> $FILE
echo ";;----  MacKeyMode ( command-C, command-V, etc. )"         >> $FILE
echo "(if  (require 'mac-key-mode nil 't)   (mac-key-mode 1) )"  >> $FILE
echo ";;----  TAB for Completion"                                >> $FILE
echo "(defvar cocoa5-tab-is-completion 't)"                      >> $FILE 
echo                                                             >> $FILE
echo ";;------[ CoCoA-5 settings end ]------------------"        >> $FILE
echo "...done"

echo "[1;97;46m"
echo "    The CoCoA-5 settings have been installed successfully:   ";
echo "    They are saved in   $HOME/.emacs                ";
echo "[1;34;49m"
echo "Run ConfigEmacs again if you move the CoCoA folder";
echo "Now you may quit  \"Terminal\"  ";
echo "[0;m"
echo;
echo;
