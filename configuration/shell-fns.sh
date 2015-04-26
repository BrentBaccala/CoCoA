#!/bin/bash

# Various shell functions used in some of the Makefiles and other
# scripts included with CoCoALib.

# Copyright 2006 John Abbott.
# You are free to use any part of this code in your own programs.


echounderline()
{
  echo "$*"
  echo "$*" | tr "\040-\377" "[-*]"
}

echobox()
{
  mesg="***  $*  ***"
  dashes=`echo "$mesg" | tr "\040-\377" "[-*]"`
  echo "$dashes"
  echo "$mesg"
  echo "$dashes"
}

echoerror()
{
  mesg="*****  $*  *****"
  equals=`echo "$mesg" | tr "\040-\377" "[=*]"`
  echo "$equals"
  echo "$mesg"
  echo "$equals"
}
