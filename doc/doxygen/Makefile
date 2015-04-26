# makefile for CoCoALib/doc/doxygen directory

COCOA_ROOT=../..
include $(COCOA_ROOT)/configuration/autoconf.mk

.PHONY: default
default: develdoc userdoc

.PHONY: develdoc
develdoc:
	@./MakeUnifiedHeader.sh "$(VERSION)"
	doxygen settings/devel-settings
	@echo "*** output in DeveloperDox/ ***"

.PHONY: userdoc
userdoc:
	@./MakeUnifiedHeader.sh "$(VERSION)"
	doxygen settings/user-settings
	@echo "*** output in UserDox/* ***"

.PHONY: clean
clean:
	/bin/rm -rf DeveloperDox/* UserDox/* DeveloperDox.tgz UserDox.tgz
	/bin/rm -f settings/html-header.html

.PHONY: publish
publish: clean develdoc userdoc
	tar czf DeveloperDox.tgz DeveloperDox
	tar czf UserDox.tgz UserDox
	scp UserDox.tgz DeveloperDox.tgz www:../cocoa/Sites/cocoalib/doc/BleedingEdge/

# for official release
.PHONY: publish-released
publish-released: clean develdoc userdoc
	tar czf DeveloperDox.tgz DeveloperDox
	tar czf UserDox.tgz UserDox
	scp UserDox.tgz DeveloperDox.tgz www:../cocoa/Sites/cocoalib/doc/

# Next few lines are for RCS header/log
# $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/doc/doxygen/Makefile,v 1.4 2014/01/30 08:37:33 abbott Exp $
# $Log: Makefile,v $
# Revision 1.4  2014/01/30 08:37:33  abbott
# Summary: Improved clean target
# Author: JAA
#
# Revision 1.3  2011/05/03 09:39:24  abbott
# All Makefiles now offer both "clean" and "veryclean" as targets.
#
# Revision 1.2  2008/09/19 15:52:09  bigatti
# -- changed: "make clean" removes also DeveloperDox.tgz and UserDox.tgz
#
# Revision 1.1.1.1  2007/03/09 15:16:11  abbott
# Imported files
#
# Revision 1.5  2007/03/08 16:55:06  cocoa
# Changed name of "range" function to "SymbolRange".
#
# Revision 1.4  2006/11/29 16:05:23  cocoa
# -- minor fix
#
# Revision 1.3  2006/11/28 23:03:33  cocoa
# -- header is now generated with version number and date
#
# Revision 1.2  2006/11/27 13:13:14  cocoa
# -- minor changes
#
# Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
# Imported files
#
# Revision 1.1.1.1  2005/10/17 10:46:53  cocoa
# Imported files
#
# Revision 1.3  2005/09/29 08:12:45  cocoa
# -- removed reference to common.mki
#
# Revision 1.2  2005/05/06 14:16:26  cocoa
# -- changed: publish only copies tar.gz files
#
# Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
# Imported files
#
# Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
# Imported files
#
