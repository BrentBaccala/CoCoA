# Makefile for CoCoALib root directory.

# This line will cause make to complain if you run make
# before running configure.
include configuration/autoconf.mk

.PHONY: default
default:
	@echo "======================================================="
	@echo "Compiling CoCoALib-$(VERSION)"
	@echo "  CXX = $(CXX)"
	@echo "  CXXFLAGS = $(CXXFLAGS)"
	@echo "  CXXFLAGS_DEFINES = $(CXXFLAGS_DEFINES)"
	@echo "======================================================="
	@$(MAKE) all


.PHONY: install
install:
	@echo "INSTALL TARGET for CoCoALib NOT YET IMPLEMENTED"
	@exit 1
	$(INSTALL) -d include/CoCoA/*.H $(INSTALL_DIR)/include/CoCoA
	$(INSTALL) -d lib/*             $(INSTALL_DIR)/lib


.PHONY: all
all: check doc examples cocoa5 server


# This target will be "built" only if you try to run make before
# having run configure (which creates autoconf.mk).
configuration/autoconf.mk: configuration/version
	@if [ -f configuration/autoconf.mk ]; \
	 then \
	   echo; \
	   echo "***** ERROR: Version has changed: you must run ./configure --again *****"; \
	   echo; \
	   exit 1; \
	 fi
	@source configuration/shell-fns.sh; \
	 echoerror "ERROR: Cannot build CoCoALib: not yet configured."
	@echo
	@echo "====================================================="
	@echo "***            R U D E   M E S S A G E            ***"
	@echo "====================================================="
	@echo "***  You must run configure before running make.  ***"
	@echo "***  Please read the file README carefully.       ***"
	@echo "====================================================="
	@exit 1


.PHONY: unified-header
unified-header:
	@cd include/CoCoA; $(MAKE) -s

# This target will print out some harmless warning messages if
# you manage to delete the dependencies file in some subdirectory.
.PHONY: dependencies
dependencies: unified-header
	@cd src; $(MAKE) -s dependencies


# Note:  we are not not using "make -C" for compatibility with Solaris make
.PHONY: library
library: dependencies
	@cd src; $(MAKE) -s library

# Just an alias for library
.PHONY: lib
lib: library

# Just an alias for library
.PHONY: cocoalib
cocoalib: library


.PHONY: check
check: library
	@cd src; $(MAKE) -s check


.PHONY: interpreter
interpreter: library
	@cd src/CoCoA-5; $(MAKE) -s CoCoAInterpreter

.PHONY: cocoa5
cocoa5: library
	@cd src/CoCoA-5; \
	 $(MAKE) -s check-prerequisites 2> /dev/null; \
	 if [ $$? = 0 ]; \
	 then \
	   $(MAKE) -s; \
	 else \
	   echo "*************************************************"; \
	   echo "*** Skipping CoCoA-5 because BOOST is missing ***"; \
	   echo "*************************************************"; \
	 fi

.PHONY: server
server: library
	@cd src/server; $(MAKE) -s all


.PHONY: benchmarks
benchmarks: server
	@cd src/server/benchmarks; $(MAKE) -s benchmarks


.PHONY: doc
doc:
	@cd doc; $(MAKE) alldoc


.PHONY: examples
examples: library
	@source configuration/shell-fns.sh; echounderline "Compiling the CoCoALib example programs..."
	@cd examples; $(MAKE) -s executables index.html
	@source configuration/shell-fns.sh; echounderline "Compilation of CoCoALib example programs completed; and index built."



.PHONY: clean clean-local clean-subdirs
clean: clean-local clean-subdirs
	@echo "Cleaned CoCoALib/"

clean-local:
	@/bin/rm -f ./*~ ./.*~  ./.\#*
	@/bin/rm -rf lib/

clean-subdirs:
	@cd src;      $(MAKE) -s clean
	@cd examples; $(MAKE) -s clean
	@cd doc;      $(MAKE) -s clean
	@cd configuration; /bin/rm -f  ./*~  ./.*~  ./.\#*

.PHONY: veryclean veryclean-subdirs
veryclean: clean-local veryclean-subdirs
	@echo "Verycleaned CoCoALib/"

veryclean-subdirs:
	@cd src;       $(MAKE) -s veryclean
	@cd examples;  $(MAKE) -s veryclean
	@cd doc;       $(MAKE) -s veryclean
	@cd include/CoCoA; $(MAKE) -s veryclean
	@cd configuration; /bin/rm -f  ./*~  ./.*~  ./.\#*

.PHONY: distclean
distclean: veryclean
	@cd configuration; /bin/rm -rf autoconf.mk ExternalLibs/


# Next few lines are for RCS header/log
# $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/Makefile,v 1.40 2014/09/02 10:57:27 abbott Exp $
# $Log: Makefile,v $
# Revision 1.40  2014/09/02 10:57:27  abbott
# Summary: Changed distclean so that it removes also configuration/ExternalLibs/
# Author: JAA
#
# Revision 1.39  2014/07/28 14:42:10  abbott
# Summary: Fairly major update: improved "clean", better general organization
# Author: JAA
#
# Revision 1.38  2014/06/13 11:58:48  abbott
# Summary: Now always check whether library.H needs rebuilding
# Author: JAA
#
# Revision 1.37  2014/04/17 12:19:14  abbott
# Summary: Changed trigger for rebuilding dependencies
# Author: JAA
#
# Revision 1.36  2014/03/14 10:59:20  abbott
# Summary: clean and veryclean targets now delete .*~ files too
# Author: JAA
#
# Revision 1.35  2014/03/05 10:42:52  bigatti
# -- undid last cvs
#
# Revision 1.33  2014/01/16 16:08:37  abbott
# Added "aliases" lin and cocoalib for "library" (all 3 do the same thing)
#
# Revision 1.32  2013/06/03 11:58:33  abbott
# Added unified-header as prerequisite for dependencies.
#
# Revision 1.31  2012/10/18 11:18:31  abbott
# Improved error message which is displayed when version has changed.
#
# Revision 1.30  2012/10/01 13:53:11  abbott
# The target "library" now builds unified-header before the dependency files
# (o/w get errors after make veryclean)
#
# Revision 1.29  2012/09/27 15:25:59  abbott
# Added printing of CXXFLAGS_DEFINES.
#
# Revision 1.28  2012/08/05 12:24:00  abbott
# Added visible message when building of CoCoA-5 is skipped.
#
# Revision 1.27  2012/08/02 16:32:23  abbott
# Target cocoa5 now checks prerequisites, and gives helpful message if there's a problem.
#
# Revision 1.26  2012/03/23 16:48:59  abbott
# Changed name of "documentation" phony target into "doc".
#
# Revision 1.25  2011/10/18 20:33:01  abbott
# Revised use of timestamp file for determining when to rebuild dependencies automatically.
#
# Revision 1.24  2011/09/30 12:57:32  bigatti
# -- added cocoa5 to all target
#
# Revision 1.23  2011/09/22 15:51:16  abbott
# Simplified documentation target; now the check for txt2tags is
# performed inside doc/Makefile.
#
# Revision 1.22  2011/08/27 20:42:44  abbott
# Corrected message which is printed when version has changed.
#
# Revision 1.21  2011/07/29 13:58:07  bigatti
# -- generate html and tex doc only after checking for txt2tags
#
# Revision 1.20  2011/07/20 09:04:20  abbott
# New approach to generating Makefile_dependencies: affects almost all Makefiles.
#
# Revision 1.19  2011/07/19 16:26:09  bigatti
# -- changed: "make server" only makes server (and not CoCoA-5)
#
# Revision 1.18  2011/06/01 14:59:34  bigatti
# -- added cocoa5 target
#
# Revision 1.17  2011/05/25 14:31:29  abbott
# Now remakes documentation by default.
#
# Revision 1.16  2011/05/20 16:57:17  abbott
# Corrected extract from autoconf.mk when version has changed
# (consequential to change where configure puts date into autoconf.mk).
#
# Revision 1.15  2011/05/20 16:52:51  abbott
# Added the "interpreter" target.
#
# Revision 1.14  2011/05/03 09:39:24  abbott
# All Makefiles now offer both "clean" and "veryclean" as targets.
#
# Revision 1.13  2011/03/10 17:03:15  bigatti
# -- "all" now makes examples after completing tests
#
# Revision 1.12  2010/10/08 14:11:41  abbott
# Makefile cleaning:
#   (a) proper handling of recursive make,
#   (b) better organized targets (to make embedded shell scripts simpler)
#
# Revision 1.11  2010/10/07 15:41:22  abbott
# Replaced explicit recursive calls to "make" bu calls to "$(MAKE)".
#
# Revision 1.10  2010/06/29 15:10:02  abbott
# Changed "indexes" into "index" following target name change in examples/Makefile.
#
# Revision 1.9  2010/02/09 10:33:57  abbott
# Added missing semicolon (though it seemed to work OK even without it).
#
# Revision 1.8  2009/06/04 17:21:48  abbott
# Make some of the targets a bit simpler (to avoid attempting to compile twice
# the library in some instances).
#
# Revision 1.7  2009/05/21 09:47:59  abbott
# Changed target "server": il now compiles "all" in the directory src/
#
# Revision 1.6  2009/05/20 14:32:32  abbott
# Changed "library" target so that it compiles only the library
# (previously it compiled the server too).
#
# Revision 1.5  2009/03/18 17:04:33  bigatti
# -- when configuration is needed it prints second line of autoconf.mk
#
# Revision 1.4  2008/12/16 15:46:06  abbott
# Changed egrep into fgrep
#
# Revision 1.3  2008/12/16 10:14:13  bigatti
# -- changed makefiles for compatibility with Solaris make (no "-C" option)
#
# Revision 1.2  2007/03/09 18:39:04  abbott
# Removed useless target.
#
# Revision 1.1.1.1  2007/03/09 15:16:10  abbott
# Imported files
#
#
