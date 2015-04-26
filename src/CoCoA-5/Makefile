# Makefile for CoCoALib/src/CoCoA-5 directory.

CXXFLAGS_CUSTOM=-g

COCOA_ROOT=../..
include $(COCOA_ROOT)/configuration/autoconf.mk
CWD=src/CoCoA-5/

COCOA5_DOC_DIR=CoCoAManual
COCOA5_DOC_XML=$(COCOA5_DOC_DIR)/CoCoAHelp.xml
COCOA5_WORDLIST=$(COCOA5_DOC_DIR)/wordlist.txt
COCOA5_DOC_AUX_DIR=$(COCOA5_DOC_DIR)/aux-files
COCOA5_DOC_HTML_DIR=$(COCOA5_DOC_DIR)/html
COCOA5_DOC_TEX_DIR=$(COCOA5_DOC_DIR)/tex

SRCS=AST.C Lexer.C Main.C Interpreter.C LineProviders.C Parser.C \
     CoCoALibSupplement.C \
     BuiltInFunctions.C BuiltInOneLiners.C \
     BuiltInFunctions-Frobby.C \
     BuiltInFunctions-GSL.C \
     BuiltInFunctions-Normaliz.C \
     OnlineHelp.C \
     VersionInfo.C Banner.C
#     CompilationDate.C


# Default target
.PHONY: all
all: cocoa5 check wordlist
	@if [ $(HAVE_QMAKE) = "no" ]; \
	 then \
	   echo; \
	   echo "*** Skipping building the GUI because \"qmake\" is absent ***"; \
	   echo; \
	 else \
	   $(MAKE) gui; \
	 fi

OBJS=$(SRCS:.C=.o)

$(OBJS): Makefile  $(COCOA_ROOT)/configuration/autoconf.mk

CompilationDate.o: $(OBJS)  CompilationDate.C
	@echo "Compiling CompilationDate"
	$(COMPILE5) -c -o CompilationDate.o CompilationDate.C

COMPILE5=$(COMPILE) $(CXXFLAGS_CUSTOM)

DEPEND_FILE=Makefile_dependencies

# Rule for compiling C++ code in *.C files into *.o object files
%.o: %.C
	@echo "Compiling `basename $@`"
	$(COMPILE5) -c -o $@ $<

# Rule for compiling and linking C++ code in *.C files
%: %.C
	@echo "Compiling `basename $@`"
	$(COMPILE5) -o $@ $< $(LDLIBS)


CoCoAInterpreter: CompilationDate.o  $(COCOA_LIB)
	@echo "CoCoAInterpreter: linking everything together"
	$(COMPILE5) $(OBJS) CompilationDate.o -o CoCoAInterpreter $(LDLIBS)
	@#Next lines remove the *.dSYM directory which the Apple compiler creates
	@#but we remove it only if the compilation flag "-g" is absent.
	@AppleDir=CoCoAInterpreter.dSYM ; \
	echo " " $(CXXFLAGS) " " $(CXXFLAGS_CUSTOM) " " | fgrep " -g " >/dev/null; \
	if [ $$? -eq 1 -a -d "$$AppleDir" ] ; \
	then \
	  /bin/rm -rf "$$AppleDir"; \
	fi


.PHONY: cocoa5
cocoa5: check-prerequisites-verbose
	@$(MAKE) CoCoAInterpreter

.PHONY: check-prerequisites
check-prerequisites:
	@test "$(HAVE_BOOST)" = "yes"

.PHONY: check-prerequisites-verbose
check-prerequisites-verbose:
	@if [ "$(HAVE_BOOST)" = "no" ]; \
	 then \
	   echo "ERROR: Need BOOST libraries to build CoCoAInterpreter (& GUI)"; \
	   echo "ERROR: Reconfigure CoCoA specifying location of BOOST"; \
	   exit 1; \
	 fi


.PHONY: gui
gui: check-prerequisites-verbose $(OBJS) $(COCOA_LIB) QCodeEdit/Makefile C5Makefile
	@echo "Making QCodeEdit"
	cd QCodeEdit; $(MAKE)
	@echo "Making the CoCoA GUI..."
	$(MAKE) -f C5Makefile
	./make-gui-finish.sh

QCodeEdit/Makefile: $(COCOA_LIB) QCodeEdit/QCodeEdit.pro.in make-qcodeeditmakefile.sh
	$(MAKE) veryclean-subdirs
	./make-qcodeeditmakefile.sh

C5Makefile: $(COCOA_LIB)  C5.pro.in  make-c5makefile.sh
	if [ -f C5Makefile ]; then $(MAKE) -f C5Makefile distclean; fi
	./make-c5makefile.sh $(CXXFLAGS_DEFINES) $(LDLIBS)

.PHONY: check
check: cocoa5
	@cd tests; $(MAKE) -s check

.PHONY: ManExamples
ManExamples: cocoa5 OnlineHelp.o AllExamples.o
	$(COMPILE5) OnlineHelp.o AllExamples.o -o PrintManExamples $(LDLIBS)
	/bin/rm -rf PrintManExamples.dSYM # Apple's compiler creates this dir
	@echo "creating ManExamples-in.cocoa5 ...";
	 ./PrintManExamples > $(COCOA_ROOT)/src/CoCoA-5/ManExamples-in.cocoa5
	@echo "...done.  running examples...";
	./CoCoAInterpreter < ManExamples-in.cocoa5;

.PHONY: wordlist
wordlist: OnlineHelp.o wordlist.o
	$(COMPILE5) OnlineHelp.o wordlist.o -o PrintWordlist $(LDLIBS)
	@/bin/rm -rf PrintWordlist.dSYM # Apple's compiler creates this dir
	@echo "creating wordlist..."; ./PrintWordlist > $(COCOA5_WORDLIST)
	@/bin/rm PrintWordlist

wordlist.o: $(COCOA_ROOT)/configuration/autoconf.mk

.PHONY : doc
doc: texdoc htmldoc

.PHONY : htmldoc
htmldoc:
	mkdir -p $(COCOA5_DOC_HTML_DIR)
	/bin/rm -rf $(COCOA5_DOC_HTML_DIR)/*
	cd $(COCOA5_DOC_HTML_DIR); java -jar ../aux-files/saxon.jar  ../../$(COCOA5_DOC_XML) ../aux-files/GUI_help.xsl
#	cp -r $(COCOA5_DOC_AUX_DIR)/GUI-extra-files/* $(COCOA5_DOC_HTML_DIR)/.
#	cd packages/man/generated/; tar cfz doc.tgz doc
#	rm -rf doc
#	cp ../man-xml/generated/doc.tgz .
#	tar xfz doc.tgz

.PHONY : texdoc
texdoc:
	mkdir -p $(COCOA5_DOC_TEX_DIR)
	/bin/rm -rf $(COCOA5_DOC_TEX_DIR)/*
	cd $(COCOA5_DOC_TEX_DIR); java -jar ../aux-files/saxon.jar  ../../$(COCOA5_DOC_XML) ../aux-files/TeX.xsl > CoCoAManual.tex
	cp -r $(COCOA5_DOC_AUX_DIR)/TeX-extra-files/* $(COCOA5_DOC_TEX_DIR)/.
	@cd $(COCOA5_DOC_TEX_DIR); \
	echo "Pass 1 of 3"; pdflatex -halt-on-error CoCoAManual > /dev/null; \
	echo "Pass 2 of 3"; pdflatex -halt-on-error CoCoAManual > /dev/null; \
	echo "Final pass";  pdflatex -halt-on-error CoCoAManual

# OptCoCoAInterpreter: $(SRCS)
# 	$(CXX) -o OptCoCoAInterpreter $(C_FLAGS) $(OPT_FLAGS) $(INCLUDE5) *.C $(LDLIBS)

.PHONY: library
library:
	@(cd $(COCOA_ROOT); $(MAKE) library)

.PHONY: clean clean-local clean-subdirs
clean: clean-local clean-subdirs
	@echo "Cleaned CoCoALib/$(CWD)"

clean-local:
	/bin/rm -rf $(OBJS)  wordlist.o  ./*.dSYM  $(DEPEND_FILE).old  ./*~  ./.*~  ./.\#*
	/bin/rm -f  C5Makefile C5.pro
	/bin/rm -f  moc_C5.cpp qrc_C5.cpp ui_Console.h ui_Debugger.h ui_MainWindow.h ui_SourceEditor.h

clean-subdirs:
	cd $(COCOA5_DOC_DIR); /bin/rm -f ./*~  ./.\#*
	cd emacs; /bin/rm -f ./*~  ./.\#*
	cd packages; /bin/rm -f ./*~  ./.\#*
	cd tests; /bin/rm -f ./*~  ./.\#*
	cd QCodeEdit/; if [ -f "Makefile" ]; then $(MAKE) -s -o QCodeEdit.pro distclean; fi
	/bin/rm -f QCodeEdit/QCodeEdit.pro QCodeEdit/Makefile
	/bin/rm -rf QCodeEdit/.build/
	if [ -d Debug ]; then /bin/rm -rf Debug/; mkdir Debug/; fi

.PHONY: veryclean veryclean-local veryclean-subdirs
veryclean: veryclean-local  veryclean-subdirs
	@/bin/rm -f QCodeEdit/libqcodeedit.a QCodeEdit/Makefile
	@/bin/rm -f ManExamples-in.cocoa5 PLOT-* CoCoAPlot
	@echo "Verycleaned CoCoALib/$(CWD)"

veryclean-local: clean-local
	/bin/rm -rf C5.app Debug
	/bin/rm -f  CoCoAInterpreter  OptCoCoAInterpreter  $(DEPEND_FILE)
	/bin/rm -f  PrintWordlist PrintManExamples
	touch "$(DEPEND_FILE)"  # Trick, otherwise make complains after a veryclean

veryclean-subdirs: clean-subdirs
	/bin/rm -f $(COCOA5_WORDLIST)
	/bin/rm -rf $(COCOA5_DOC_TEX_DIR)
	/bin/rm -rf $(COCOA5_DOC_HTML_DIR)


# This target checks whether we should really build DEPEND_FILE
.PHONY: dependencies
dependencies:
	@if [ -f "$(COCOA_HDR)" -a \! "$(DEPEND_FILE)" -nt "$(COCOA_HDR)" ] ; \
	then \
	  $(MAKE) create-dependencies; \
	fi

# This is the target which really does rebuild the DEPEND_FILE
.PHONY: create-dependencies
create-dependencies:
	if [ -f $(DEPEND_FILE) ] ; \
	then \
	  mv -f $(DEPEND_FILE) $(DEPEND_FILE).old ; \
	fi
	@echo "Rebuilding dependencies in <CoCoALib>/$(CWD)"
	@echo "##################################################################"   >> $(DEPEND_FILE)
	@echo "#        ---  CoCoALib version $(VERSION)  ---"                       >> $(DEPEND_FILE)
	@echo "# Dependencies for directory $(CWD)"                                  >> $(DEPEND_FILE)
	@echo "# DO NOT EDIT: this file is automatically generated; to regenerate"   >> $(DEPEND_FILE)
	@echo "# run \`make dependencies' in the CoCoA-5 directory."                 >> $(DEPEND_FILE)
	@echo "##################################################################"   >> $(DEPEND_FILE)
	@echo                                                                        >> $(DEPEND_FILE)
	if [ $(HAVE_BOOST) = "no" ] ; \
	then \
	  echo "********************"; \
	  echo "*** BOOST absent ***   ==> not rebuilding CoCoA-5 dependencies"; \
	  echo "********************"; \
	  exit 0; \
	fi; \
	for srcfile in $(SRCS); \
	do \
	  objfile=`basename $$srcfile .C`.o; \
	  $(CXX) -I $(COCOA_ROOT)/include $(GMP_INCLUDE) $(BOOST_INCLUDE) -MM "$$srcfile" -MQ "$$objfile" >> $(DEPEND_FILE); \
	  echo >> $(DEPEND_FILE); \
	done


$(DEPEND_FILE):
	@touch $(DEPEND_FILE)  # to avoid infinite recursion in make!


include $(DEPEND_FILE)

# Next few lines are for RCS header/log
# $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/CoCoA-5/Makefile,v 1.97 2014/09/05 16:12:24 abbott Exp $
# $Log: Makefile,v $
# Revision 1.97  2014/09/05 16:12:24  abbott
# Summary: Missed a couple of changes CoCoAInterpreter --> cocoa5; now done
# Author: JAA
#
# Revision 1.96  2014/09/05 16:10:29  abbott
# Summary: Now rebuilds CompilationDate only when necessary; CoCoAInterpreter cleaner; new target cocoa5
# Author: JAA
#
# Revision 1.95  2014/09/01 12:09:27  abbott
# Summary: Strengthened clean-subirs; removed redundancy from veryclean-subdirs
# Author: JAA
#
# Revision 1.94  2014/08/01 11:54:08  bigatti
# -- fixed veryclean
#
# Revision 1.93  2014/08/01 11:51:43  bigatti
# == changed ManExamples (and veryclean)
#
# Revision 1.92  2014/08/01 11:32:02  bigatti
# -- changed ManExamples
#
# Revision 1.91  2014/07/30 16:56:43  bigatti
# ++ added "make doc"
#
# Revision 1.90  2014/07/28 14:55:07  abbott
# Summary: Improved targets (very)clean-local; new scheme for building dependencies
# Author: JAA
#
# Revision 1.89  2014/07/22 11:15:32  bigatti
# -- changed output file for target ManExamples
#
# Revision 1.88  2014/06/28 14:55:26  abbott
# Summary: Revised in line with new configuration scripts
# Author: JAA
#
# Revision 1.87  2014/06/17 10:17:57  abbott
# Summary: depend_file now depends on library.H (instead of version)
# Author: JAA
#
# Revision 1.86  2014/06/13 12:06:44  abbott
# Summary: Corrected name of "CoCoAManual/aux" subdir
# Author: JAA
#
# Revision 1.85  2014/05/14 13:16:21  bigatti
# -- ManExamples nowcreates ManExamples-out
#
# Revision 1.84  2014/05/14 07:37:29  bigatti
# -- commented out line for manual about old gui
#
# Revision 1.83  2014/04/22 12:47:47  abbott
# Summary: Added capability to increase max stack size when starting CoCoAInterpreter
# Author: JAA
#
# Revision 1.82  2014/04/17 12:22:04  abbott
# Summary: Changed clean-local, veryclean-local; new approach for rebuilding dependencies
# Author: JAA
#
# Revision 1.81  2014/04/04 16:41:12  bigatti
# Summary: makefile gui: gmpxx must be last!
#
# Revision 1.80  2014/03/21 12:01:39  abbott
# Summary: Corrected (& simplified) defn of LIBS5
# Author: JAA
#
# Revision 1.79  2014/03/14 11:03:31  abbott
# Summary: clean and veryclean targets now delete *~ and .*~ files too
# Author: JAA
#
# Revision 1.78  2014/03/07 14:17:57  abbott
# Summary: Corrected dependencies for C5Makefile and QCodeEdit/Makefile
# Author: JAA
#
# Revision 1.77  2014/03/05 11:15:23  bigatti
# -- undoing wrong commit
#
# Revision 1.75  2014/01/29 18:55:34  abbott
# Summary: Improved clean-subdirs target; wordlist.txt target now deletes PrintWordList executable.
# Author: JAA
#
# Revision 1.74  2014/01/29 18:00:12  abbott
# Summary: Target veryclean-local now deletes wordlist (generated file)
# Author: JAA
#
# Revision 1.73  2014/01/29 15:46:25  abbott
# Summary: Moved some rm command from veryclean-local to clean-local
# Author: JAA
#
# Revision 1.72  2014/01/28 09:49:23  abbott
# Corrected veryclean-local & veryclean-subdirs targets.
#
# Revision 1.71  2013/06/26 17:37:26  bigatti
# -- added BuiltInFunctions-Frobby.C
#
# Revision 1.70  2013/05/14 14:24:47  abbott
# Revised to be compatible with new BOOST libs script.
#
# Revision 1.69  2013/05/08 09:43:14  abbott
# Added dependency of all object files on the Makefile itself.
#
# Revision 1.68  2013/03/15 17:47:46  abbott
# Further change to make sure Qt dependent code is rebuilt when autoconf.mk changes.
#
# Revision 1.67  2013/03/15 11:03:36  abbott
# Added extra line to force rebuilding of QCodeEdit stuff when autoconf.mk changes.
#
# Revision 1.66  2013/03/12 15:16:57  abbott
# Attempt to fix a "wrong architecture" error message.
#
# Revision 1.65  2013/03/07 17:49:34  abbott
# Added (first impl of) new command SourceRegion.
# Had to hack uglily readAndExecute (in Interpreter.C).
#
# Revision 1.64  2012/10/24 14:05:59  abbott
# Made veryclean-local remove more files.
#
# Revision 1.63  2012/10/02 12:12:40  abbott
# Simplified defn of COMPILE5.
# Removed defn of INCLUDE5.
# Removed commented out code referring to ZIP library.
#
# Revision 1.62  2012/10/01 13:59:30  abbott
# Updated command for generating dependency file.
#
# Revision 1.61  2012/09/27 15:29:09  abbott
# Simplified and updated the rule for building C5Makefile.
#
# Revision 1.60  2012/08/06 09:45:41  abbott
# Removed an unwanted space in a printed message.
#
# Revision 1.59  2012/08/06 09:20:07  abbott
# Improved layout of message saying we're not building GUI.
# Added wordlist.o to files to remove when cleaning.
#
# Revision 1.58  2012/08/06 09:14:53  abbott
# Added a check for HAVE_QMAKE so we give the user a helpful message instead of
# an error when trying to builf GUI but QMAKE is absent.
#
# Revision 1.57  2012/08/05 12:26:03  abbott
# Added check-prerequisites and check-prerequisites-verbose targets
# (replaced the old targets boost-hdr-check and boost-lib-check)
#
# Revision 1.56  2012/08/02 16:35:35  abbott
# Added check-prerequisites target, so main Makefile can give a helpful mesg if there's a problem.
# Added a check for BOOST when rebuilding dependencies -- if absent, the dependency file is left
# (largely) empty.
#
# Revision 1.55  2012/07/19 17:03:23  bigatti
# -- passing GMPXX_LIB as last argument to make-c5makefile
#   (for linker problem with Normaliz)
#
# Revision 1.54  2012/07/16 10:20:24  bigatti
# -- added GSL and Normaliz arguments to "make-c5makefile" call
#
# Revision 1.53  2012/07/11 09:05:06  bigatti
# -- added BuiltInFunctions-Normaliz and BuiltInFunctions-GSL
#
# Revision 1.52  2012/07/05 14:24:25  bigatti
# -- added BuiltInOneLiners.C
#
# Revision 1.51  2012/07/03 12:43:20  abbott
# Added specific path for PrintWordlist
#
# Revision 1.50  2012/06/26 08:52:35  bigatti
# -- make target wordlist less verbose
#
# Revision 1.49  2012/06/26 06:52:49  bigatti
# -- removed experimental COCOA_5_OBJ_DIR
#
# Revision 1.48  2012/06/26 06:44:48  bigatti
# -- added texdoc target
# -- fixed CompilationDate proper compilation also for GUI
# -- added wordlist generation
#
# Revision 1.47  2012/06/18 10:04:22  bigatti
# -- added VersionInfo.C
#
# Revision 1.46  2012/06/07 16:07:05  bigatti
# -- added make for wordlist and htmldoc
#
# Revision 1.45  2012/04/11 09:55:59  abbott
# Makefile now performs boost-lib-check for CoCoAInterpreter
# (previsouly was only boost-hdr-check).
#
# Revision 1.44  2012/04/10 18:55:55  abbott
# Several changes to accomodate using BOOST filesystem fns in Main.C to check
# whether the default packages directory exists.  Needed to modify the script
# which generates C5Makefile as the extra BOOST libs must be passed as params.
#
# Revision 1.43  2012/03/16 17:04:53  abbott
# Corrected clean-subdirs target; it now leaves subdir Debug/ as empty but existing.
#
# Revision 1.42  2012/03/16 16:22:02  abbott
# Improved clean & clean-subdirs targets.
#
# Revision 1.41  2011/10/05 10:00:45  bigatti
# -- fixed make gui for non-standard installation of boost
#
# Revision 1.40  2011/10/03 16:07:57  abbott
# Corrected two minor bugs.
#
# Revision 1.39  2011/09/30 14:29:12  abbott
# Major revision to gui target recipe.
# Removed cruft from start of file.
#
# Revision 1.38  2011/09/30 13:02:33  bigatti
# -- modified target all: no longer fails if boost not present (just
#    prints a big warning)
# -- removed now useless configuration for boost
#
# Revision 1.37  2011/09/28 12:44:52  abbott
# Moved default target up -- correcting a bug I introduced :-(
#
# Revision 1.36  2011/09/21 14:48:11  abbott
# OBJS now depend on autoconf.mk ==> rebuild them if the compilation flags change.
#
# Revision 1.35  2011/09/12 16:49:36  bigatti
# -- added GMP_INC_DIR as 4th arg for make-gui.sh
#
# Revision 1.34  2011/08/17 12:00:48  abbott
# Changed the following targets clean-subdirs, veryclean, veryclean-local
# (cleaning works better in QCodeEdit directory).
#
# Revision 1.33  2011/07/21 17:34:35  abbott
# Consequential change for modified make-gui.sh script.
#
# Revision 1.32  2011/07/20 13:55:06  bigatti
# -- added some "@"
#
# Revision 1.31  2011/07/20 13:02:32  abbott
# Removed the word "version" from logging message when rebuilding dependencies.
#
# Revision 1.30  2011/07/20 09:55:51  abbott
# Corrected infinite recursion bug (when Makefile_dependencies doesn't exist).
#
# Revision 1.29  2011/07/20 09:04:19  abbott
# New approach to generating Makefile_dependencies: affects almost all Makefiles.
#
# Revision 1.28  2011/07/18 09:03:47  abbott
# Removed some pointless (& harmless?) double quotes.
# Added missing ./ before call to make-gui.sh
#
# Revision 1.27  2011/07/16 19:53:43  abbott
# New files and change to Makefile to allow C5 GUI to built automatically
# (as if it were a normal "make" target)
#
# Revision 1.26  2011/07/15 16:56:46  abbott
# Now gets BOOST location from autoconf.mk (i.e. as set by configure script).
#
# Revision 1.25  2011/05/18 15:22:42  bigatti
# -- fixed dependencies on cocoalib
#
# Revision 1.24  2011/05/16 19:02:20  lagorio
# Workaround per gcc bug
#
# Revision 1.23  2011/05/13 17:02:05  abbott
# Moved definition of OBJS; added two .PHONY declarations.
#
# Revision 1.22  2011/05/10 16:45:59  bigatti
# -- added  "make check"
#
# Revision 1.21  2011/05/03 09:39:24  abbott
# All Makefiles now offer both "clean" and "veryclean" as targets.
#
# Revision 1.20  2011/02/15 09:39:34  bigatti
# -- added ManExamples target
#
# Revision 1.19  2011/01/27 15:52:10  bigatti
# -- added -DNDEBUG
#
# Revision 1.18  2011/01/19 16:36:58  bigatti
# -- removed dependency for BuiltInFunctions.o (I didn't do what I meant)
#
# Revision 1.17  2011/01/18 14:44:48  bigatti
# -- added #include CoCoALibSupplement.H
# -- added dependency of BuiltInFunctions from cocoalib.a
#
# Revision 1.16  2010/11/22 17:36:25  bigatti
# -- added "make library"
#
# Revision 1.15  2010/10/11 14:19:27  bigatti
# -- fixed boost path for bigatti
#
# Revision 1.14  2010/10/08 22:05:35  abbott
# Removed pointless CoCoALib prefix in definitions of CWD in Makefiles.
#
# Revision 1.13  2010/10/07 15:41:22  abbott
# Replaced explicit recursive calls to "make" by calls to "$(MAKE)".
#
# Revision 1.12  2010/08/31 14:57:22  bigatti
# -- added OnlineHelp.C to sources
#
# Revision 1.11  2010/08/31 13:14:58  bigatti
# -- commented out compilation of zip library
#
# Revision 1.10  2010/08/23 20:19:22  lagorio
# *** empty log message ***
#
# Revision 1.9  2010/08/02 16:12:03  bigatti
# -- added code for zip support
#
# Revision 1.8  2010/05/07 14:19:01  bigatti
# -- added BuiltInFunctions.C
#
# Revision 1.7  2010/04/22 13:40:14  abbott
# Now inherits most compilation flags from autoconf.mk
# Make clean works properly.
# Object files now depend one the library.
#
# Revision 1.6  2010/04/22 09:30:17  bigatti
# -- added CompilationDate.C and forced its compilation for
#    CoCoAInterpreter target
#
# Revision 1.5  2010/04/21 16:08:06  bigatti
# -- more clever now: recompile only what is needed and with dependencies
#
# Revision 1.4  2010/04/13 15:24:44  bigatti
# -- updated after new "compile" file (with true interpreter starting!)
#
# Revision 1.3  2009/09/13 11:12:06  lagorio
# *** empty log message ***
#
# Revision 1.2  2009/09/13 08:41:13  lagorio
# *** empty log message ***
#
# Revision 1.1  2009/09/03 14:35:44  bigatti
# -- first import
#
