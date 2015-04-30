# Makefile for CoCoALib/src/AlgebraicCore/ directory.

COCOA_ROOT=../..
include $(COCOA_ROOT)/configuration/autoconf.mk
CWD=src/AlgebraicCore/

DEPEND_FILE=Makefile_dependencies
SUBDIRS=TmpFactorDir TmpHilbertDir

FROBBY_FILES=ExternalLibs-Frobby.C
GSL_FILES=ExternalLibs-GSL.C
NORMALIZ_FILES=ExternalLibs-Normaliz.C

SRCS=time.C assert.C utils.C  utils_gmp.C  ProgressReporter.C  BuildInfo.C ThreadsafeCounter.C \
     error.C MemPool.C GlobalManager.C \
     bool3.C MachineInt.C BigInt.C BigRat.C IntOperations.C \
     convert.C  NumTheory.C  ToString.C  FloatApprox.C \
     random.C \
     RingElemInput.C symbol.C OpenMath.C OpenMathXML.C TeX.C \
     SmallFpImpl.C SmallFpLogImpl.C SmallFpDoubleImpl.C DUPFp.C \
     ring.C ideal.C FieldIdeal.C RingHom.C CanonicalHom.C apply.C \
     RingFp.C RingFpLog.C RingFpDouble.C RingZZ.C RingQQ.C RingTwinFloat.C \
     QuotientRing.C FractionField.C PolyRing.C SparsePolyRing.C \
     degree.C PPOrdering.C  PPMonoid.C  PPMonoidHom.C \
     DivMask.C OrdvArith.C PPWithMask.C \
     PPMonoidEv.C PPMonoidEvOv.C PPMonoidOv.C PPMonoidSparse.C \
     QBGenerator.C  QuasiPoly.C \
     matrix.C  MatrixOperations.C  DenseMatrix.C  MatrixView.C \
     MatrixForOrdering.C  MatrixSpecial.C\
     DistrMPolyClean.C        RingDistrMPolyClean.C \
     DistrMPolyInlPP.C   RingDistrMPolyInlPP.C \
     DistrMPolyInlFpPP.C RingDistrMPolyInlFpPP.C \
     UPoly.C  DenseUPolyRing.C       RingDenseUPolyClean.C   DenseUPolyClean.C \
     GBEnv.C \
     TmpGPoly.C TmpGPair.C TmpGReductor.C TmpGRStats.C TmpGOperations.C \
     SugarDegree.C reduce.C \
     ReductionCog.C geobucket.C \
     module.C FGModule.C FreeModule.C submodule.C QuotientModule.C \
     ModuleOrdering.C \
     ApproxPts.C ApproxPts2.C TmpLDP.C \
     RingWeyl.C \
     IdealOfPoints.C \
     TmpHilbert.C TmpToric.C \
     factor.C  TmpF5.C  TmpF5Mat.C  TmpIsTree.C \
     TmpMonomialIdeal.C TmpPPVector.C \
     TmpMayerVietorisTree.C  TmpMonomialFunctions.C\
     DynamicBitset.C \
     TmpJBDatastructure.C TmpJBSets.C TmpJBEnv.C \
     TmpMorseElement.C  TmpMorsePaths.C  TmpMorseGraph.C  TmpResolutionMinimization.C \
     BenchmarkToolkit.C \
     $(FROBBY_FILES) $(GSL_FILES) $(NORMALIZ_FILES)

# (see below) OBJS depend on $(COCOA_ROOT)/configuration/autoconf.mk
OBJS=$(SRCS:.C=.o)

# The default target compiles debug_new.C and leak_checker.C as
# well as all the files listed in SRCS.
default:
	@echo "*** $(CWD)Makefile: default target ***"
	@cd $(COCOA_ROOT); $(MAKE) library

.PHONY: all
all: $(COCOA_LIB) debug_new.o leak_checker

$(COCOA_LIB): check-source-files $(OBJS) $(SUBDIRS)
	mkdir -p `dirname $(COCOA_LIB)`
	for subdir in $(SUBDIRS); \
	do \
	  echo "Compiling directory $(CWD)$$subdir/"; \
	  touch "$$subdir"; \
	  ( cd "$$subdir"; $(MAKE) -s ); \
	  if [ $$? -ne 0 ]; then echo "***** Compilation of $(CWD)$$subdir/ FAILED *****"; exit 1; fi; \
	done
	echo "doing ar for AlgebraicCore"
	ar -rcuS $(COCOA_LIB) $(OBJS)

$(OBJS):  $(COCOA_ROOT)/configuration/autoconf.mk

.PHONY: check-source-files
check-source-files:
	@./CheckSourceFiles.sh $(SRCS)


#######################################################
# Special compilation targets with extra CPP symbols.
BuildInfo.o: BuildInfo.C $(COCOA_ROOT)/include/CoCoA/BuildInfo.H
	@echo "Compiling BuildInfo.o (with special preprocessor flags)"
	$(COMPILE) -c  -DCOCOA_VERSION="\"$(VERSION)\""  -DCOCOA_CXX="\"$(CXX)\""  -DCOCOA_CXXFLAGS="\"$(CXXFLAGS)\""  -DCOCOA_CXXFLAGS_DEFINES="\"$(CXXFLAGS_DEFINES)\""  -o BuildInfo.o BuildInfo.C

IdealOfPoints.o: IdealOfPoints.C $(COCOA_ROOT)/include/CoCoA/IdealOfPoints.H
	@echo "Compiling IdealOfPoints --  keep your fingers crossed"
	$(COMPILE) -ITmpFactorDir/ -c IdealOfPoints.C

factor.o: factor.C $(COCOA_ROOT)/include/CoCoA/factor.H
	@echo "Compiling factor --  keep your fingers crossed"
	$(COMPILE) -ITmpFactorDir/ -ITmpFactorDir/multivariate/ -c factor.C

TmpHilbert.o: TmpHilbert.C $(COCOA_ROOT)/include/CoCoA/TmpHilbert.H
	@echo "Compiling TmpHilbert --  keep your fingers crossed"
	$(COMPILE) -ITmpHilbertDir/ -c TmpHilbert.C

debug_new.o: debug_new.C $(COCOA_ROOT)/include/CoCoA/debug_new.H
	@echo "Compiling debug_new.o (not part of CoCoALib)"
	$(COMPILE) -c -o debug_new.o debug_new.C

leak_checker: leak_checker.C
	@echo "Compiling leak_checker (not part of CoCoALib)"
	$(COMPILE) -o leak_checker leak_checker.C
	/bin/rm -rf leak_checker.dSYM


.PHONY: clean clean-local clean-subdirs
clean: clean-local clean-subdirs
	@echo "Cleaned CoCoALib/$(CWD)"

clean-local:
	/bin/rm -f $(OBJS) a.out core  "$(DEPEND_FILE).old"
	/bin/rm -f  ./*~  ./.*~  ./.\#*

clean-subdirs:
	for subdir in $(SUBDIRS); do ( cd $$subdir; $(MAKE) -s clean ); done

.PHONY: veryclean veryclean-local veryclean-subdirs
veryclean: veryclean-local veryclean-subdirs
	@echo "Verycleaned CoCoALib/$(CWD)"

veryclean-local: clean-local
	/bin/rm -f  "$(DEPEND_FILE)"  leak_checker  debug_new.o
	touch "$(DEPEND_FILE)"  # Trick, otherwise make complains after a veryclean

veryclean-subdirs:
	for subdir in $(SUBDIRS); do ( cd $$subdir; $(MAKE) -s veryclean ); done


# This target checks whether we should really build DEPEND_FILE
.PHONY: dependencies
dependencies:
	@if [ -f "$(COCOA_HDR)" -a \! "$(DEPEND_FILE)" -nt "$(COCOA_HDR)" ] ; \
	then \
	  $(MAKE) create-dependencies; \
	fi
	for subdir in $(SUBDIRS); do ( cd $$subdir; $(MAKE) -s dependencies ); done

# This is the target which really does rebuild the DEPEND_FILE
.PHONY: create-dependencies
create-dependencies:
	@if [ -f "$(DEPEND_FILE)" ] ; \
	then \
	  mv -f "$(DEPEND_FILE)" "$(DEPEND_FILE).old" ; \
	fi
	@echo "Rebuilding dependencies in <CoCoALib>/$(CWD)"
	@echo "##################################################################"   >> $(DEPEND_FILE)
	@echo "#        ---  CoCoALib version $(VERSION)  ---"                       >> $(DEPEND_FILE)
	@echo "# Dependencies for directory $(CWD)"                                  >> $(DEPEND_FILE)
	@echo "# DO NOT EDIT: this file is automatically generated; to regenerate"   >> $(DEPEND_FILE)
	@echo "# run \`make dependencies' in the CoCoALib root directory."           >> $(DEPEND_FILE)
	@echo "##################################################################"   >> $(DEPEND_FILE)
	@echo                                                                        >> $(DEPEND_FILE)
	for srcfile in $(SRCS); \
	do \
	  objfile=`basename $$srcfile .C`.o; \
	  $(CXX) $(CXXFLAGS) $(CXXFLAGS_DEFINES) $(INCLUDE) $(GMP_INCLUDE) -ITmpFactorDir -ITmpFactorDir/linalg -ITmpFactorDir/multivariate -MM "$$srcfile" -MQ "$$objfile" >> $(DEPEND_FILE); \
	  echo >> $(DEPEND_FILE); \
	done

$(DEPEND_FILE):
	@touch $(DEPEND_FILE)  # to avoid infinite recursion in make!


include $(DEPEND_FILE)


# Next few lines are for RCS header/log
# $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/Makefile,v 1.100 2014/07/28 15:43:50 abbott Exp $
# $Log: Makefile,v $
# Revision 1.100  2014/07/28 15:43:50  abbott
# Summary: Corrected printed message
# Author: JAA
#
# Revision 1.99  2014/07/28 14:47:49  abbott
# Summary: Improved targets (very)clean-local; new scheme for building dependencies
# Author: JAA
#
# Revision 1.98  2014/07/14 10:03:49  abbott
# Summary: Added new QuasiPoly.C file
# Author: JAA
#
# Revision 1.97  2014/06/25 10:35:08  abbott
# Summary: Added missing CXXFLAGS_DEFINES when building dependency file
# Author: JAA
#
# Revision 1.96  2014/06/17 10:08:41  abbott
# Summary: depend_file now depends on library.H (instead of version)
# Author: JAA
#
# Revision 1.95  2014/06/13 12:05:14  abbott
# Summary: Added new (phony) target check-source-files; added two new source files
# Author: JAA
#
# Revision 1.94  2014/05/15 12:32:38  abbott
# Summary: Removed file io.C (no longer needed)
# Author: JAA
#
# Revision 1.93  2014/04/22 13:26:20  abbott
# Summary: Added new file ProgressReporter
# Author: JAA
#
# Revision 1.92  2014/04/17 13:41:16  abbott
# Summary: renamed MatrixViews to MatrixView
# Author: JAA
#
# Revision 1.91  2014/04/11 15:16:32  abbott
# Summary: Renamed TmpFactor.C to factor.C
# Author: JAA
#
# Revision 1.90  2014/04/11 13:33:28  abbott
# Summary: Replaced decimal.C by ToString.C
# Author: JAA
#
# Revision 1.89  2014/04/10 15:45:20  abbott
# Summary: Added FloatApprox.C to sources; updated Makefile_dependencies (& veryclean-local)
# Author: JAA
#
# Revision 1.88  2014/03/14 11:02:24  abbott
# Summary: clean and veryclean targets now delete .*~ files too
# Author: JAA
#
# Revision 1.87  2014/03/07 14:20:33  abbott
# Summary: Added UPoly.C
# Author: JAA
#
# Revision 1.86  2014/03/05 11:13:12  bigatti
# -- undoing wrong commit
#
# Revision 1.84  2014/01/30 15:51:29  bigatti
# -- added RingElemInput.C
#
# Revision 1.83  2013/06/27 16:52:54  abbott
# Added Mario Albert's new resolution/morse impl (still temporary code).
#
# Revision 1.82  2013/06/05 17:45:57  abbott
# Added TmpLDP.C
#
# Revision 1.81  2013/06/03 09:13:35  bigatti
# renamed ModuleTermOrdering into ModuleOrdering
#
# Revision 1.80  2013/05/27 13:09:05  abbott
# Consequential changes from moving all server-related files to src/server/.
#
# Revision 1.79  2013/04/11 15:08:00  abbott
# Major update to IdealOfPoints: first steps towards revitalizing the old BM code.
#
# Revision 1.78  2013/03/15 11:00:50  abbott
# Added check for exponent overflow when powering a PP.
# Merged PPMonoidEv and PPMonoidEvZZ implementations into a single file.
# Implemented new interface for pseudo-ctors for PPMonoidEv which uses a "flag"
# to say whether exponents are big or not.
#
# Revision 1.77  2013/01/23 15:52:05  bigatti
# added BenchmarkToolkit.[CH]
#
# Revision 1.76  2013/01/21 15:30:36  abbott
# Renamed files called BuchbergerMoeller* to IdealOfPoints*.
#
# Revision 1.75  2012/12/04 20:11:18  abbott
# Replaced the 4 source files for random generators, by the single new source file.
#
# Revision 1.74  2012/11/23 17:27:34  abbott
# Added DUPFp and BuchbergerMoeller to standard sources.
#
# Revision 1.73  2012/10/19 16:24:05  bigatti
# -- added      TmpJBDatastructure.C TmpJBSets.C TmpJBEnv.C   (by Mario Albert)
#
# Revision 1.72  2012/10/02 12:56:55  abbott
# Corrected special rule for compiling BuildInfo.
#
# Revision 1.71  2012/10/01 13:58:41  abbott
# Updated command for generating dependency file.
#
# Revision 1.70  2012/05/29 14:56:15  abbott
# Separated ThreadsafeCounter from symbol; also employed it in ring.C.
#
# Revision 1.69  2012/05/28 12:35:28  bigatti
# -- removed personal *.C file
#
# Revision 1.68  2012/05/28 09:18:21  abbott
# Created IntOperations which gathers together all operations on
# integers (both big and small).  Many consequential changes.
#
# Revision 1.67  2012/04/23 09:47:24  abbott
# Added CXXFLAGS and INCLUDE to the command for generating the dependencies file.
#
# Revision 1.66  2012/03/16 15:40:12  abbott
# Merged contents of NumTheoryQQ (continued fraction functions) into NumTheory.
# Merged the doc too.
#
# Revision 1.65  2012/02/10 17:08:44  abbott
# Changed RingZ.C and RingQ.C into RingZZ.C and RingQQ.C
#
# Revision 1.64  2011/11/09 16:27:34  bigatti
# -- renamed MachineInteger --> MachineInt
#
# Revision 1.63  2011/09/23 15:37:24  abbott
# Moved dependency declaration for $(OBJS)
#
# Revision 1.62  2011/09/23 13:20:35  bigatti
# -- QQ.C renamed into BigRat.C
#
# Revision 1.61  2011/08/14 15:52:17  abbott
# Changed ZZ into BigInt (phase 1: just the library sources).
#
# Revision 1.60  2011/08/12 15:21:26  abbott
# Added BigInt impl (derived from ZZ); not used by anyone yet.
#
# Revision 1.59  2011/07/21 09:31:22  bigatti
# -- updated inclusion of NORMALIZ
#
# Revision 1.58  2011/07/20 13:02:32  abbott
# Removed the word "version" from logging message when rebuilding dependencies.
#
# Revision 1.57  2011/07/20 09:55:51  abbott
# Corrected infinite recursion bug (when Makefile_dependencies doesn't exist).
#
# Revision 1.56  2011/07/20 09:04:20  abbott
# New approach to generating Makefile_dependencies: affects almost all Makefiles.
#
# Revision 1.55  2011/06/27 13:44:21  bigatti
# -- added TmpMonomialIdeal
#
# Revision 1.54  2011/06/23 16:07:13  abbott
# Added incomplete but compiling PPMonoidSparse: first prototype,
# simple rather than fast!
#
# Revision 1.53  2011/05/09 14:27:48  bigatti
# -- added TmpToric.C
#
# Revision 1.52  2011/05/03 09:39:24  abbott
# All Makefiles now offer both "clean" and "veryclean" as targets.
#
# Revision 1.51  2011/03/22 20:17:46  abbott
# Removed SmallPrime.C from SRCS.
#
# Revision 1.50  2011/03/21 17:36:39  bigatti
# -- added ExternalLibs-GSL
#
# Revision 1.49  2011/03/01 15:23:43  abbott
# New version to mark addition of MantissaAndExponent, FloatStr.
#
# Revision 1.48  2011/02/10 15:30:59  bigatti
# -- added MatrixSpecial
#
# Revision 1.47  2010/12/17 16:03:15  abbott
# Added new class RandomSource.
#
# Revision 1.46  2010/11/17 15:52:53  abbott
# Removed GMPAllocator.C from list of source files.
#
# Revision 1.45  2010/11/05 14:44:45  bigatti
# -- renamed TmpFrobby.H into ExternalLibs-Frobby.H
# -- renamed TmpNormaliz.H into ExternalLibs-Normaliz.H
#
# Revision 1.44  2010/10/08 22:05:35  abbott
# Removed pointless CoCoALib prefix in definitions of CWD in Makefiles.
#
# Revision 1.43  2010/10/08 14:11:41  abbott
# Makefile cleaning:
#   (a) proper handling of recursive make,
#   (b) better organized targets (to make embedded shell scripts simpler)
#
# Revision 1.42  2010/10/08 13:38:18  bigatti
# -- added TmpNormalizTypes.C
#
# Revision 1.41  2010/10/07 15:41:22  abbott
# Replaced explicit recursive calls to "make" by calls to "$(MAKE)".
#
# Revision 1.40  2010/10/01 16:49:39  bigatti
# -- added NORMALIZ_FILES
#
# Revision 1.39  2010/07/09 17:04:02  abbott
# Added PPMonoid homs.
#
# Revision 1.38  2010/06/29 15:18:08  abbott
# Added RandomZZStream.
#
# Revision 1.37  2010/03/30 15:33:46  bigatti
# -- added DynamicBitset
#
# Revision 1.36  2010/03/23 15:37:00  bigatti
# -- added GBEnv.C
# -- moved "ar" for fasted compilation on Anna's PPC
#
# Revision 1.35  2010/03/23 14:59:55  bigatti
# -- removed RingTwinFloatTZ.C
#
# Revision 1.34  2010/03/18 16:35:08  abbott
# Added NumTheoryQQ files for number theory operations involving QQs.
#
# Revision 1.33  2010/02/16 10:40:09  abbott
# Removed something that shouldn't have been there (well, not yet).
#
# Revision 1.32  2010/02/16 10:19:29  abbott
# Added new class RandomLongStream with example and test.
#
# Revision 1.31  2010/02/03 16:55:19  bigatti
# -- added TmpMonomialFunctions
#
# Revision 1.30  2010/02/03 12:40:15  bigatti
# -- added RegisterServerOpsUser.C
#
# Revision 1.29  2009/12/03 17:26:34  abbott
# Renamed EFGModule to FGModule.
# Renamed ModuleBase member fns  myInit -> myNew, myKill -> myDelete.
# Removed some cruft (old code that was not used by anyone).
#
# Revision 1.28  2009/07/02 16:29:42  abbott
# Added new class QQ to represent rational numbers.
# Consequent change to the Makefile.
#
# Revision 1.27  2009/06/05 12:14:55  abbott
# Major change:
#   created new files NumTheory.H/C  which contain basic number theory operations
#   removed several basic number theory operations from ZZ.H/C
#   removed gcd from MachineInteger.H/C
#   changed names of some basic fns:
#      IsPPrime -> IsProbPrime
#      invmod -> InvMod    (changed signature too)
#      powermod -> PowerMod  (changed signature too)
#   added new fns
#      NextProbPrime & PrevProbPrime
#   consequent changes to other code and tests and examples
#
# Revision 1.26  2009/05/20 14:29:16  abbott
# Moved indexing of CoCoALib into the Makefile of the parent directory (namely src/).
#
# Revision 1.25  2009/05/14 09:47:02  abbott
# Changed "ar" flags to make library construction a bit faster (I hope).
#
# Revision 1.24  2009/02/20 13:27:52  bigatti
# -- renamed sugar.[CH] --> SugarDegree.[CH]
#
# Revision 1.23  2009/02/09 13:59:59  bigatti
# -- new abstract class for computing sugar (tested but not yet used)
#
# Revision 1.22  2009/01/30 17:11:14  bigatti
# -- frobby integration
#
# Revision 1.21  2008/12/16 10:14:13  bigatti
# -- changed makefiles for compatibility with Solaris make (no "-C" option)
#
# Revision 1.20  2008/11/23 18:58:32  abbott
# Major overhaul to preprocessing and SOI/NBM code.
# Split SOI/NBM off into a separate file.
# Preprocessing is now "rational" (but internally guided by C++ doubles).
# SOI/NBM now each have 3 similar interfaces: one purely rational, one for
# input which is represented as doubles, and one which converts the input
# to RingTwinFloat values and produces a result which is over some RingTwinFloat
# (the precision is increased automatically until an answer is obtained).
#
# Revision 1.19  2008/05/30 12:49:43  abbott
# Added an aesthetic space.
#
# Revision 1.18  2008/05/27 16:23:54  bigatti
# -- added TmpPPVector and TmpMayerVietorisTree
#
# Revision 1.17  2008/04/21 11:23:11  abbott
# Separated functions dealing with matrices and PPOrderings into a new file.
# Added matrix norms, and completed adjoint.
#
# Revision 1.16  2008/04/18 15:35:57  abbott
# (long overdue) Major revision to matrices
#
# Revision 1.15  2007/11/20 10:01:26  bigatti
# -- change: TmpF5.C --> TmpF5Mat.C  (by Alberto Arri)
# -- updated and improved test-F5.C
#
# Revision 1.14  2007/10/05 15:28:55  bigatti
# -- added abstract class DenseUPolyRing for representing dense
#    univariate polynomials
# -- added concrete class RingDenseUPolyClean, the cleanest
#    implementation
# -- just for testing, still horribly inefficient and incomplete
#
# Revision 1.13  2007/09/24 14:05:31  abbott
# Added check that subdirectory compilation succeeded; previously ignored subdir errors.
#
# Revision 1.12  2007/06/21 21:29:47  abbott
# Changed name of RingFloat into RingTwinFloat.
#
# Revision 1.11  2007/06/06 15:16:48  abbott
# Added new RandomBitStream class (now based upon GMP's random generator).
# Consequential changes to Makefiles, etc.  There's even doc and an example!
#
# Revision 1.10  2007/05/21 12:57:28  abbott
# New class for passing machine integers as args; includes some simple
# operations on machine integers (cmp, gcd, IsNegative,...).  Operations
# between ZZ and machine integers modified to use the new class.  Inexact
# integer division (of a ZZ) by a negative value now triggers an error;
# new error for reporting inexact integer division by a negative value.
#
# Revision 1.9  2007/05/14 16:36:34  bigatti
# -- removed .C files for Dortmund library (Kaspar)
#
# Revision 1.8  2007/05/03 10:37:13  abbott
# Updated Makefile.
#
# Revision 1.7  2007/04/27 14:58:05  bigatti
# -- content of CoCoAServer.C split into dedicated files
# -- new registration mechanism
#
# Revision 1.6  2007/03/23 18:38:43  abbott
# Separated the "convert" functions (and CheckedCast) into their own files.
# Many consequential changes.  Also corrected conversion to doubles.
#
# Revision 1.5  2007/03/09 20:58:36  abbott
# Final retouches: modified ar flags (again), more important change to main Makefile.
#
# Revision 1.4  2007/03/09 19:09:19  abbott
# Changed flags to ar -- still have trouble with numerical stuff vanishing
# from libcocoa.a when something in AlgebraicCore has to be recompiled.
#
# Revision 1.3  2007/03/09 18:56:56  bigatti
# -- added Tmp prefix to Groebner related files
#
# Revision 1.2  2007/03/09 18:21:56  abbott
# Changed ar flags: use r instead of q (which causes trouble if you partially
# recompile).
#
# Revision 1.1.1.1  2007/03/09 15:16:11  abbott
# Imported files
#
#
