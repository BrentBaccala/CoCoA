# Makefile for CoCoALib/examples directory

COCOA_ROOT=..
include $(COCOA_ROOT)/configuration/autoconf.mk
CWD=examples/

CXXFLAGS += --std=c++11

SRCS=ex-empty.C \
      ex-c++.C \
      ex-AlexanderDual.C \
      ex-ApproxPts1.C \
      ex-BigInt1.C ex-BigInt2.C ex-BigInt3.C ex-BigRat1.C \
      ex-bool3.C \
      ex-BuildInfo.C \
      ex-convert1.C \
      ex-DivMask1.C  ex-DivMask2.C \
      ex-DynamicBitset1.C ex-DynamicBitset2.C \
      ex-factor1.C ex-frobby1.C \
      ex-geobucket1.C \
      ex-GMPAllocator1.C ex-GMPAllocator2.C \
      ex-hilbert1.C \
      ex-Janet1.C ex-Janet2.C ex-Janet3.C \
      ex-matrix1.C ex-matrix2.C ex-matrix3.C ex-matrix4.C ex-module1.C ex-module2.C \
      ex-MorseGraph.C \
      ex-MVT.C \
      ex-NF.C ex-Normaliz1.C ex-Normaliz2.C ex-NumTheory1.C ex-NumTheory2.C \
      ex-OrderingGrading1.C \
      ex-PolyInput1.C ex-PolyInput2.C ex-PolyIterator1.C ex-PolyIterator2.C \
      ex-PolyRing1.C ex-PolyRing2.C ex-PolyRing3.C ex-PolyRing4.C \
      ex-PPMonoidElem1.C ex-PPMonoidElem2.C \
      ex-PPMonoidHom1.C \
      ex-PPVector1.C \
      ex-PPWithMask1.C ex-PPWithMask2.C \
      ex-ProgressReporter1.C \
      ex-QuotientBasis.C \
      ex-RandomSource1.C ex-RandomSource2.C \
      ex-RandomBool1.C ex-RandomLong1.C ex-RandomBigInt1.C \
      ex-ToString1.C ex-ToString2.C \
      ex-ring1.C ex-ring2.C \
      ex-RingElem1.C \
      ex-RingFp1.C ex-RingFp2.C \
      ex-RingHom1.C ex-RingHom2.C ex-RingHom3.C ex-RingHom4.C ex-RingHom5.C \
      ex-RingQQ1.C \
      ex-RingTwinFloat1.C ex-RingTwinFloat2.C ex-RingTwinFloat3.C \
      ex-RingWeyl1.C ex-RingWeyl2.C ex-RingWeyl3.C ex-RingWeyl4.C ex-RingWeyl5.C \
      ex-RingZZ1.C \
      ex-SmallFp1.C \
      ex-symbol1.C ex-symbol2.C \
      ex-ToString1.C  ex-ToString2.C \
      ex-UtilsTemplate1.C \
      ex-VectorOperations1.C

ALL_C=$(wildcard *.C)

EXECS=$(SRCS:.C=)
ALL_EXECS=$(ALL_C:.C=)

.PHONY: default
default: 
	@echo "*** CoCoALib/examples/Makefile: default target ***"
	@(cd $(COCOA_ROOT); $(MAKE) examples)

.PHONY: library
library:
	@(cd $(COCOA_ROOT); $(MAKE) library)

# This target should be made only after the CoCoA library has been compiled;
# normally it would be called by the command "make examples" in the CoCoALib root
# directory.
.PHONY: executables
executables: $(EXECS)


.PHONY: check
check: $(EXECS)
	@echo "Running all examples in CoCoALib/examples/.  Please wait."
	@failures=""; \
	 for prog in $(EXECS); \
	 do \
	   /bin/rm -rf errors; \
	   echo "Running $$prog"; \
	   if [ -f "$$prog.in" ]; \
	   then \
	     ./$$prog < $$prog.in > /dev/null 2> errors; \
	   else \
	     ./$$prog < /dev/null > /dev/null 2> errors; \
	   fi; \
	   if [ $$? -ne 0 -o -s errors ]; then failures="$$failures $$prog"; fi; \
	   /bin/rm -rf errors; \
	 done; \
	 source $(COCOA_ROOT)/configuration/shell-fns.sh; \
	 if [ -n "$$failures" ]; then echobox "These examples failed: $$failures"; exit 1; fi; \
	 echobox "Good news: all examples ran successfully."


.PHONY: ValgrindExamples
ValgrindExamples: executables
	@./ValgrindExamples.sh $(EXECS)


.PHONY: clean clean-local
clean: clean-local
	@echo "Cleaned CoCoALib/$(CWD)"

clean-local:
	@/bin/rm -rf $(EXECS) ./*.o core a.out gmon.out
	@/bin/rm -rf ./*.dSYM ./*~ ./.*~ ./.\#*
	@/bin/rm -rf index-tmp.html

.PHONY: veryclean veryclean-local
veryclean: veryclean-local
	@echo "Verycleaned CoCoALib/$(CWD)"

veryclean-local: clean-local
	@/bin/rm -f $(ALL_EXECS)
	@/bin/rm -f CopyInfo index.html


index.html: $(SRCS) CopyInfo index-preamble
	@echo "Creating index.html for directory CoCoALib/examples/"
	@/bin/rm -f index.html index-tmp.html
	@cat index-preamble                                   >> index-tmp.html
	@echo "<table>"                                       >> index-tmp.html;
	@for F in $(SRCS); do \
	  echo "<!--=====================================-->" >> index-tmp.html;\
	  echo "<tr bgcolor=\"#ddddff\">"                     >> index-tmp.html;\
	  echo "<th valign=top><a href=\"$$F\" name=\"$$F\">$$F</a></th>" >> index-tmp.html; \
	  echo "<td>"                                         >> index-tmp.html; \
          ./CopyInfo "ShortDescription" < "$$F"               >> index-tmp.html; \
	  if [ $$? -ne 0 ]; then echo; echo "***ERROR*** Problem with ShortDescription in file $$F"; echo; exit 1; fi; \
	  echo "</td></tr><tr><td colspan=2><pre>"            >> index-tmp.html; \
          ./CopyInfo "LongDescription"  < "$$F"               >> index-tmp.html; \
	  if [ $$? -ne 0 ]; then echo; echo "***ERROR*** Problem with LongDescription in file $$F"; echo; exit 1; fi; \
	  echo "</pre></td></tr>"                             >> index-tmp.html; \
        done
	@echo "</table></body></html>"                        >> index-tmp.html;
	@/bin/mv -f index-tmp.html index.html


$(ALL_EXECS): $(COCOA_LIB)


#############################################################################
# Next few lines are for RCS header/log
# $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/Makefile,v 1.72 2014/07/31 15:59:59 abbott Exp $
# $Log: Makefile,v $
# Revision 1.72  2014/07/31 15:59:59  abbott
# Summary: Replaced ex-io by ex-VectorOperations1
# Author: JAA
#
# Revision 1.71  2014/07/28 14:43:23  abbott
# Summary: Improved the cleaning targets
# Author: JAA
#
# Revision 1.70  2014/07/17 14:49:51  abbott
# Summary: Added ex-matrix4.C
# Author: JAA
#
# Revision 1.69  2014/07/14 15:32:29  abbott
# Summary: Added new example for UtilsTemplate
# Author: JAA
#
# Revision 1.68  2014/07/11 12:12:50  abbott
# Summary: Updated clean-local target
# Author: JAA
#
# Revision 1.67  2014/07/11 11:07:20  abbott
# Summary: Added better error handling to the generation of index.html
# Author: JAA
#
# Revision 1.66  2014/07/03 06:15:50  bigatti
# -- added ex-PPVector1
#
# Revision 1.65  2014/06/13 12:01:53  abbott
# Summary: Changed make variable EXAMPLES into SRCS; replaced phony "index" by real "index.html"
# Author: JAA
#
# Revision 1.64  2014/06/04 10:12:33  bigatti
# -- added ex-ToString1 ex-ToString2
#
# Revision 1.63  2014/05/20 09:02:48  abbott
# Summary: Renamed ex-decimal1 to ex-ToString1; added ex-ToString2
# Author: JAA
#
# Revision 1.62  2014/04/22 14:09:33  abbott
# Summary: Includes new example for ProgressReporter
# Author: JAA
#
# Revision 1.61  2014/04/17 08:39:23  bigatti
# -- added ex-PolyRing4.C (pseudo-constructors)
#
# Revision 1.60  2014/04/09 14:17:13  bigatti
# -- aggiunto ex-DynamicBitset2
#
# Revision 1.59  2014/03/14 11:01:16  abbott
# Summary: clean and veryclean targets now delete .*~ files too
# Author: JAA
#
# Revision 1.58  2014/01/30 15:46:02  bigatti
# -- added ex-PolyInput2
#
# Revision 1.57  2013/06/27 16:46:37  abbott
# Added example for Mario Albert's resolution/morse code.
#
# Revision 1.56  2013/05/31 13:34:28  abbott
# Inserts index-preamble at the start of index.html
#
# Revision 1.55  2013/05/29 14:41:36  abbott
# Added ex-NumTheory2.C (which had been forgotten).
#
# Revision 1.54  2013/05/28 13:29:25  abbott
# Added ex-PolyRing3.C
#
# Revision 1.53  2013/05/28 07:08:54  bigatti
# -- added ex-module2.C
#
# Revision 1.52  2013/05/27 12:55:04  abbott
# Added new example ex-SmallFp1.C
#
# Revision 1.51  2013/03/26 14:54:33  abbott
# Added new example for conversion fns.
#
# Revision 1.50  2012/12/05 11:01:59  abbott
# Renamed the examples for random sequence generators.
#
# Revision 1.49  2012/10/19 16:29:10  bigatti
# -- added ex-Janet 1,2,3 (by Mario Albert)
#
# Revision 1.48  2012/10/05 06:46:21  bigatti
# -- added ex-geobucket1
#
# Revision 1.47  2012/07/04 12:16:22  abbott
# Added ex-matrix3.C
#
# Revision 1.46  2012/05/29 16:34:52  abbott
# Added ex-bool3.C
#
# Revision 1.45  2012/05/11 10:12:42  bigatti
# -- moved ex-C++ as second example (after ex-empty)
#
# Revision 1.44  2012/05/11 10:09:16  bigatti
# -- added ex-c++.C
#
# Revision 1.43  2012/05/10 14:42:17  abbott
# Added new example for anonymous symbols.
#
# Revision 1.42  2012/03/30 10:36:47  abbott
# Renamed ex-BigIntPrime1 to ex-BigInt3.
#
# Revision 1.41  2012/02/10 13:27:01  bigatti
# -- changed ex-RingQ/Z1  --> ex-RingQQ/ZZ1
#
# Revision 1.40  2012/02/03 10:36:26  bigatti
# -- added ex-RandomSource2
#
# Revision 1.39  2011/12/23 15:30:30  bigatti
# -- added ex-matrix2
#
# Revision 1.38  2011/08/26 10:20:11  bigatti
# -- renamed zz->BigInt, QQ->BigRat
#
# Revision 1.37  2011/08/02 13:27:47  abbott
# Added new example about use of decimal conversion fns.
#
# Revision 1.36  2011/07/19 16:21:23  bigatti
# -- added ex-frobby1
#
# Revision 1.35  2011/05/16 13:11:58  abbott
# Target veryclean target now deletes ALL_EXECS (and not just the "official" examples)
#
# Revision 1.34  2011/05/13 17:18:15  abbott
# Fixed clean-local target in Makefile
#
# Revision 1.33  2011/05/03 12:12:04  abbott
# Very minor "cosmetic" changes.
#
# Revision 1.32  2011/05/03 09:39:24  abbott
# All Makefiles now offer both "clean" and "veryclean" as targets.
#
# Revision 1.31  2010/12/17 16:00:49  abbott
# Added new example for RandomSource.
#
# Revision 1.30  2010/10/22 09:13:39  abbott
# Added new GMPAllocator example (and renamed existing one).
# Added a check for output on cerr/clog -- any output means the example "failed".
#
# Revision 1.29  2010/10/08 22:05:36  abbott
# Removed pointless CoCoALib prefix in definitions of CWD in Makefiles.
#
# Revision 1.28  2010/10/08 11:00:31  bigatti
# -- added "wildcard" for convenient compilation of non-official
#    examples (warning: working only with gmake?)
# -- added ex-Normaliz2.C
#
# Revision 1.27  2010/10/07 12:57:36  bigatti
# -- change: make --> $(MAKE)
#
# Revision 1.26  2010/07/16 15:47:34  abbott
# Added simple example for PPMonoidHoms.
#
# Revision 1.25  2010/06/29 15:17:25  abbott
# Added new example for RandomZZStream
#
# Revision 1.24  2010/05/28 16:02:02  bigatti
# -- added ex-DynamicBitset1.C
# -- aligned structure with tests/Makefile
# -- fixed sorting
#
# Revision 1.23  2010/03/30 15:23:23  bigatti
# -- removed indexes from default target
#
# Revision 1.22  2010/03/22 10:17:41  abbott
# Added new example for ZZ: actually renamed ex-ZZ1 to ex-ZZ2, & added new ex-ZZ1
#
# Revision 1.21  2010/03/18 17:36:12  bigatti
# -- added ex-PolyRing2.C
# -- fixed minor mis-sorting
#
# Revision 1.20  2010/03/11 16:33:07  bigatti
# -- fixed internal links in index.html
#
# Revision 1.19  2010/03/11 15:36:49  bigatti
# -- sorted examples alphabetically (except ex-empty.C)
#
# Revision 1.18  2010/03/11 15:32:47  bigatti
# -- new: creation of index.html
#
# Revision 1.17  2010/03/03 10:41:50  abbott
# Added example for basic number theory functions.
#
# Revision 1.16  2010/02/16 10:19:29  abbott
# Added new class RandomLongStream with example and test.
#
# Revision 1.15  2010/02/09 10:34:50  abbott
# Added missing semicolon (though it seemed to work OK even without it).
#
# Revision 1.14  2010/02/04 10:41:43  bigatti
# -- changes: default "make" makes indexes
#
# Revision 1.13  2010/02/01 22:37:43  abbott
# Added new examples for DivMask and PPWithMask
#
# Revision 1.12  2009/10/29 19:10:11  abbott
# Now uses echobox if an error occurs when running the examples.
#
# Revision 1.11  2009/07/08 12:26:53  abbott
# Added floor and ceil functions for QQs.
# Added example program for QQs.
# Minor correction to convert.C; minor cleaning to ex-ZZ1.C
#
# Revision 1.10  2009/07/02 16:22:42  abbott
# Better clean target.
# Added ex-factor1 example to build list.
#
# Revision 1.9  2008/12/16 10:14:13  bigatti
# -- changed makefiles for compatibility with Solaris make (no "-C" option)
#
# Revision 1.8  2008/12/12 11:32:02  abbott
# Updated Makefiles to make the new test/example for symbol visible.
#
# Revision 1.7  2008/11/19 09:21:47  bigatti
# -- added ex-PPMonoidElem2
#
# Revision 1.6  2008/07/16 10:03:54  bigatti
# -- added ex-MVT.C
#
# Revision 1.5  2007/06/21 21:29:47  abbott
# Changed name of RingFloat into RingTwinFloat.
#
# Revision 1.4  2007/06/06 15:41:08  abbott
# Fixed typo in name of ex-RandomBitStream1.C
#
# Revision 1.3  2007/06/06 15:16:48  abbott
# Added new RandomBitStream class (now based upon GMP's random generator).
# Consequential changes to Makefiles, etc.  There's even doc and an example!
#
# Revision 1.2  2007/05/14 16:35:36  bigatti
# -- removed examples for Dortmund library
#
# Revision 1.1.1.1  2007/03/09 15:16:11  abbott
# Imported files
#
