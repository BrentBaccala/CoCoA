# Makefile for CoCoALib/doc/ directory.

CWD=doc/

SRC_DIR=txt
TEX_DIR=tex
TASK_DIR=CoCoALib-tasks
HTML_DIR=html
T2T_DIR=aux-txt2tags
DOC_SRCS = $(wildcard $(SRC_DIR)/*.txt)
DOC_HTML = $(patsubst $(SRC_DIR)/%.txt,$(HTML_DIR)/%.html,$(DOC_SRCS))
TASKS_OUT = CoCoALib-tasks.html CoCoALib-CompletedTasks.html

.SUFFIXES:

# Rule for compiling *.txt files into *.html in HTML_DIR
$(HTML_DIR)/%.html: $(SRC_DIR)/%.txt
	@echo "Generating $*.html"
	@txt2tags -t html -o $@ -i $<  >/dev/null


.PHONY : alldoc
alldoc: htmldoc texdoc tasktables

$(TEX_DIR):
	@mkdir -p $(TEX_DIR)

$(TEX_DIR)/CoCoALib.tex: $(DOC_SRCS) $(TEX_DIR) $(T2T_DIR)/DocTeX.t2t
	@echo "Generating CoCoALib.tex..."
	@txt2tags -o $(TEX_DIR)/CoCoALib.tex $(T2T_DIR)/DocTeX.t2t >/dev/null

$(TEX_DIR)/cocoalib-doc.sty: $(TEX_DIR) $(T2T_DIR)/cocoalib-doc.sty
	@cp -p $(T2T_DIR)/cocoalib-doc.sty $(TEX_DIR)


CoCoALib.pdf: $(DOC_SRCS) $(T2T_DIR)/DocTeX.t2t $(T2T_DIR)/cocoalib-doc.sty
	@which txt2tags >/dev/null 2>&1 ; \
	 if [ $$? -ne 0 ]; \
	 then \
	   echo; \
	   echo "***WARNING***  txt2tags missing ==> not creating/updating PDF documentation"; \
	   echo; \
	   exit; \
	 fi; \
	 make $(TEX_DIR)/CoCoALib.tex $(TEX_DIR)/cocoalib-doc.sty; \
	 echo "Doing LaTeX --> PDF"; \
	 cd $(TEX_DIR); \
	 echo "Starting pass 1 of 3"; pdflatex -halt-on-error CoCoALib.tex >/dev/null; \
	 if [ $$? -ne 0 ]; then echo "!!!LaTeX error!!! -- see file tex/CoCoALib.log"; exit 1; fi; \
	 echo "Starting pass 2 of 3"; pdflatex -halt-on-error CoCoALib.tex >/dev/null; \
	 echo "Starting final pass"; pdflatex -halt-on-error CoCoALib.tex >/dev/null; \
	 mv CoCoALib.pdf ..

.PHONY: texdoc
texdoc: CoCoALib.pdf


.PHONY: htmldoc
htmldoc: $(HTML_DIR)/cocoalib-doc.css $(HTML_DIR)/index.html $(HTML_DIR)/index-alpha.html
	@make -q $(DOC_HTML); \
	 if [ $$? -ne 0 ]; \
	 then \
	   which txt2tags >/dev/null 2>&1; \
	   if [ $$? -ne 0 ]; \
	   then \
	     echo; \
	     echo "***WARNING***  txt2tags missing ==> not creating/updating HTML documentation"; \
	     echo; \
	     exit; \
	   fi; \
	   $(MAKE) -s $(DOC_HTML); \
	 fi


$(HTML_DIR):
	@mkdir -p $(HTML_DIR)

$(HTML_DIR)/cocoalib-doc.css: $(HTML_DIR) $(T2T_DIR)/cocoalib-doc.css
	@cp -p $(T2T_DIR)/cocoalib-doc.css $(HTML_DIR)

$(HTML_DIR)/index-alpha.html: $(HTML_DIR) $(T2T_DIR)/index-alpha.html
	@cp -p $(T2T_DIR)/index-alpha.html $(HTML_DIR)

$(HTML_DIR)/index.html: $(HTML_DIR) $(T2T_DIR)/index.html
	@cp -p $(T2T_DIR)/index.html $(HTML_DIR)

.PHONY: tasktables
tasktables:
	@cd $(TASK_DIR); $(MAKE) -s


.PHONY: clean clean-local clean-subdirs
clean: clean-local clean-subdirs
	@echo "Cleaned CoCoALib/$(CWD)"

clean-local:
	@/bin/rm -f ./*~ ./.*~  ./.\#*
	@/bin/rm -rf $(TEX_DIR)

clean-subdirs:
	@cd $(TASK_DIR); $(MAKE) -s clean
	@cd $(SRC_DIR); /bin/rm -f ./*~  ./.*~  ./.\#*
	@cd $(T2T_DIR); /bin/rm -f ./*~  ./.*~  ./.\#*

.PHONY: veryclean veryclean-local veryclean-subdirs
veryclean: veryclean-local veryclean-subdirs
	@echo "Verycleaned CoCoALib/$(CWD)"

veryclean-local: clean-local
	@/bin/rm -f $(HTML_DIR).tgz  CoCoALib.pdf $(TASKS_OUT)
	@/bin/rm -rf $(HTML_DIR)

veryclean-subdirs:
	@cd $(TASK_DIR); $(MAKE) -s veryclean
	@cd $(SRC_DIR); /bin/rm -f ./*~  ./.*~  ./.\#*
	@cd $(T2T_DIR); /bin/rm -f ./*~  ./.*~  ./.\#*


# Next few lines are for RCS header/log
# $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/doc/Makefile,v 1.32 2014/07/28 14:42:33 abbott Exp $Log
# $Log: Makefile,v $
# Revision 1.32  2014/07/28 14:42:33  abbott
# Summary: Improved the cleaning targets
# Author: JAA
#
# Revision 1.31  2014/03/14 10:59:55  abbott
# Summary: clean and veryclean targets now delete *~ and .*~ files too
# Author: JAA
#
# Revision 1.30  2014/01/29 17:57:05  abbott
# Summary: Rectified silly mistake
# Author: JAA
#
# Revision 1.29  2014/01/29 17:56:09  abbott
# Summary: Improved veryclean-local target
# Author: JAA
#
# Revision 1.28  2013/10/16 14:12:58  abbott
# Added check that pdflatex has run without error; prints mesg & bails if error occurs.
#
# Revision 1.27  2012/10/15 08:52:26  abbott
# Added clean-local to veryclean-local target.
#
# Revision 1.26  2012/05/30 11:00:54  bigatti
# -- added more dependencies
#
# Revision 1.25  2012/02/02 16:02:40  bigatti
# -- removed creation of index-new.html (had been renamed)
#
# Revision 1.24  2011/12/23 14:51:44  bigatti
# -- copying also the two new index pages
#
# Revision 1.23  2011/10/18 12:03:44  abbott
# Redirected stderr when calling "which" (twice)
#
# Revision 1.22  2011/09/23 13:28:21  abbott
# Corrected unified command script for making CoCoALib.pdf
# (removed some "cd" comands which are no longer needed).
#
# Revision 1.21  2011/09/22 16:22:29  abbott
# Minor correction: the CoCoALib.pdf rule is now a single shell script.
#
# Revision 1.20  2011/09/22 15:52:07  abbott
# Considerably changed, and cleaned.
# This Makefile now prints out a "warning" message if some documentation
# needs to be updated but txt2tags is absent.
#
# Revision 1.19  2011/08/31 09:37:18  bigatti
# -- undid last change about target tasktables
#
