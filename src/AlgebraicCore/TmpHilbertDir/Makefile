# Makefile for CoCoALib/src/AlgebraicCore/TmpHilbertDir/ directory

COCOA_ROOT=../../..
include $(COCOA_ROOT)/configuration/autoconf.mk
CWD=src/AlgebraicCore/TmpHilbertDir/
DEPEND_FILE=Makefile_dependencies

SRCS=AnnaUtils.C IVectors.C eterms.C poincare.C rum.C unipoly.C TermList.C \
     toric.C \
     TmpPoincareCPP.C

OBJS=$(SRCS:.C=.o)

# Default target, must be the first one!
.PHONY: default
default: $(OBJS)
	ar -rcuS $(COCOA_LIB) $(OBJS)

.PHONY: clean clean-local
clean: clean-local
	@echo "Cleaned CoCoALib/$(CWD)"

clean-local:
	/bin/rm -f $(OBJS)  "$(DEPEND_FILE).old"  ./*~

.PHONY: veryclean
veryclean: clean-local
	/bin/rm -f  "$(DEPEND_FILE)"
	touch "$(DEPEND_FILE)"  # Trick, otherwise make complains after a veryclean
	@echo "Verycleaned CoCoALib/$(CWD)"


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
	@if [ -f "$(DEPEND_FILE)" ] ; \
	then \
	  mv -f "$(DEPEND_FILE)" "$(DEPEND_FILE).old" ; \
	fi
	@echo "Rebuilding dependencies in <CoCoALib>/$(CWD)"
	@echo "##################################################################"   >> $(DEPEND_FILE)
	@echo "# Dependencies for directory $(CWD)"                                  >> $(DEPEND_FILE)
	@echo "#              for CoCoALib version $(VERSION)"                       >> $(DEPEND_FILE)
	@echo "# DO NOT EDIT: this file is automatically generated; to regenerate"   >> $(DEPEND_FILE)
	@echo "# run \`make dependencies' in the CoCoALib root directory."           >> $(DEPEND_FILE)
	@echo "##################################################################"   >> $(DEPEND_FILE)
	@echo                                                                        >> $(DEPEND_FILE)
	for srcfile in $(SRCS); \
	do \
	  objfile=`basename $$srcfile .c`.o; \
	  $(CXX) -I$(COCOA_ROOT)/include $(GMP_INCLUDE) -I. -MM "$$srcfile" -MQ "$$objfile" >> $(DEPEND_FILE); \
	  echo >> $(DEPEND_FILE); \
	done


$(DEPEND_FILE):
	@touch $(DEPEND_FILE)  # to avoid infinite recursion in make!


include $(DEPEND_FILE)

$(OBJS): $(COCOA_ROOT)/configuration/autoconf.mk
