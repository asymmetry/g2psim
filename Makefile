###############################################################################
#   This Makefile shows how to compile all C++, C and Fortran
#   files found in $(SRCDIR) directory.
#   Linking is done with g++. Need to have $ROOTSYS defined
###############################################################################

###############################################################################
# System Info
MYOS        := $(shell uname)
ARCH        := $(shell uname -m)
SHELL       := /bin/bash

###############################################################################
# Project Info
LIBFILE     := libG2PSim.so
VERSION     := 2.0.0

SRCDIR      := src
INCDIR      := include
OBJDIR      := obj.$(ARCH)

MODELLIST   := HRSTrans G2PPhys

DICTFILE    := g2psim_dict.cc

MODELLIST   := $(strip $(MODELLIST))
INCDIRS     := $(INCDIR) $(MODELLIST)
LINKDEFFILE := $(subst dict.cc,linkdef.h,$(DICTFILE))
DICTHEADER  := $(addsuffix .h,$(basename $(DICTFILE)))
DICTOBJ     := $(OBJDIR)/$(addsuffix .o,$(basename $(DICTFILE)))

###############################################################################
# Compilers
AR          := ar
CC          := gcc
CXX         := g++
FC          := gfortran
LD          := g++

###############################################################################
# Flags
ifeq ($(ARCH),i686)
    MODE    := -m32
else
    MODE    := -m64
endif
INCFLAGS    := $(addprefix -I,$(INCDIRS))
CFLAGS      := -Wall -fPIC -O3 -g $(MODE) $(INCFLAGS)
CXXFLAGS    := -Wall -fPIC -O3 -g $(MODE) $(INCFLAGS)
FFLAGS      := -Wall -fPIC -O3 -g $(MODE) $(INCFLAGS)
ifeq ($(MYOS),Darwin) #in Darwin, do not use -fno-leading-underscore
    FFLAGS  += -fno-second-underscore -fno-automatic -fbounds-check \
               -funroll-all-loops -fdollar-ok -ffixed-line-length-none \
               -fno-range-check
else
    FFLAGS  += -fno-leading-underscore -fno-second-underscore \
               -fno-automatic -fbounds-check -funroll-all-loops \
               -fdollar-ok -ffixed-line-length-none -fno-range-check
endif
CPPFLAGS    :=
LDFLAGS     := -O3 -g $(MODE)

###############################################################################
# Generate file lists
# source files
FSOURCES    := $(wildcard $(SRCDIR)/*.[Ff])
CSOURCES    := $(wildcard $(SRCDIR)/*.[Cc])
CSOURCES    += $(wildcard $(SRCDIR)/*.[Cc][Cc])
CSOURCES    += $(wildcard $(SRCDIR)/*.[Cc][XxPp][XxPp])
SOURCES     := $(FSOURCES) $(CSOURCES)

# header files
HEADERS     := $(wildcard $(INCDIR)/*.[Hh])
HEADERS     += $(wildcard $(INCDIR)/*.[Hh][Hh])

# obj and dep files
OBJS        := $(addsuffix .o,$(basename $(SOURCES)))
OBJS        := $(patsubst $(SRCDIR)/%.o,$(OBJDIR)/%.o,$(OBJS))
DEPS        := $(subst .o,.d,$(OBJS))

###############################################################################
# Libraries
LIBS        := $(foreach mod,$(MODELLIST),$(join -L -l,$(mod) $(mod)))
SYSLIBS     := -lstdc++ -lm -lgfortran
LIBS        += $(SYSLIBS)

# ROOT
ifdef ROOTSYS
    ROOTCFLAGS := $(shell root-config --cflags)
    ROOTLIBS := $(shell root-config --libs) -lMinuit

    CFLAGS  += $(ROOTCFLAGS)
    CXXFLAGS += $(ROOTCFLAGS)
    LIBS    += $(ROOTLIBS)
else
    $(error $$ROOTSYS environment variable not defined)
endif

# LIBConfig
ifdef LIBCONFIG
    ifneq ($(LIBCONFIG),/usr)
        LCONFCFLAGS += -I$(LIBCONFIG)/include
    endif
    LCONFLIBS := -L$(LIBCONFIG)/lib -lconfig

    INCFLAGS += $(LCONFCFLAGS)
    CFLAGS  += $(LCONFCFLAGS)
    CXXFLAGS += $(LCONFCFLAGS)
    LIBS    += $(LCONFLIBS)
else
    $(error $$LIBCONFIG environment variable not defined)
endif

###############################################################################
# Pre-defined commands
define make-dependency-cc
echo "Making dependency for file $< ......"
set -e; [ -d $(OBJDIR) ] || mkdir -p $(OBJDIR); \
$(CC) -MM $(CPPFLAGS) $(CFLAGS) $< > $@.$$$$; \
sed 's,\($*\)\.o[ :]*,$(OBJDIR)/\1.o $@ : ,g' < $@.$$$$ > $@; \
rm -f $@.$$$$
endef

define make-dependency-cxx
echo "Making dependency for file $< ......"
set -e; [ -d $(OBJDIR) ] || mkdir -p $(OBJDIR); \
$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) $< > $@.$$$$; \
sed 's,\($*\)\.o[ :]*,$(OBJDIR)/\1.o $@ : ,g' < $@.$$$$ > $@; \
rm -f $@.$$$$
endef

define compile-cc
echo "Compiling $< ......"
$(CC) $(CFLAGS) -c $< -o $@
endef

define compile-cxx
echo "Compiling $< ......"
$(CXX) $(CXXFLAGS) -c $< -o $@
endef

define compile-fc
echo "Compiling $< ......"
set -e; [ -d $(OBJDIR) ] || mkdir -p $(OBJDIR); \
$(FC) $(FFLAGS) -c $< -o $@
endef

###############################################################################
# Targets
VPATH       := $(SRCDIR)

.SUFFIXES:
.SUFFIXES: .c .C .cc .CC .cpp .cxx .f .F
.SUFFIXES: .h .hh .H .HH
.SUFFIXES: .o .d

.PHONY: all lib libanalyzer
.PHONY: $(MODELLIST) libanalyzer.so
.PHONY: clean distclean
.PHONY: test help env

###############################################################################
all: lib libanalyzer

lib: $(LIBFILE).$(VERSION)

libanalyzer: libanalyzer.so

###############################################################################
$(LIBFILE).$(VERSION): $(MODELLIST) $(OBJS) $(DICTOBJ)
	@$(LD) -shared $(LDFLAGS) -o $@ $(OBJS) $(DICTOBJ) $(LIBS)
	@ln -sf $@ $(LIBFILE)
	@cd pyg2psim; ln -sf ../$@ $(LIBFILE); cd ..
	@echo "Linking $(LIBFILE) ...... done!"

$(DICTFILE): $(HEADERS) $(LINKDEFFILE)
	@echo "Generating dictionary $@ ......"
	@$(ROOTSYS)/bin/rootcint -f $@ -c $(INCFLAGS) $^

$(DICTOBJ): $(DICTFILE)
	@$(compile-cxx)

libanalyzer.so:
	@$(MAKE) -s -C libanalyzer
	@ln -sf libanalyzer/$@ $@
	@cd pyg2psim; ln -sf ../$@ $@; cd ..
	@echo "Linking $@ ...... done!"

$(MODELLIST):
	@$(MAKE) -s -C $@

###############################################################################
# Make dependencies
$(OBJDIR)/%.d: %.c
	@$(make-dependency-cc)

$(OBJDIR)/%.d: %.C
	@$(make-dependency-cxx)

$(OBJDIR)/%.d: %.cc
	@$(make-dependency-cxx)

$(OBJDIR)/%.d: %.CC
	@$(make-dependency-cxx)

$(OBJDIR)/%.d: %.cpp
	@$(make-dependency-cxx)

$(OBJDIR)/%.d: %.cxx
	@$(make-dependency-cxx)

ifneq ($(DEPS),)
    sinclude $(DEPS)
endif

###############################################################################
# Compile source files
$(OBJDIR)/%.o: %.c
	@$(compile-cc)

$(OBJDIR)/%.o: %.C
	@$(compile-cxx)

$(OBJDIR)/%.o: %.cc
	@$(compile-cxx)

$(OBJDIR)/%.o: %.CC
	@$(compile-cxx)

$(OBJDIR)/%.o: %.cpp
	@$(compile-cxx)

$(OBJDIR)/%.o: %.cxx
	@$(compile-cxx)

$(OBJDIR)/%.o: %.f
	@$(compile-fc)

$(OBJDIR)/%.o: %.F
	@$(compile-fc)

###############################################################################
# Clean
clean:
	@rm -f $(LIBFILE) $(LIBFILE).$(VERSION)
	@rm -f $(DICTFILE) $(DICTHEADER)
	@rm -f libanalyzer.so
	@rm -f $(OBJDIR)/*
	@rm -f *~ *# */*~ */*#

distclean: clean
	@for model in $(MODELLIST); do \
	    if [[ -d $$model ]]; then \
	        make distclean -s -C $$model; \
	    fi; \
	done;
	@make distclean -s -C libanalyzer

test:
	@echo "MYOS:     $(ARCH) $(MYOS)"
	@echo "CFLAGS:   $(CFLAGS)"
	@echo "CXXFLAGS: $(CXXFLAGS)"
	@echo "FFLAGS:   $(FFLAGS)"
	@echo "LDFLAGS:  $(LDFLAGS)"
	@echo "CSOURCES: $(CSOURCES)"
	@echo "FSOURCES: $(FSOURCES)"
	@echo "HEADERS:  $(HEADERS)"
	@echo "OBJS:     $(OBJS)"
	@echo "DEPS:     $(DEPS)"
	@echo "SYSLIBS:  $(SYSLIBS)"
	@echo "LIBS:     $(LIBS)"

help: test

env: test
