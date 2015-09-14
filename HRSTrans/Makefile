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
LIBFILE     := libHRSTrans.a
VERSION     := 1.4.2

SRCDIR      := .
INCDIR      := .
OBJDIR      := obj.$(ARCH)

MODELLIST   := HRSTransSTD G2PTrans484816 G2PTrans400016 G2PTrans484816OLD \
               G2PTrans484816R00 G2PTrans400016OLD GDHTransSTD GDHTransLargeX0

MODELLIST   := $(strip $(MODELLIST))
INCDIRS     := $(INCDIR) $(MODELLIST)

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

SUBOBJS     := $(foreach n,$(MODELLIST),$(patsubst $(n)%,$(n)/$(OBJDIR)%,\
	       $(wildcard $(n)/*.[CcFf])))
SUBOBJS     += $(foreach n,$(MODELLIST),$(patsubst $(n)%,$(n)/$(OBJDIR)%,\
               $(wildcard $(n)/*.[Cc][Cc])))
SUBOBJS     += $(foreach n,$(MODELLIST),$(patsubst $(n)%,$(n)/$(OBJDIR)%,\
               $(wildcard $(n)/*.[Cc][XxPp][XxPp])))
SUBOBJS     := $(addsuffix .o, $(basename $(SUBOBJS)))

###############################################################################
# Libraries

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

.PHONY: all lib
.PHONY: $(MODELLIST)
.PHONY: clean distclean
.PHONY: test help env

###############################################################################
all: lib

lib: $(LIBFILE).$(VERSION)

###############################################################################
$(LIBFILE).$(VERSION): $(MODELLIST) $(OBJS)
	@$(AR) cr $(LIBFILE).$(VERSION) $(OBJS) $(SUBOBJS) # &>/dev/null
	@ln -sf $(LIBFILE).$(VERSION) $(LIBFILE)
	@echo "Linking $(LIBFILE) ...... done!"

$(MODELLIST):
	@make -s -C $@

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
	@rm -f $(OBJDIR)/*
	@rm -f *~ *# */*~ */*#

distclean: clean
	@for model in $(MODELLIST); do \
	    if [[ -d $$model ]]; then \
	        make distclean -s -C $$model; \
	    fi; \
	done;

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
	@echo "SUBOBJS:  $(SUBOBJS)"
	@echo "DEPS:     $(DEPS)"

help: test

env: test
