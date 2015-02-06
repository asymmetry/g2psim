###############################################################################
#   This Makefile shows how to compile all C++, C and Fortran
#   files found in $(SRCDIR) directory.
#   Linking is done with g++. Need to have $ROOTSYS defined
###############################################################################

###############################################################################
# System Info
MYOS        := $(shell uname)
ARCH        := $(shell uname -m)
SHELL       := /bin/sh

###############################################################################
# Project Info (only objs)
SRCDIR      := .
INCDIR      := ..
OBJDIR      := obj.$(ARCH)

INCDIRS     := $(INCDIR)

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
CFLAGS      := -Wall -fPIC -Wno-unused-variable -O3 -g $(MODE) $(INCFLAGS)
CXXFLAGS    := -Wall -fPIC -Wno-unused-variable -O3 -g $(MODE) $(INCFLAGS)
FFLAGS      := -Wall -fPIC -Wno-unused-variable -Wno-intrinsic-shadow \
               -Wno-unused-label -fno-align-commons -O3 -g $(MODE) $(INCFLAGS)
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

.PHONY: all
.PHONY: clean distclean
.PHONY: test help env

###############################################################################
all: $(OBJS)

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
	@rm -f $(OBJDIR)/*
	@rm -f *~ *# */*~ */*#

distclean: clean

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

help: test

env: test
