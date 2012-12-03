# By Jixie Zhang
# Makefile for fortran program, you can specify your executable by command line: 
# "make PROG=your_exe"
# Will look for program file itself under current dir, make sure remove
# all unrelavant source file out of current dir
#
MYOS := $(shell uname -s)
ARCH := $(shell arch)
KEY  := $(ARCH)
ifeq ($(ARCH),x86_64)
  ifeq ($(COMPILER_BIT),32)
    KEY := $(ARCH)_m32
  endif
endif
#####################################################
PROG   :=HRSTransport
SRCDIR := .
INCDIR := .
#OBJDIR := obj.$(KEY)
OBJDIR := obj.$(ARCH)
MYOS   := $(shell uname)
#####################################################
CC = g++
FF = f77
LD = f77
FFLAGS = -I$(INCDIR) -O2
CFLAGS = -I$(INCDIR) -O2
AR = ar r

#for 64-bit mac
ifeq ($(MYOS),Darwin)
    FF = gfortran
    LD = gfortran
    FFLAGS += -m64
    CFLAGS += -m64
    AR = libtool -static -o
endif

ifeq ($(ARCH),x86_64)
  ifeq ($(COMPILER_BIT),32)
    FFLAGS += -m32
    CFLAGS += -m32
  else
    FF = gfortran
    LD = gfortran
    FFLAGS += -m64
    CFLAGS += -m64
  endif
endif

#fortran options are crutial, some of them should not be included here
#depends on what options have been used in the cernlib 
ifeq ($(FF),gfortran)
  ifeq ($(MYOS),Darwin)
    FFLAGS  += -fno-second-underscore \
 -fno-automatic -fbounds-check -fno-range-check\
 -funroll-all-loops -fdollar-ok -ffixed-line-length-none \
 -Wunused -Wuninitialized 
# -finit-local-zero #this one is not recgonized now
#-Wno-globals , this one is not recgoganized
#-fno-leading-underscore 
#in Darwin, do not use -fno-leading-underscore
  else
     FFLAGS += -fno-leading-underscore -fno-second-underscore \
 -fno-automatic -fbounds-check \
 -funroll-all-loops -fdollar-ok -ffixed-line-length-none \
 -Wunused -Wuninitialized \
 -fno-range-check
  endif 
else
  FFLAGS  += -fno-leading-underscore -fno-second-underscore \
 -fno-automatic -fbounds-check \
 -funroll-all-loops -fdollar-ok -ffixed-line-length-none \
 -Wunused -Wuninitialized \
 -finit-local-zero -Wno-globals 
endif

ifdef DEBUG
 FFLAGS += -g	
 CFLAGS += -g	
endif


################################ Libraries
ifndef SYSLIBS
 SYSLIBS  = -lnsl -lX11 -lstdc++
endif

#for SunOS, cernlib 2001 need libsocket.a
ifeq ($(MYOS),SunOS)
 SYSLIBS  += -lsocket
endif

ifndef CERNLIBS
 #this program used 2001 cernlib /apps/cernlib/sun4_solaris8/2001 in jlabs2,
 ifndef CERN_ROOT
  CERN_ROOT := /apps/cernlib/i386_rhel3/2003
 endif
 CERNLIBS = -L$(CERN_ROOT)/lib -lpdflib804 -lmathlib -lphtools \
 -lgeant321 -lpawlib -lgraflib -lgrafX11 -lpacklib -lkernlib 
endif

#if you are a clas member ......
ifdef USE_CLAS_LIB
 ifndef CLASLIBS
  CLASLIBS = -L$(CLAS_LIB) -lbankdefs -lc_bos_io -lfputil  -lrecutl \
  -lbos -lfpack -lcc -lmapmanager -lclasutil
 endif
else
 CLASLIBS =  # use this # to set a space  
endif

LIBS = $(SYSLIBS) $(CLASLIBS) $(CERNLIBS)
#===========================

###################################
ifndef PROG
#  NFILES = $(shell ls -1 *.[fF] | wc -l | sed s/\ //g)   
#  ifeq ($(NFILES),1)
#    PROG = $(shell grep -i program *.[fF] | awk '$$1=="program"&&$$3==""{print $$2}' | sed s/\ //g) 
#  else
#    PROG = $(shell grep -i program *.[fF] | awk '$$2=="program"&&$$4==""{print $$3}' | sed s/\ //g)
#  endif
  PROG = $(shell grep -i program *.[fF] | awk '$$2=="program"&&$$4==""{print $$3}' | sed s/\ //g)
endif
ifeq ($(PROG),) 
  PROG = $(shell grep -i program *.[fF] | awk '$$1=="program"&&$$3==""{print $$2}' | sed s/\ //g)  
endif

PROG := $(shell echo ${PROG} | sed s/\ //g)
EXE   = $(PROG)_$(KEY)
ifdef DEBUG
EXE  :=$(EXE)_debug
endif

SOURCES  := $(wildcard $(SRCDIR)/*.[fF])
SOURCES  += $(wildcard $(SRCDIR)/*.cc)
SOURCES  += $(wildcard $(SRCDIR)/*.cpp)
#OBJS     = $(shell echo $(SOURCES) | sed s/f\ /o\ /g)
OBJS    := $(addsuffix .o, $(basename $(SOURCES)))
OBJS    := $(patsubst $(SRCDIR)/%.o,$(OBJDIR)/%.o,$(OBJS))

###################################
$(OBJDIR)/%.o:%.f
	$(FF) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o:%.F
	$(FF) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o:%.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)/%.o:%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

###################################
.PHONY: lib exe clean

all:  $(OBJDIR) lib #$(EXE)

$(OBJDIR):
	@if [ ! -d $(OBJDIR) ] ; then mkdir -p $(OBJDIR) ;fi

$(EXE): $(OBJDIR) $(OBJS) 
	$(FF) -o $(EXE) $(OBJS) $(LIBS) $(FFLAGS)
	@echo "Linking $(EXE) ... done!"

#this part is to support CLAS_PACK make mechanism
exe: lib
	$(FF) -o $(EXE) $(OBJS) -L$(OBJDIR) -l$(PROG) $(LIBS) $(FFLAGS)
	@echo "Linking $(EXE) ... done!"

#this part is to support CLAS_PACK make mechanism
lib: $(OBJDIR) $(OBJS)
#	@rm -f $(OBJDIR)/lib$(PROG).a
	@$(AR) $(OBJDIR)/lib$(PROG).a $(OBJS)
#	@if [ -e $(OBJDIR)/lib$(PROG).a ] ; then echo "Finish creating target lib"  $(OBJDIR)/lib$(PROG).a ;fi
#	@libtool -static  -o $(OBJDIR)/lib$(PROG).a $(OBJS)
#	@ar r $(OBJDIR)/lib$(PROG).a $(OBJS)

clean:
	rm -f $(OBJS) $(OBJDIR)/lib${PROG}.a core *~ *#

delete:
	rm -f $(EXE)

help:
	@echo PROG = ${PROG}
	@echo SOURCES = ${SOURCES}
	@echo OBJS = ${OBJS}
	@echo EXE = ${EXE}
	@echo LIBS = ${LIBS}
	@echo FFLAGS = ${FFLAGS}
	@echo CLASLIBS = ${CLASLIBS}
