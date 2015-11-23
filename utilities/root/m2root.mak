# automatic settings for root

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

# Linux with egcs, gcc 2.9x, gcc 3.x (>= RedHat 5.2)
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC -I$(ROOTSYS)/include
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
#

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

# sources
SRCS		= m2root.C
OBJS		= m2root.o

# program
PROGRAM		= m2root

all:	$(PROGRAM)
	@echo "all done"

$(PROGRAM):	$(OBJS)
		@echo "Linking $(PROGRAM) ..."
		$(LD) $(LDFLAGS) $(OBJS) $(LIBS) $(GLIBS) -o $(PROGRAM)
		@echo "done"

clean:;		@rm -f $(OBJS) core
###
m2root.o:	m2root.h

