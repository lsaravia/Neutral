# This file is automatically generated by RHIDE Version 1.4.7
# created from within RHIDE
vpath_src=.. ../../randlib/src ../../fortify ../../SpatialAnalysis/mfsba ../../SpatialAnalysis
vpath %.c $(vpath_src)
vpath %.cpp $(vpath_src)
vpath %.hpp $(vpath_src)
vpath %.h   $(vpath_src)

# The X11 base dir on your system
X11BASE=/usr/X11R6
# Add directories with X11 include files here
X11INCS=-I$(X11BASE)/include
# put X11 required libraries and directories here
X11LIBS=-L$(X11BASE)/lib -lX11

SDLDEFS = -D__XWIN__

I_DIRS=-I.. -I../../randlib/src -I../../fortify -I../../SpatialAnalysis/mfsba -I../../SpatialAnalysis

#P_DEFS=-DGRAPHICS -DEXP_DISP 
#P_DEFS=-DGRAPHICS -DUNIFORM_DISP
P_DEFS=-DGRAPHICS -DEXP_DISP -DHIERARCHICAL_CONT
#P_DEFS=-DGRAPHICS -DPOWER_DISP -DRANGE_CHECKING   
#P_DEFS=-DGRAPHICS -DPOWER_DISP -DHIERARCHICAL_CONT

CXXFLAGS = -O2 -Wall -std=gnu++0x $(I_DIRS) $(X11INCS)  $(SDLDEFS) $(P_DEFS)
#CXXFLAGS = -g -Wall -std=gnu++0x $(I_DIRS) $(X11INCS)  $(SDLDEFS) $(P_DEFS)


L = -lm -lgrx20S -lSDL $(X11LIBS)


O=bgi.o cabase.o IpsNeutral.o mfSBA.o RWFile.o\
	IpsNeutralAux.o IpsNeutralStats.o MainIPS.o 


MAIN_TARGET=ipsNeutral
all: $(O)
	g++ -o $(MAIN_TARGET) $(O) $(L)

clean:
	rm $(MAIN_TARGET) *.o 

	
all:

cabase.o: cabase.cpp cabase.hpp makefile

bgi.o: bgi.cpp makefile

IpsNeutral.o: IpsNeutral.cpp makefile IpsNeutral.h ran.hpp

IpsNeutralAux.o: IpsNeutralAux.cpp makefile IpsNeutral.h

MainIPS.o: MainIPS.cpp makefile IpsNeutral.h

IpsNeutralStats.o: IpsNeutralStats.cpp makefile IpsNeutral.h

mfSBA.o: mfSBA.cpp makefile 

RWFile.o: RWFile.cpp makefile 
