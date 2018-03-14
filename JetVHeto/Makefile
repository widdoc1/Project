PROG =	libjetvheto.a

ALLPROG = 	

# This will be used one day...
ALLPROGSRC =	

ALLPROGOBJS =	

SRCS =	consts_mod.f90 qcd_mod.f90 emsn_tools_mod.f90 rad_tools_mod.f90 resummation_mod.f90 jetvheto_interface.f90

POSTSRCS =	

OBJS =	consts_mod.o qcd_mod.o emsn_tools_mod.o rad_tools_mod.o resummation_mod.o jetvheto_interface.o

POSTOBJS =
POSTLIB =

LIBS =

CC = cc
CFLAGS = -O
FC = gfortran
FFLAGS = -fno-second-underscore -fno-automatic -ff2c -O0 -g 
F90 = gfortran
F90FLAGS = -fno-second-underscore -fno-automatic -ff2c -O0 -g 
LDFLAGS =  


all: libjetvheto.a

libjetvheto.a: $(OBJS)
	ar cru libjetvheto.a $(OBJS)
	ranlib libjetvheto.a

ALL: libjetvheto.a  $(ALLPROG)

libclean:
	rm -f  $(ALLPROGOBJS) $(OBJS) $(POSTOBJS) 

clean:
	rm -f  $(ALLPROGOBJS) $(OBJS) $(POSTOBJS) *.mod *.d *.o

realclean:
	rm -f libjetvheto.a $(ALLPROG) $(ALLPROGOBJS) $(OBJS) $(POSTOBJS) *.mod *.d *.o

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90 
	$(F90) $(F90FLAGS) -c $<

# a couple of lines to help with the combination of automatic
# generation of banner.f90 and the multi-directory build.
types_mod.o:
consts_mod.o: types_mod.o
qcd_mod.o: types_mod.o consts_mod.o 
rad_tools_mod.o: types_mod.o consts_mod.o qcd_mod.o
emsn_tools_mod.o: types_mod.o consts_mod.o qcd_mod.o rad_tools_mod.o
resummation_mod.o: types_mod.o consts_mod.o qcd_mod.o emsn_tools_mod.o rad_tools_mod.o
jetvheto_interface.o: types_mod.o consts_mod.o resummation_mod.o
