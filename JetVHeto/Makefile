PROG =	libjetvheto.a

ALLPROG = 	

# This will be used one day...
ALLPROGSRC =	

ALLPROGOBJS =	

SRCS =	qcd.f90 emsn_tools.f90 rad_tools.f90 resummation.f90 jetvheto_interface.f90 types.f90 warnings_and_errors.f90

POSTSRCS =	

OBJS =	qcd.o emsn_tools.o rad_tools.o resummation.o jetvheto_interface.o types.o warnings_and_errors.o

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
types.o:
qcd.o: types.o  
rad_tools.o: types.o qcd.o warnings_and_errors.o
emsn_tools.o: types.o  qcd.o rad_tools.o
resummation.o: types.o  qcd.o emsn_tools.o rad_tools.o warnings_and_errors.o 
jetvheto_interface.o: types.o resummation.o
warnings_and_errors.o:
