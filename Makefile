#### @file
#### @author IWB
#### @brief and a makefile

HOST=UW

F77 = gfortran
F90 = gfortran 
F95 = gfortran
FFLAGS = -lstdc++  -O2 -fopenmp  -fdefault-real-8
# those will be more like platform independent but we can put them here
#########################################################################
#
MAIN  = opt libmbpol.a
#
OBJ90=conj.o calc_grad.o
#
default: $(MAIN) $(TEST)
.PHONY:  clean
#
# MAKE THE HOST SPECIFIC OBJECT LIST 
OBJ77H = $(patsubst %.o, %_$(HOST).o,$(OBJ77))
OBJ90H = $(patsubst %.o, %_$(HOST).o,$(OBJ90))

OBJECTS = $(OBJ77H) $(OBJ90H) 


%_$(HOST).o: %.f95
	$(F95) -c $(FFLAGS) $(DEFS)  -o $@ $<

%_$(HOST).o: %.f
	$(F90) -c $(FFLAGS) $(DEFS)  -o $@ $<

%_$(HOST).o: %.F
	$(F77) -c $(FFLAGS) $(DEFS)  -o $@ $<

$(MAIN): $(OBJ77H) $(OBJ90H)  opt.f 
	 $(F95) opt.f libmbpol.a $(FFLAGS) $(OBJ90H) $(OBJ77H) -o $(MAIN) 


clean:
	rm *.o	
