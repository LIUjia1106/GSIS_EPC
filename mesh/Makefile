PROG = MSH

FCOMP = ifx
FLAG =  -O3 -g
# LIB_LIST = -qmkl #-L $(PREFIX) -lpardiso500-INTEL1301-X86-64
LIB_LIST = -qmkl  -L /opt/intel/oneapi/mkl/2024.2/lib/

OBJ = Mesh.o

$(PROG) : $(OBJ)
	$(FCOMP) $(OBJ) $(LIB_LIST) -o $(PROG)

$(PREFIX)%.o : %.f90
	$(FCOMP) $(FLAG) -c $< -o $@

#$(PREFIX)%.mod : %.mod
#	mv %.mod -f -t $(PREFIX)

clean :
	rm *.o $(PROG)


