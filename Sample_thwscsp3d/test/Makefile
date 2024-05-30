
include ../make.inc

F90 = gfortran

SRC= laplace_3d_test5.f
EXES=$(subst .f,.exe, $(SRC))

all : $(EXES)

$(EXES) : $(SRC)

%.exe : %.f
	rm -f $@
	$(F90) $< -o $@ -L../lib -l fishpack
	./$@

clean :
	rm -f $(EXES)
