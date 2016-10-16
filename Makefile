FFLAGS=-fPIC -Wall -O2 -fmax-errors=5
LDFLAGS=-L . -lfandg
LIB=libfandg.a

SRC= stark.f90 \
     gates.f90 \
	 fg_a0.f90 \
	 fg_a1.f90 \
	 fg_a3_2.f90 \
	 fg_a2.f90 \
	 fg_a0ex.f90 \
	 fg_a1ex.f90 \
	 fg_a3_2ex.f90 \
	 fg_a2ex.f90

OBJS=$(SRC:.f90=.o)

.PHONY: all clean module

all: main module

module: stark.pyf ${LIB}
	f2py -c $^ -lgfortran

main: main.o ${LIB}
	gfortran -g -o $@ $^ ${LDFLAGS}

${LIB}: ${OBJS}
	ar rcs $@ $^

%.o : %.f90
	gfortran -g -o $@ $(FFLAGS) -c $<

clean:
	rm -f $(OBJS) main.o main *.so *.a
