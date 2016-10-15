FFLAGS=-fPIC -Wall -O2 -fmax-errors=5
LDFLAGS=-L . -lfandg
LIB=libfandg.a

SRC= stark.f90 \
     gates.f90 \
	 fg5.f90 \
	 fg5ex.f90 \
	 fg6.f90 \
	 fg6ex.f90 \
	 fg7.f90 \
	 fg7ex.f90 \
	 fg8.f90 \
	 fg8ex.f90 \
	 fg9.f90 \
	 fg9ex.f90 \
	 fg10.f90 \
	 fg10ex.f90 \
	 fg11.f90 \
	 fg11ex.f90 \
	 fg12.f90 \
	 fg12ex.f90 \
	 fg13.f90 \
	 fg13ex.f90 \
	 fg14.f90 \
	 fg14ex.f90 \
	 fg15.f90 \
	 fg15ex.f90 \
# 	 fg16.f90 \
	 fg16ex.f90 \
	 fg17.f90 \
	 fg17ex.f90 \
	 fg18.f90 \
	 fg18ex.f90 \
	 fg19.f90 \
	 fg19ex.f90 \
	 fg20.f90 \
	 fg20ex.f90

OBJS=$(SRC:.f90=.o)

.PHONY: all clean module

all: main module

module: stark.pyf ${LIB}
	f2py -c $^

main: main.o ${LIB}
	gfortran -g -o $@ $^ ${LDFLAGS}

${LIB}: ${OBJS}
	ar rcs $@ $^

%.o : %.f90
	gfortran -g -o $@ $(FFLAGS) -c $<

clean:
	rm -f $(OBJS) main.o main *.so *.a
