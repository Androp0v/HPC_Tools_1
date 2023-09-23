# Default Lapacke: Openblas at CESGA
#LDLIBS=-lopenblas

# Other systems (my Debian boxes, for example)
#LDLIBS=-llapacke

# Intel MKL at CESGA
# Module needed: imkl
# => module load openblas
# LDLIBS for intel compiler: icx (module needed: intel)
# Just invoke make like this: make CC=icx
#LDLIBS=-qmkl=sequential -lmkl_intel_lp64

# macOS + Homebrew
LDFLAGS=-L/opt/homebrew/opt/openblas/lib
CPPFLAGS=-I/opt/homebrew/opt/openblas/include
dgesv: dgesv.o timer.o main.o /opt/homebrew/opt/openblas/lib/libopenblas.a
	$(CC) dgesv.o timer.o main.o $(LDFLAGS) -lopenblas -o dgesv

# dgesv: dgesv.o timer.o main.o

clean:
	$(RM) dgesv *.o *~
