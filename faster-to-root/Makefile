
#
# WARNING : g++ and root-system are requiered to compile this example
#


#
# If pkg-config isn't installed on your system, comment the following lines and
# set the fasterac flags as indicated with your own paths:
#
# FASTERAC_CFLAGS = -I/usr/include
# FASTERAC_LIBS   = -L/usr/lib -lfasterac
#
FASTERAC_CFLAGS = $(shell pkg-config  --cflags libfasterac)
FASTERAC_LIBS   = $(shell pkg-config  --libs   libfasterac)
#

ROOT_CFLAGS     = $(shell root-config --cflags)
ROOT_LIBS       = $(shell root-config --libs) -ltbb     # TODO tbb


CC        = g++
CFLAGS    = ${FASTERAC_CFLAGS} ${ROOT_CFLAGS}
LIBS      = ${FASTERAC_LIBS}   ${ROOT_LIBS}

SRCEXE    = $(shell ls *.C)
EXE       = $(SRCEXE:.C=)

all : $(EXE)

$(EXE): $(SRCEXE)
	${CC} $@.C -o $@ ${CFLAGS} ${LIBS}

rootlogon :
	@echo "---"

clean :
	rm -f *.o
	rm -f *faster2tree*.root
	rm -f $(EXE)


