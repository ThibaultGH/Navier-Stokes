CC      = g++
FLAGS   = -G -Wall
LDFLAGS =
LIBRARY = -I/usr/local/include/eigen3 -I/usr/include/suitesparse -lumfpack -lcholmod -lccolamd -lcolamd -lcamd -lamd -lsuitesparseconfig -lblas
INCLUDE =

EXEC = main
SRCS = main.cpp function.cpp
OBJS = main.o function.o


all : $(EXEC)

$(EXEC) : $(OBJS)
		$(CC) -g $(LDFLAGS) -o $@ $^ $(LIBRARY)

%.o : %.cpp $(HEADERS)
		$(CC) -g $(CFLAGS) -c -o $@ $< $(LIBRARY)

.PHONY: clean, mrproper

clean :
		rm -rf *.o

mrproper: clean
		rm -rf $(EXEC)
