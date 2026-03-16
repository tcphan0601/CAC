SHELL = /bin/sh

#-------------------------------------------------------------------
# compiler/linker settings
# specify flags and libraries needed for your compiler

CC =           mpicxx
CCFLAGS =      -g -O3
SHFLAGS =      -fPIC
#DEPFLAGS =     -M

LINK =         mpicxx
LINKFLAGS =    -g -O
LINKLIB =      -lpthread
LIB =          Tecplot/lib/libstdc++.so.6

TECIOLIB =     Tecplot/lib/libtecio.so
EXTRAINCLUDE = -I Tecplot/include
SIZE =         size


EXE =          ./cac.2.4.09
SRC =          $(wildcard *.cpp)
INC =          $(wildcard *.h)
OBJ =          $(SRC:%.cpp=obj/%.o)

# Link target

$(EXE):	$(OBJ)
	$(LINK) $(LINKFLAGS) $(OBJ) $(TECIOLIB) $(LIB) -o $(EXE)
	$(SIZE) $(EXE)

# Compilation rules

obj/%.o:%.cpp
	@ mkdir -p obj
	$(CC) $(CCFLAGS) ${EXTRAINCLUDE} $(SHFLAGS) -c $< -o $@

#%.d:%.cpp
#	$(CC) $(CCFLAGS) $(DEPFLAGS) $(DEPFLAGS) $< > $@

# Individual dependencies

#DEPENDS = $(OBJ:.o=.d)
#sinclude $(DEPENDS)

clean:
	rm -rf obj ./cac.2.4.09

clean-all:

	rm -rf obj ./cac.2.4.09
#	rm -f *.gch
	rm -f dump_tecplot.cpp dump_tecplot.h write_tecplot.cpp write_tecplot.h
	rm -f style_command.h style_dump.h
	cp Tecplot/style_dump_orgin.h ./style_dump.h
	cp Tecplot/style_command_orgin.h ./style_command.h

# Adding Tecplot output
tecplot:
	cp Tecplot/dump_tecplot.cpp ./
	cp Tecplot/dump_tecplot.h ./
	cp Tecplot/write_tecplot.cpp ./
	cp Tecplot/write_tecplot.h ./
	cp Tecplot/style_command.h ./
	cp Tecplot/style_dump.h ./
