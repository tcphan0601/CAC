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
LIB =
SIZE =         size


EXE =          ../cac.2.4.0-a1
SRC =          $(wildcard *.cpp)
INC =          $(wildcard *.h)
OBJ =          $(SRC:.cpp=.o)

# Link target

$(EXE): $(OBJ)
	$(LINK) $(LINKFLAGS) $(OBJ) $(LIB) -o $(EXE)
	$(SIZE) $(EXE)

# Compilation rules

%.o:%.cpp
	$(CC) $(CCFLAGS) $(SHFLAGS) -c $<

#%.d:%.cpp
#	$(CC) $(CCFLAGS) $(DEPFLAGS) $(DEPFLAGS) $< > $@

# Individual dependencies

#DEPENDS = $(OBJ:.o=.d)
#sinclude $(DEPENDS)

clean:
	rm -rf *.o ../cac.2.4.0-a1 
