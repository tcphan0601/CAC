SHELL = /bin/sh

#-------------------------------------------------------------------
# Compiler/linker settings
#-------------------------------------------------------------------

CC =           mpicxx
CCFLAGS =      -g -O3
SHFLAGS =      -fPIC
LINK =         mpicxx
LINKFLAGS =    -g -O
LIB =
SIZE =         size

#-------------------------------------------------------------------
# DEEPMD Package Settings
#-------------------------------------------------------------------

DEEPMD_PATH =
DEEPMD_LIB =
DEEPMD_INC =

ifneq ($(wildcard .deepmd_enabled),)
  CCFLAGS += -DDEEPMD_ENABLED
  DEEPMD_PATH = -LDEEPMD/lib
  DEEPMD_LIB = -Wl,--no-as-needed -ldeepmd_cc -Wl,-rpath=DEEPMD/lib
  DEEPMD_INC = -IDEEPMD/include/
endif


#-------------------------------------------------------------------
# TECPLOT Package Settings
#-------------------------------------------------------------------

TECIOLIB =
TECPLOT_INC =

ifneq ($(wildcard .tecio_enabled),)
  CCFLAGS += -DTECIO_ENABLED
  TECPLOT_LIB = TECPLOT/lib/libtecio.so
  TECPLOT_INC = -I TECPLOT/include
endif

#-------------------------------------------------------------------
# Build Targets
#-------------------------------------------------------------------

EXE =          ./cac.3.2.0
SRC =          $(wildcard *.cpp)
INC =          $(wildcard *.h)
OBJ =          $(SRC:%.cpp=obj/%.o)

all: $(EXE)

# Link target
$(EXE):	$(OBJ)
	$(LINK) $(LINKFLAGS) $(OBJ) $(TECPLOT_LIB) $(LIB) $(DEEPMD_PATH) $(DEEPMD_LIB) -o $(EXE)
	$(SIZE) $(EXE)

# Compilation rules
obj/%.o:%.cpp
	@ mkdir -p obj
	$(CC) $(CCFLAGS) ${TECPLOT_INC} $(DEEPMD_INC) $(SHFLAGS) -c $< -o $@

#-------------------------------------------------------------------
# Package Management Targets (LAMMPS Style)
#-------------------------------------------------------------------

# Installs the TECPLOT package by creating symbolic links to the
# necessary source files.
yes-tecplot:
	@echo "--- Installing package TECPLOT ---"
	@touch .tecio_enabled
	@rm -f obj/dump_custom.o obj/write_tecplot.o
	@echo "Package TECPLOT installed. Please re-build your executable."

# Uninstalls the TECPLOT package by removing the linked files and
# deleting the #include directives from the style header files.
no-tecplot:
	@echo "--- Uninstalling package TECPLOT ---"
	@rm -f .tecio_enabled obj/dump_custom.o obj/write_tecplot.o
	@echo "Package tecplot uninstalled."

# Installs the DEEPMD package
yes-deepmd:
	@echo "--- Installing package DEEPMD ---"
	@ln -sf DEEPMD/pair_base.cpp .
	@ln -sf DEEPMD/pair_base.h .
	@ln -sf DEEPMD/pair_deepmd.cpp .
	@ln -sf DEEPMD/pair_deepmd.h .
	@ln -sf DEEPMD/deepmd_version.h .
	@grep -q -F '#include "pair_deepmd.h"' style_pair.h || echo '#include "pair_deepmd.h"' >> style_pair.h
	@grep -q -F '#include "pair_base.h"' style_pair.h || echo '#include "pair_base.h"' >> style_pair.h
	@touch .deepmd_enabled
	@echo "Package DEEPMD installed. Please re-build your executable."

# Uninstalls the DEEPMD package
no-deepmd:
	@echo "--- Uninstalling package DEEPMD ---"
	@rm -f pair_base.cpp pair_base.h pair_deepmd.cpp pair_deepmd.h deepmd_version.h
	@sed -i '/#include "pair_deepmd.h"/d' style_pair.h
	@sed -i '/#include "pair_base.h"/d' style_pair.h
	@rm -f .deepmd_enabled
	@echo "Package DEEPMD uninstalled."

#-------------------------------------------------------------------
# Cleanup Targets
#-------------------------------------------------------------------

clean:
	@rm -rf obj ./cac.3.2.0

# A full clean uninstalls all packages and removes all build artifacts.
clean-all:
	@$(MAKE) no-tecplot
	@$(MAKE) no-deepmd
	@$(MAKE) clean

# Specify targets that are not files
.PHONY: all yes-tecplot no-tecplot yes-deepmd no-deepmd clean clean-all
