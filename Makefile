#On Midway, first load modules:
#module load gcc/4.8 cuda/8.0 mkl/2017

exe     = amideImap
CXX     = g++

HOST=$(firstword $(subst -, ,$(shell hostname)))
ifeq ($(HOST),midway)
	SOFTDIR=/home/strong1/install
else
	SOFTDIR=/usr/local
endif
## $(info $$SOFTDIR is [${SOFTDIR}])

INC     = -I$(SOFTDIR)/include/xdrfile
FLAGS	= -O2 -march=native -Wall -std=c++11 -MMD -MP -fopenmp -DUSEOMP
DEBUGFLAGS = -g -std=c++11 -Wall -fopenmp -MMD -MP #-DUSEOMP #-DDEBUG -pg
LIBS    = -lxdrfile -lm -ldl -lgomp -lpthread
LIBDIRS = -L$(SOFTDIR)/lib

#MMD and MP generate .d files
OBJ_DIR = OBJ
SRCS  = main.cpp timer.cpp input.cpp groFile.cpp itpFile.cpp charges.cpp traj.cpp calcFreq.cpp compare.cpp compareEnergy.cpp compareDipole.cpp

OBJS := $(SRCS:%.cpp=$(OBJ_DIR)/%.o)
DEPS := $(OBJS:%.o=%.d)
OBJS_DEBUG := $(SRCS:%.cpp=$(OBJ_DIR)/%_d.o)
DEPS_DEBUG := $(OBJS:%.o=%_d.d)

exed := "$(exe)_debug"

cpu: $(exe)
debug: $(exed)
all: cpu debug

#link
$(exe): $(OBJS)
	$(CXX) $(FLAGS) -o $(exe) $^ $(LIBDIRS) $(LIBS)
$(exed): $(OBJS_DEBUG)
	$(CXX) $(DEBUGFLAGS) -o $(exed) $^ $(LIBDIRS) $(LIBS)

#compile
$(OBJ_DIR)/%.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(INC) $(FLAGS) -c -o $@ $<
$(OBJ_DIR)/%_d.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(INC) $(DEBUGFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	rm -rf "$(OBJ_DIR)"

-include $(DEPS)

#with help from https://spin.atomicobject.com/2016/08/26/makefile-c-projects/
