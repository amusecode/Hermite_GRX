# standard amuse configuration include
# config.mk will be made after ./configure has run
#AMUSE_DIR?=/disks/strw3/por/amuse-svn
#AMUSE_DIR?=/disks/strw3/por/amuse-svn
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
export AMUSE_DIR

-include $(AMUSE_DIR)/config.mk

MPICXX   ?= mpicxx

CFLAGS   += -Wall -g -std=c++11 -I$(AMUSE_DIR)/lib/stopcond -pthread
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS) 

OBJS = interface.o

CODELIB = src/libhermitepn.a

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: hermitepn_worker

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ hermitepn_worker worker_code.cc
	$(RM) -f *~
	make -C src clean

$(CODELIB):
	make -C src all

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py HermitePNInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H -i amuse.community.interface.stopping_conditions.StoppingConditionInterface interface.py HermitePNInterface -o $@

hermitepn_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $< -o $@ $(OBJS) $(CODELIB) -L./src -L$(AMUSE_DIR)/lib/stopcond -lstopcond -lhermitepn

interface.o: interface.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<
