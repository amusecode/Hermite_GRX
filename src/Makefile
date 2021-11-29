CPP_FILES := $(wildcard *.cpp)
OBJ_FILES := $(CPP_FILES:.cpp=.o)
CC_FLAGS := -std=c++11 -Wall -O4 -I$(AMUSE_DIR)/lib/stopcond
LD_FLAGS := -pg

AR = ar ruv
RANLIB = ranlib

all: libhermite_grx.a

libhermite_grx.a: $(OBJ_FILES)
	rm -f $@
	$(AR) $@ $(OBJ_FILES)
	$(RANLIB) $@

%.o: %.cpp
	g++ $(CC_FLAGS) -c -o $@ $<

clean:
	rm -f *.o *.a *~
