hn = $(shell hostname)
username = $(shell whoami)

CPPSETTINGS+= -std=c++11
CPPFLAGS+= -O3 -Wno-all -Wno-deprecated -w -Wno-reorder #-march=native
# CPPFLAGS+= -g -O0 -Wno-all -Wno-deprecated -w -Wno-reorder #-march=native

COMMON_OBJS = \
	Vector3.o \
	Matrix3.o \
	SpeedFunctions.o \
	MathFunctions.o

VPATH := ../../source/
INC=-I$(VPATH)
.DEFAULT_GOAL := roofline

#Linking stage
roofline: $(COMMON_OBJS) roofline_benchmark.cpp
	$(CC) $(CPPSETTINGS) $(LIB) $(INC) $^ -o $@
	@echo done

#Linking stage
benchmark: $(COMMON_OBJS) roofline_benchmark.cpp
	$(CC) $(CPPSETTINGS) $(LIB) -DBENCHMARK_RESULTS $(INC) $^ -o $@
	@echo done

#Compiling stage
%.o: %.cpp
	$(CC) $(CPPSETTINGS) $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f *.o
	rm -f *.s
	rm -f roofline
