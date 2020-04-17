OMPFLAG = -g -std=c++0x -w -fopenmp -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -s -O3
CFLAGS = -Wall -c ${OMPFLAG}
LFLAGS = -Wall ${OMPFLAG} -lm


PRJ	= $(shell pwd)/
LIB	= ${PRJ}lib/
BIN	= ${PRJ}bin/
CPPSRC = ${PRJ}src/cpp/
CSRC	= ${PRJ}src/c/
# C includes
INC_C1	= ${PRJ}include/c
INC_C		= -I${INC_C1}

# CPP includes
INC_CPP1 = ${PRJ}include/cpp
INC_CPP	= -I${INC_CPP1}

# INC_ROOT = -I$(ROOTSYS)/include/root
INC_ROOT = ""

all: c_lib cpp_lib HiSIF.o Final_Combine proc

dumper:
	g++ -Wno-write-strings ${INC_CPP} ${INC_C} src/cpp/dumper_test.cpp ${CSRC}r_read.c ${CSRC}r_write.c ${CSRC}funcErr.c -o dumpTest
new_all:
	#g++ ${LFLAGS} ${INC_ROOT} ${INC_C} ${INC_CPP} ${CSRC}*.c ${CPPSRC}*.cpp -o HiSIF `root-config --glibs`
	g++ ${LFLAGS} ${INC_C} ${INC_CPP} ${CSRC}*.c ${CPPSRC}*.cpp -o HiSIF

# combine all of the .o files, and run the root configuration with graphic elements
Final_Combine:
	#g++ ${INC_ROOT} ${INC_C} ${INC_CPP} -o ${BIN}HiSIF ${LIB}*.o `root-config --glibs`
	g++ ${INC_C} ${INC_CPP} -o ${BIN}HiSIF ${LIB}*.o

HiSIF: HiSIF.o
	g++ -L ${LFLAGS} ${INC_CPP} ${LIB}HiSIF.o -o HiSIF -lboost_filesystem -lboost_system -lboost_filesystem -lboost_system -static -static-libgcc

HiSIF.o: 
	# cd src/cpp; g++ ${INC_CPP} ${INC_ROOT} ${CFLAGS} HiSIF.cpp -o ${LIB}HiSIF.o
	cd src/cpp; g++ ${INC_CPP} ${CFLAGS} HiSIF.cpp -o ${LIB}HiSIF.o

#HiSIF.o: HiSIF.cpp Pair_Frag_2D_Peak_Info.hpp Frag_Info_With_Cutting_Site.hpp lib.hpp
#	g++ -I ~/lib/include lib ${CFLAGS} HiSIF.cpp -o HiSIF.o

cpp_lib:
	#cd src/cpp; g++ ${INC_ROOT} ${INC_CPP} ${CFLAGS} -c extractChrNum.cpp chrSort.cpp Get_Chr_Len.cpp Parser_Enzyme_Cutting_Site.cpp Map_Chr_Len.cpp PoissMix.cpp FindQT.cpp Get_Prox_Ligation_Frags.cpp HiSIF.cpp ; mv *.o ${LIB}
	cd src/cpp; g++ ${INC_CPP} ${CFLAGS} -c extractChrNum.cpp chrSort.cpp Get_Chr_Len.cpp Parser_Enzyme_Cutting_Site.cpp Map_Chr_Len.cpp PoissMix.cpp FindQT.cpp fragCopy.cpp Get_Prox_Ligation_Frags_File.cpp Pair_Frags_Interaction_Freq_File.cpp Gen_Ran_Dis_File.cpp HiSIF.cpp ; mv *.o ${LIB}

c_lib:
	cd src/c; g++ ${INC_C} ${CFLAGS} -c readline.c r_read.c r_write.c funcErr.c getchrlen.c addDirEnd.c getfilelist.c ;mv *.o ${LIB}

proc:
	#make c_lib
	cd src/c; gcc ${INC_C} -o proc readline.c r_write.c isdirectory.c getChrNum.c addDirEnd.c funcErr.c getfilelist.c frag_proc.c readInteractingRegions.c processesMain.c; mv proc ${PRJ}bin

clean:
	rm -rf lib/*.o HiSIF
	rm -rf bin*
