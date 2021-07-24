CXXFLAGS=-DNDEBUG -std=c++11 -O3 -w -fPIC
CC=g++
OBJ=src/ParseCommandLine.o src/common_utils.o RNA_class/RNA.o RNA_class/thermodynamics.o src/DynProgArray.o src/extended_double.o src/forceclass.o src/outputconstraints.o src/pfunction.o src/rna_library.o src/structure.o src/substructure.o src/TProgressDialog.o
OBJ2=CSSR/cssr_struct.o CSSR/PDBParser.o 
LDFLAGS=-lstdc++

exe/CSSR: ${OBJ} CSSR/CSSR.o ${OBJ2}
	${CC} -o exe/CSSR ${CXXFLAGS} CSSR/CSSR.o ${OBJ} ${OBJ2} ${LDFLAGS}
exe/ProbablePairRR: ${OBJ} ProbablePairRR/ProbablePairRR.o
	${CC} -o exe/ProbablePairRR ${CXXFLAGS} ${OBJ} ProbablePairRR/ProbablePairRR.o ${LDFLAGS}

CSSR/CSSR.o: CSSR/CSSR.cpp
	${CC} -c -o $@ ${CXXFLAGS} CSSR/CSSR.cpp
CSSR/PDBParser.o: CSSR/PDBParser.cpp
	${CC} -c -o $@ ${CXXFLAGS} CSSR/PDBParser.cpp
CSSR/cssr_struct.o: CSSR/cssr_struct.cpp
	${CC} -c -o $@ ${CXXFLAGS} CSSR/cssr_struct.cpp
ProbablePairRR/ProbablePairRR.o: ProbablePairRR/ProbablePairRR.cpp
	${CC} -c -o $@ ${CXXFLAGS} ProbablePairRR/ProbablePairRR.cpp
src/ParseCommandLine.o: src/ParseCommandLine.cpp
	${CC} -c -o $@ ${CXXFLAGS} src/ParseCommandLine.cpp
src/common_utils.o: src/common_utils.cpp
	${CC} -c -o $@ ${CXXFLAGS} src/common_utils.cpp
RNA_class/RNA.o: RNA_class/RNA.cpp
	${CC} -c -o $@ ${CXXFLAGS} RNA_class/RNA.cpp
RNA_class/thermodynamics.o: RNA_class/thermodynamics.cpp
	${CC} -c -o $@ ${CXXFLAGS} RNA_class/thermodynamics.cpp
src/DynProgArray.o: src/DynProgArray.cpp
	${CC} -c -o $@ ${CXXFLAGS} src/DynProgArray.cpp
src/extended_double.o: src/extended_double.cpp
	${CC} -c -o $@ ${CXXFLAGS} src/extended_double.cpp
src/forceclass.o: src/forceclass.cpp
	${CC} -c -o $@ ${CXXFLAGS} src/forceclass.cpp
src/outputconstraints.o: src/outputconstraints.cpp
	${CC} -c -o $@ ${CXXFLAGS} src/outputconstraints.cpp
src/pfunction.o: src/pfunction.cpp
	${CC} -c -o $@ ${CXXFLAGS} src/pfunction.cpp
src/rna_library.o: src/rna_library.cpp
	${CC} -c -o $@ ${CXXFLAGS} src/rna_library.cpp
src/structure.o: src/structure.cpp
	${CC} -c -o $@ ${CXXFLAGS} src/structure.cpp
src/substructure.o: src/substructure.cpp
	${CC} -c -o $@ ${CXXFLAGS} src/substructure.cpp
src/TProgressDialog.o: src/TProgressDialog.cpp
	${CC} -c -o $@ ${CXXFLAGS} src/TProgressDialog.cpp
all: exe/CSSR exe/ProbablePairRR
clean:
	rm ${OBJ} exe/CSSR exe/ProbablePairRR
