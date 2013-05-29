# path to source files
SRC = src

main:
	g++ -Wno-deprecated -O3 -o orderparams ${SRC}/main.cpp ${SRC}/conncomponents.cpp ${SRC}/orderparameter.cpp ${SRC}/opfunctions.cpp ${SRC}/readwrite.cpp ${SRC}/qlmfunctions.cpp ${SRC}/gyration.cpp ${SRC}/diagonalize.cpp ${SRC}/qdata.cpp -l gsl -l blas
