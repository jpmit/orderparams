# path to source files
SRCDIR = src
OBJDIR = src
CXX = g++
CXXFLAGS = -Wno-deprecated -O3
LDFLAGS = 
LDLIBS = -l gsl -l blas
OBJS = $(addprefix $(OBJDIR)/, main.o conncomponents.o orderparameter.o \
         opfunctions.o readwrite.o qlmfunctions.o gyration.o diagonalize.o \
         qdata.o)

all: orderparams

orderparams: $(OBJS)
	g++ $(LDFLAGS) -o test/testorderparams $(OBJS) $(LDLIBS)

main.o : main.cpp readwrite.h box.h orderparameter.h gyration.h qdata.h

conncomponents.o : conncomponents.cpp particle.h box.h

orderparameter.o : orderparameter.cpp orderparameter.h box.h particle.h \
                   constants.h opfunctions.h conncomponents.h qlmfunctions.h \
                   gyration.h utilityfunctions.h readwrite.h

opfunctions.o : opfunctions.cpp constants.h

readwrite.o : readwrite.cpp particle.h

qlmfunctions.o : qlmfunctions.cpp constants.h particle.h box.h opfunctions.h \
                 constants.h

gyration.o : gyration.cpp particle.h box.h conncomponents.h utilityfunctions.h

diagonalize.o : diagonalize.cpp

qdata.o : qdata.cpp qdata.h box.h particle.h qlmfunctions.h constants.h \
          conncomponents.h orderparameter.h utilityfunctions.h gyration.h

clean:
	rm -f $(OBJDIR)/*.o
