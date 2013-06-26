# path to source files
SRCDIR = src
OBJDIR = src
CXX = g++
CXXFLAGS = -Wno-deprecated -g
LDFLAGS = 
LDLIBS = -l gsl -l blas
OBJS = $(addprefix $(OBJDIR)/, main.o conncomponents.o orderparameter.o \
         opfunctions.o readwrite.o qlmfunctions.o gyration.o diagonalize.o \
         qdata.o particlesystem.o orderparameters.o)

all: orderparams

orderparams: $(OBJS)
	g++ $(LDFLAGS) -o test/testorderparams $(OBJS) $(LDLIBS)

main.o : main.cpp particlesystem.h orderparameters.h qdata.h constants.h

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

particlesystem.o : particlesystem.cpp particlesystem.h readwrite.h box.h

orderparameters.o : orderparameters.cpp

clean:
	rm -f $(OBJDIR)/*.o
