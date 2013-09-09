# TODO: simplify this Makefile

# path to source files
SRCDIR = src
OBJDIR = src
CXX = g++
CXXFLAGS = -O3
LDFLAGS = 
LDLIBS = -l gsl -l blas
OBJS = $(addprefix $(OBJDIR)/, main.o conncomponents.o \
         opfunctions.o readwrite.o qlmfunctions.o gtensor.o diagonalize.o \
         qdata.o particlesystem.o orderparameters.o)
LDOBJS = $(addprefix $(OBJDIR)/, ldtool.o conncomponents.o opfunctions.o \
           readwrite.o qlmfunctions.o qdata.o particlesystem.o)

all: orderparams

orderparams: $(OBJS)
	g++ $(LDFLAGS) -o orderparams $(OBJS) $(LDLIBS)

ldtool: $(LDOBJS)
	g++ $(LDFLAGS) -o ldtool $(LDOBJS)

main.o : main.cpp particlesystem.h orderparameters.h qdata.h constants.h \
         utility.h gtensor.h

conncomponents.o : conncomponents.cpp typedefs.h particle.h box.h

opfunctions.o : opfunctions.cpp constants.h

readwrite.o : readwrite.cpp particle.h

qlmfunctions.o : qlmfunctions.cpp constants.h particle.h box.h opfunctions.h

gtensor.o : gtensor.cpp particlesystem.h particle.h box.h \
            conncomponents.h utility.h diagonalize.h gtensor.h

diagonalize.o : diagonalize.cpp

qdata.o : qdata.cpp qdata.h box.h particle.h qlmfunctions.h constants.h \
          conncomponents.h utility.h typedefs.h

particlesystem.o : particlesystem.cpp particlesystem.h readwrite.h box.h \
                   compile.h

orderparameters.o : orderparameters.cpp constants.h qlmfunctions.h \
                    qdata.h gtensor.h orderparameters.h

ldtool.o : ldtool.cpp

clean:
	rm -f $(OBJDIR)/*.o
