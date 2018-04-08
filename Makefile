CXX=	clang++

# compiler flags.
CXXFLAGS=	-I/usr/local/include/eigen3
CXXFLAGS+=	-std=c++11
CXXFLAGS+=	-Ofast -g2 -mtune=native
LDFLAGS=	-lc++

CXXFLAGS+=	-DACC_LDOUBLE

CLEANFILES= *.o *.dSYM tools

all:	tools

clean:
	@rm -rf ${CLEANFILES}

tools.o:	tools.cc

