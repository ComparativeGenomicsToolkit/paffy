rootPath = .
include ${rootPath}/include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

stPafDependencies = ${LIBDEPENDS}
stPafLibs = ${LDLIBS}

all: all_libs all_progs
all_libs: ${LIBDIR}/stPaf.a
all_progs: all_libs ${BINDIR}/stPafTests ${BINDIR}/paffy

sonLib:
	mkdir -p ${LIBDIR} ${BINDIR}
	cd submodules/sonLib && PKG_CONFIG_PATH=${CWD}/lib/pkgconfig:${PKG_CONFIG_PATH} ${MAKE}
	mkdir -p ${BINDIR} ${LIBDIR} ${INCLDIR}
	rm -rf submodules/sonLib/bin/*.dSYM
	ln -f submodules/sonLib/lib/sonLib.a ${LIBDIR}/sonLib.a
	ln -f submodules/sonLib/lib/cuTest.a ${LIBDIR}/cuTest.a
	ln -f submodules/sonLib/lib/sonLib.a ${LIBDIR}/libsonLib.a

stPafDependencies = ${sonLibDir}/sonLib.a ${sonLibDir}/cuTest.a

${sonLibDir}/sonLib.a : sonLib

${sonLibDir}/cuTest.a : sonLib

${LIBDIR}/stPaf.a : ${libSources} ${libHeaders}  ${stPafDependencies}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -c ${libSources}
	${AR} rc stPaf.a *.o
	${RANLIB} stPaf.a
	mv stPaf.a ${LIBDIR}/

${BINDIR}/paffy : paffy_main.c ${LIBDEPENDS} ${stPafDependencies} ${commonPafLibs} ${libHeaders}
	${CC} ${CPPFLAGS} ${CFLAGS} -o ${BINDIR}/paffy paffy_main.c ${libSources} ${commonPafLibs} ${LDLIBS}

${BINDIR}/stPafTests : ${libTests} ${LIBDIR}/stPaf.a ${stPafDependencies}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/stPafTests ${libTests} ${libSources} ${LIBDIR}/stPaf.a ${stCafLibs} ${LDLIBS}

clean :
	cd submodules/sonLib && ${MAKE} clean
	rm -rf *.o ${BINDIR} ${LIBDIR}

test : all
	${BINDIR}/stPafTests
