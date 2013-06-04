include ../../inc/common.mk

L += $(MYSQLLIBS) -lm -lz
MYLIBDIR = ../../lib/${MACHTYPE}
MYLIBS =  ${MYLIBDIR}/jkhgap.a ${MYLIBDIR}/jkweb.a

MYF = iteres

O = generic.o stat.o filter.o nearby.o density.o cpg.o cpgstat.o cpgfilter.o genomecov.o medip.o from_kent.o $(MYF).o

all: ${O} $(MYLIBS)
	${CC} ${COPT} -o $(MYF) $O ${MYLIBS} $L
	cp $(MYF) ~/bin/
clean:
	rm -f iteres *.o
