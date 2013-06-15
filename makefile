KENT=/opt/kent/src
SAMTOOLS=/opt/samtools

CC=gcc
COPT= -O -g
CFLAGS= -Wall -Werror -Wformat -Wimplicit -Wreturn-type -Wuninitialized
DFLAGS= -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE
INCLUDES= -I$(KENT)/inc -I$(KENT)/hg/inc -I$(SAMTOOLS)
MYSQLLIB= -rdynamic -L/usr/lib64/mysql -lmysqlclient -lz -lcrypt -lnsl -lm -lssl -lcrypto
L += ${MYSQLLIB} -lm -lz
MYLIBDIR = $(KENT)/lib/${MACHTYPE}
MYLIBS =  ${MYLIBDIR}/jkhgap.a ${MYLIBDIR}/jkweb.a

MYF = iteres

O = generic.o stat.o filter.o nearby.o cpgstat.o cpgfilter.o from_kent.o $(MYF).o

%.o: %.c
	${CC} ${COPT} ${CFLAGS} ${DFLAGS} $(INCLUDES) -o $@ -c $<

all: ${O} $(MYLIBS)
	${CC} ${COPT} -o $(MYF) $O ${MYLIBS} $(SAMTOOLS)/libbam.a -pthread -lssl $L
	cp $(MYF) ~/bin/
clean:
	rm -f iteres *.o
