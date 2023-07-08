RUN=convert
RUNNAME=cl.11.5
NUM=100
CLONE_LIST=rdlist
convert=convert

#JLIB=/star/data01/pwg/heshu/libs/jam_64/src
#JLIB=/star/u/hexh/libs/jam-1.90596/src
JLIB=/star/data01/pwg/rseto/v3/jam1/jam-1.90596/src

.PHONY: makeseed cls
OTHERCLEAN = cls
CLONE_DEP = makeseed

./src/mysrc/libcjam.a: ./src/mysrc/cjam.f
	cd src/mysrc && make

_build/$(convert): $(convert).cc TJam.h ./src/mysrc/libcjam.a $(JLIB)/libjam.a
	@if ! [ -d _build ]; then \
		mkdir _build; \
	fi;
	g++ -O3 -c $(CFLAGS) $(convert).cc -o _build/$(convert).o
	g++ -O3 -o _build/$(convert) _build/$(convert).o $(LDFLAGS) -L./src/mysrc -L$(JLIB) -Wl,-Bstatic -lcjam -ljam -Wl,-Bdynamic -lgfortran $(LIBS)

makeseed:
	rm -f rdlist
	cd src && python extrd.py $(NUM) > ../rdlist

cls:
	cd ./src/mysrc && make clean
