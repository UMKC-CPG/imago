SHELL := /bin/bash

all:

install:
	. /etc/environ.sh && use -e -r cmake-3.6.3 && . ../.imago/imagorc $$(pwd)/.. && \
	cd ../build/release && \
	export FC=h5fc && \
	cmake ../.. -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX="../.." && \
	make install && \
	unpackIMAGODB

clean:
	rm -rf ../build/release/*

distclean: clean
	rm -rf ../bin
	rm -rf ../share

.PHONY: all install clean distclean
