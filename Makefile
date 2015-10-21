all:
	make ccdc
	make classification

ccdc:
	cd ccdc && \
	make && \
	make install

classification:
	cd classification && \
	make && \
	make install

clean:
	cd ccdc && make clean
	cd classification && make clean

.PHONY: ccdc classification clean
