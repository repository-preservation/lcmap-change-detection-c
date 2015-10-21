ccdc:
	cd src && \
	make && \
	make install

classification:
	cd classification && make

.PHONY: ccdc classification
