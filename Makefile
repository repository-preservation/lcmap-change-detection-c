ccdc:
	cd ccdc && \
	make && \
	make install

classification:
	cd classification && make

.PHONY: ccdc classification
