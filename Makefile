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

docker:
	cd docker && make

dockerhub:
	cd docker && make publish-docker

.PHONY: ccdc classification clean docker
