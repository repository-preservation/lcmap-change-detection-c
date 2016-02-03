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

ubuntu-bash:
	docker run -i -t --entrypoint=/bin/bash usgseros/ubuntu-c-ccdc -s

debian-bash:
	docker run -i -t --entrypoint=/bin/bash usgseros/debian-c-ccdc -s

.PHONY: ccdc classification clean docker
