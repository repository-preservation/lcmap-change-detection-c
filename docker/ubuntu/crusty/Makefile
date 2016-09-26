build: Dockerfile
	docker build -t crusty .

Dockerfile:
	sed "s/DEVELOP_MODE_FLAG//g" Dockerfile.in > Dockerfile

cache-bust-Dockerfile:
	sed "s/DEVELOP_MODE_FLAG/`date`/g" Dockerfile.in > Dockerfile

rebuild-c: cache-bust-Dockerfile
	make 

shell:
	docker run -it --entrypoint /bin/bash crusty  -s

run:
	docker run crusty
