build:
	docker build -t crusty .

rebuild-c:
	sed "s/DEVELOP_MODE_FLAG/`date`/g" Dockerfile.in > Dockerfile
	make 

shell:
	docker run -it --entrypoint /bin/bash crusty  -s

run:
	docker run crusty
