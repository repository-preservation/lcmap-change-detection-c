TAG_PREFIX = usgseros
DOCKERHUB_ORG = $(TAG_PREFIX)

all: debian-c-ccdc ubuntu-c-ccdc

publish-docker: all debian-publish-c-ccdc ubuntu-publish-c-ccdc

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Common
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Note that the format of these base/common targets is:
# >-$ docker build -t <final image tag> <Dockerfile parent directory>

base-c-ccdc:
	@docker build -t $(TAG_PREFIX)/$(SYSTEM)-c-ccdc $(SYSTEM)

base-publish:
	@docker push $(DOCKERHUB_ORG)/$(REPO)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Debian
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

debian-c-ccdc:
	@SYSTEM=debian make base-c-ccdc

debian-publish-c-ccdc:
	@REPO=debian-c-ccdc make base-publish

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Ubuntu
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ubuntu-c-ccdc:
	@SYSTEM=ubuntu make base-c-ccdc

ubuntu-publish-c-ccdc:
	@REPO=ubuntu-c-ccdc make base-publish

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CentOS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# TBD

