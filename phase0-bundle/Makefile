# Phase 0 reproducible run.
IMAGE   ?= geneweaver-phase0:1.0
DB      ?= $(PWD)/../databases
ROOT    ?= $(PWD)/..
ENGINE  ?= podman

.PHONY: image fetch run shell local-fetch local-run

image:  ## build the container image
	$(ENGINE) build -t $(IMAGE) -f Containerfile .

fetch:  ## download + verify databases into $(DB)
	$(ENGINE) run --rm -v $(DB):/data:z -v $(ROOT):/work:z $(IMAGE) \
	  /work/phase0-bundle/scripts/fetch_databases.sh

run:    ## full Phase 0 inside the container
	$(ENGINE) run --rm -v $(DB):/data:z -v $(ROOT):/work:z $(IMAGE) \
	  /work/phase0-bundle/scripts/run_phase0.sh /work

shell:  ## interactive shell in the image
	$(ENGINE) run --rm -it -v $(DB):/data:z -v $(ROOT):/work:z $(IMAGE)

local-fetch:  ## same, without a container (needs diamond on PATH)
	GENEWEAVER_DB=$(DB) scripts/fetch_databases.sh

local-run:
	GENEWEAVER_DB=$(DB) scripts/run_phase0.sh $(ROOT)
