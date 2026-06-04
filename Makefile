R_SCRIPT ?= Rscript

.PHONY: clean run

run:
	$(R_SCRIPT) main.R

clean:
	rm -rf data/raw data/processed results
