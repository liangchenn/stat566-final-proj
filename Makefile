R_SCRIPT ?= Rscript

.PHONY: all clean run

all: clean run

run:
	$(R_SCRIPT) main.R

clean:
	rm -rf data/raw data/processed results
