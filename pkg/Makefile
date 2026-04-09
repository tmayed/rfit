clean:
	rm -rf renv
	rm -rf .Rprofile
	rm -rf renv.lock

install:
	R -e 'renv::init()'

.PHONY: poc

poc:
	RENV_PATH="$(PWD)/renv" Rscript poc/poc.R