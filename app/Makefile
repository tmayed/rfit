.PHONY: clean install danish danish_2m mixture_random test poc

clean:
	rm -rf renv .Rprofile renv.lock

install:
	R -e 'renv::init()'

test:
	# Run distribution tests by cd-ing into their directory first
	for f in tests/dist/test_*.R; do (cd $$(dirname $$f) && Rscript $$(basename $$f)); done
	# Run orchestrator tests
	(cd tests && Rscript test_fit_all.R)
	(cd tests && Rscript test_fit_all_2m.R)
	(cd tests && Rscript test_fit_all_3m.R)

poc_best_fit:
	(cd poc && Rscript poc_best_fit.R)

poc_best_fit_2m:
	(cd poc && Rscript poc_best_fit_2m.R)

poc_best_fit_2m_custom:
	(cd poc && Rscript poc_best_fit_2m_custom.R)

poc_best_fit_3m_custom:
	(cd poc && Rscript poc_best_fit_3m_custom.R)


poc_filter_dist:
	(cd poc && Rscript poc_filter_dist.R)

danish_fit:
	(cd poc && Rscript danish_fit.R)

danish_2m:
	(cd poc && Rscript danish_fit_2m.R)


fit_2m:
	(cd poc && Rscript fit_2m.R)

poc_fit_3m:
	(cd poc && Rscript poc_fit_3m.R)

2m_random_mixture_fit:
	(cd poc && Rscript 2m_random_mixture_fit.R)

fit_2m_w:
	(cd poc && Rscript fit_2m_w.R)

fit_2m_mean:
	(cd poc && Rscript fit_2m_mean.R)

danish_3m:
	(cd poc && Rscript danish_fit_3m.R)

3m_random_mixture_fit:
	(cd poc && Rscript 3m_random_mixture_fit.R)

poc: danish danish_2m 2m_random_mixture_fit danish_3m 3m_random_mixture_fit
