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

poc_fit_ln_gamma_dpnl:
	(cd poc && Rscript poc_fit_ln_gamma_dpnl.R)

poc_fit_ln_gamma:
	(cd poc && Rscript poc_fit_ln_gamma.R)

poc_filter_dist:
	(cd poc && Rscript poc_filter_dist.R)

poc_segment_data:
	(cd poc && Rscript poc_segment_data.R)
