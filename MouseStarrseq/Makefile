SHELL=/bin/bash

REPOSITORIES=${HOME}/Software/repositories

ledpred_branch=develop
ledpred_sha=1dc0c268178cd0e59b24a78a5813a9e654547598
ledpred_path=${REPOSITORIES}/ledpred2

install_ledpred:
	cd ${ledpred_path}; git fetch origin; git checkout ${ledpred_branch}; git pull; git reset --hard $(ledpred_sha)
	cd ${ledpred_path}; Rscript -e 'try(remove.packages("LedPred", lib="'${R_LIBS_USER}'"), silent=T)'
	cd ${ledpred_path}; Rscript -e 'devtools::with_libpaths(new = "'${R_LIBS_USER}'", devtools::install("."))'


