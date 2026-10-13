# Define variables
PKG_DIR = .

# Phony targets don't correspond to actual files, they are just action names
.PHONY: all version-patch document check build install clean

all: build install

version-patch:
	Rscript -e "usethis::use_version('patch' )"

vignettes:
	Rscript -e "devtools::build_vignettes('$(PKG_DIR)')"

document:
	Rscript -e "devtools::document('$(PKG_DIR)')"

check: document
	Rscript -e "devtools::check('$(PKG_DIR)', vignettes = FALSE)"

build: check
	Rscript -e "pkgs <- c('BiocStyle','AnnotationHub','ensembldb','org.Hs.eg.db'); \
	  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]; \
	  if (length(missing) > 0) BiocManager::install(missing)"
	Rscript -e "devtools::build('$(PKG_DIR)', vignettes = FALSE)"

install: build
	Rscript -e "devtools::install('$(PKG_DIR)')"

clean:
	@echo "Cleaning up..."
	# Add commands here to remove generated files (e.g., from data-raw scripts)