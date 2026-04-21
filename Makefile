DATADIR := /app/data
DB      := $(DATADIR)/multiomics.db
PYTHON  := python

export PYTHONUNBUFFERED := 1

.PHONY: all

all: $(DATADIR)/abundance_missingness.csv

$(DATADIR)/.clinical.done:
	$(PYTHON) scripts/01_fetch_tcga_clinical.py $(DATADIR)
	touch $@

$(DATADIR)/.rnaseq.done: $(DATADIR)/.clinical.done
	$(PYTHON) scripts/02_fetch_tcga_rnaseq.py $(DATADIR)
	touch $@

$(DATADIR)/.rppa.done: $(DATADIR)/.clinical.done
	$(PYTHON) scripts/03_fetch_tcga_rppa.py $(DATADIR)
	touch $@

$(DATADIR)/.cptac.done:
	$(PYTHON) scripts/04_fetch_cptac_ccrcc.py $(DATADIR)
	touch $@

$(DB): $(DATADIR)/.rnaseq.done $(DATADIR)/.rppa.done $(DATADIR)/.cptac.done
	$(PYTHON) scripts/05_build_database.py $(DATADIR) $(DB)

$(DATADIR)/abundance_missingness.csv: $(DB)
	$(PYTHON) scripts/06_analyze.py $(DATADIR) $(DB)
