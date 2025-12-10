SOURCE_DIR = $(abspath ca_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.out, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.out : $(SOURCE_DIR)/%.target
	Rscript ~/bin/plot-manhattan-GWAS-metabolome.R $* 2> $@
clean :
	- rm -rf *.out

