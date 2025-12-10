SOURCE_DIR = $(abspath ba_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.out, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.out : $(SOURCE_DIR)/%.target
#	~/bin/toRun_format-mvcf-GWAS-metabolome.sh $< 2> $@
#	~/bin/get_mvcf-LD-prunning.sh $< 2> $@
	Rscript ~/bin/merge-DMRs-for-GWAS-metabolome.R $< 2> $@

clean :
	- rm -rf *.out

