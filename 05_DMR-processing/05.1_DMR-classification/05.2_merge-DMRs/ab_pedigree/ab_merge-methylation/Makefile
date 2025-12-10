SOURCE_DIR = $(abspath ba_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.uniteOut, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.uniteOut : $(SOURCE_DIR)/%.target
	~/bin/merge_methylation_CG-DMR_pedigree.sh $< 2> $@
	~/bin/merge_methylation_C-DMR_pedigree.sh $< 2> $@

clean :
	- rm -rf *.uniteOut
	- rm *tsv
	- rm *bed
