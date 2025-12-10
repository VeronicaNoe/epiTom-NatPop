SOURCE_DIR = $(abspath ba_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.out, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.out : $(SOURCE_DIR)/%.target
	~/bin/get_methylation_CG-DMRs_pedigree.sh $< 2> $@
	~/bin/get_methylation_C-DMRs_pedigree.sh $< 2> $@


clean :
	- rm -rf *.out


