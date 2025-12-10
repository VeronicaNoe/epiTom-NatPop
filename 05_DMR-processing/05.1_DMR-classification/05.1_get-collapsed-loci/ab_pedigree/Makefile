SOURCE_DIR = $(abspath ba_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.loci, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.loci : $(SOURCE_DIR)/%.target
	~/bin/CG-DMR-loci_pedigree.sh $< 2> $@ # both ctxts

clean :
	- rm -rf *.loci
	- rm -rf *.bed

