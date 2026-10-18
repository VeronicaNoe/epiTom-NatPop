SOURCE_DIR = $(abspath $(wildcard *_targets))
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.uniteOut, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.uniteOut : $(SOURCE_DIR)/%.target
	~/bin/intersect_lowPriority_CG-DMR_pedigree_ONT.sh $< 2> $@
	~/bin/intersect_lowPriority_C-DMR_pedigree_ONT.sh $< 2 > $@

clean :
	- rm -rf *.uniteOut

