SOURCE_DIR = $(abspath bb_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.uniteOut, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.uniteOut : $(SOURCE_DIR)/%.target
	~/bin/intersect_allGenes.sh $< 2> $@
	~/bin/intersect_highPriority.sh $< 2> $@
	~/bin/intersect_exon-intron.sh $< 2> $@

clean :
	- rm -rf *.uniteOut

