SOURCE_DIR = $(abspath ../bb_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.uniteOut, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.uniteOut : $(SOURCE_DIR)/%.target
	~/bin/fruit_intersect_lowPriority_CG-DMR.sh $< 2> $@
	~/bin/fruit_intersect_lowPriority_C-DMR.sh $< 2 > $@

clean :
	- rm -rf *.uniteOut

