SOURCE_DIR = $(abspath ../ba_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.uniteOut, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.uniteOut : $(SOURCE_DIR)/%.target
	~/bin/merge_fruit-leaf-methylation_CG-DMR.sh $< 2> $@

clean :
	- rm -rf *.uniteOut

