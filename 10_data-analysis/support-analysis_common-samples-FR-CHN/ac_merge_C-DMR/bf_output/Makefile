SOURCE_DIR = $(abspath ../be_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.uniteOut, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.uniteOut : $(SOURCE_DIR)/%.target
	~/bin/merge_common-methylation_C-DMR.sh $< 2> $@

clean :
	- rm -rf *.uniteOut

