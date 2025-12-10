SOURCE_DIR = $(abspath ../bc_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.uniteOut, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.uniteOut : $(SOURCE_DIR)/%.target
	~/bin/get_common-weighted-methylation_C-DMRs.sh $< 2> $@

clean :
	- rm -rf *.uniteOut

