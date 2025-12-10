SOURCE_DIR = $(abspath ../ba_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.uniteOut, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.uniteOut : $(SOURCE_DIR)/%.target
	~/bin/fruit-leaf_intersect_highPriority.sh $< 2> $@

clean :
	- rm -rf *.uniteOut

