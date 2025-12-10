SOURCE_DIR = $(abspath ba_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.mergedOut, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.mergedOut : $(SOURCE_DIR)/%.target
	~/bin/merge_annotation.sh $< 2> $@

clean :
	- rm -rf *.mergedOut

