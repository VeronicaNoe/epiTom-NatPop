SOURCE_DIR = $(abspath ba_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.splitOut, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.splitOut : $(SOURCE_DIR)/%.target
	~/bin/split_by_context.sh $< 2> $@

clean :
	- rm -rf *.splitOut

