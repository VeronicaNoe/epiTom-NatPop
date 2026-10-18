SOURCE_DIR = $(abspath $(wildcard *_targets))
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.splitOut, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.splitOut : $(SOURCE_DIR)/%.target
	~/bin/review_filter_meth.sh $< 2> $@
	~/bin/review_split_by_chromosomes.sh $< 2> $@

clean :
	- rm -rf *.splitOut

