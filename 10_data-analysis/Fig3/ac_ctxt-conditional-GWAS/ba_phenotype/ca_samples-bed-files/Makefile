SOURCE_DIR = $(abspath da_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.out, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.out : $(SOURCE_DIR)/%.target
	#OLD NAME ~/bin/descriptive_stat_format-bed.sh $< 2> $@
	~/bin/merged_chr-meth-per-sample.sh $< 2> $@

clean :
	- rm -rf *.out

