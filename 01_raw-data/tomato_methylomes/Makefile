SOURCE_DIR = $(abspath ../ab_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.merged, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.merged : $(SOURCE_DIR)/%.target
	~/bin/merge_fastqfiles.sh $< 2> $@

clean :
	- rm -rf *.merged
