SOURCE_DIR = $(abspath $(wildcard *_targets))
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.out, $(FILES))
OUT_DIR    = $(CURDIR)

all : $(TARGETS)

$(CURDIR)/%.out : $(SOURCE_DIR)/%.target
	~/bin/deduplicate.sh $< $(OUT_DIR) 2> $@

clean :
	- rm -rf *.out

