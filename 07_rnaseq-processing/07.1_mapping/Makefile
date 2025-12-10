SOURCE_DIR = $(abspath aa_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.out, $(FILES))
OUT_DIR    = $(CURDIR)

all : $(TARGETS)

$(CURDIR)/%.out : $(SOURCE_DIR)/%.target
	~/bin/mapping_RNAseq_STAR.sh $< $(OUT_DIR) 2> $@

clean :
	- rm -rf *.out

