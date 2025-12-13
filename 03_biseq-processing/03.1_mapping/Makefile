SOURCE_DIR = $(abspath $(wildcard *_targets))
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.mappingOut, $(FILES))
OUT_DIR    = $(CURDIR)
REF_TYPE   = "pseudo-genome"

all : $(TARGETS)

$(CURDIR)/%.mappingOut : $(SOURCE_DIR)/%.target
	~/bin/mapping_WGBS_bismark.sh $< $(REF_TYPE) $(OUT_DIR) 2> $@

clean :
	- rm -rf *.mappingOut

