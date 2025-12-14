SOURCE_DIR = $(abspath $(wildcard *_targets))
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.mappingOUT, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.mappingOUT : $(SOURCE_DIR)/%.target
	~/bin/sort_RNAseq.sh $< 2> $@

clean :
	- rm -rf *.mappingOUT

