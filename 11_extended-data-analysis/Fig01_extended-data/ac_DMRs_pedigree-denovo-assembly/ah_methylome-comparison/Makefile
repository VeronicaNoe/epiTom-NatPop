SOURCE_DIR = $(abspath $(wildcard *_targets))
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.outDMR, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.outDMR : $(SOURCE_DIR)/%.target
	~/bin/review_toRun_denovo-pedigreeDMR.sh  $< 2> $@

clean :
	- rm -rf *.outDMR

