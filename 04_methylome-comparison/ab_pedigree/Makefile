SOURCE_DIR = $(abspath ba_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.outDMR, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.outDMR : $(SOURCE_DIR)/%.target
	~/bin/toRun_pedigreeDMR.sh $< 2> $@

clean :
	- rm -rf *.outDMR

