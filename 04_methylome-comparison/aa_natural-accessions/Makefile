SOURCE_DIR = $(abspath aa_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.outDMR, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.outDMR : $(SOURCE_DIR)/%.target
	~/bin/toRun_accDMR.sh $< 2> $@

clean :
	- rm -rf *.outDMR

