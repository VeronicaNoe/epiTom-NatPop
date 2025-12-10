SOURCE_DIR = $(abspath ba_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.out, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.out : $(SOURCE_DIR)/%.target
	~/bin/pairwise-divergence.sh $< 2> $@
	~/bin/pairwise-divergence_DMRs.sg $< 2> $@

clean :
	- rm -rf *.out

