SOURCE_DIR = $(abspath aa_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.out, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.out : $(SOURCE_DIR)/%.target
#	~/bin/review_liftoff-DMRs-SL2.5-M82.sh $< 2> $@ #done in folder ../04
	~/bin/review_dmr2m82-mapq-filter.sh $< 2> $@ #outputs in ../04

clean :
	- rm -rf *.out

