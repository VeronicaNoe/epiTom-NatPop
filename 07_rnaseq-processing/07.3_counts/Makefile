SOURCE_DIR = $(abspath aa_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.out, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.out : $(SOURCE_DIR)/%.target
	#~/bin/normalizing_RNAseq.sh $< 2> $@
	~/bin/counts_RNAseq.sh $< 2> $@

clean :
	- rm -rf *.out

