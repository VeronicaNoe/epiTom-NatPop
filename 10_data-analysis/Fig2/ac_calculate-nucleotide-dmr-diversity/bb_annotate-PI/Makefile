SOURCE_DIR = $(abspath ca_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.out, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.out : $(SOURCE_DIR)/%.target
	~/bin/annotate-DMRs_epiPI.sh $< 2> $@
	~/bin/annotate-SNPs_forPI.sh $< 2> $@

clean :
	- rm -rf *.out

