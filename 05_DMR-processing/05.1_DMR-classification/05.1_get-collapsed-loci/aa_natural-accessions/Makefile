SOURCE_DIR = $(abspath $(wildcard *_targets))
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.loci, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.loci : $(SOURCE_DIR)/%.target
	~/bin/CG-DMR-loci.sh $< 2> $@
	~/bin/C-DMR-loci.sh $< 2> $@

clean :
	- rm -rf *.loci

