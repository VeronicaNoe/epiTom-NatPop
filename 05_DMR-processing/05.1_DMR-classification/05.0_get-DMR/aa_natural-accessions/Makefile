SOURCE_DIR = $(abspath aa_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.DMR.out, $(FILES))
SAMPLE_DIR = "/mnt/disk2/vibanez/04_methylome-comparison/aa_natural-accessions"
OUT_DIR    = $(CURDIR)

all : $(TARGETS)

$(CURDIR)/%.DMR.out : $(SOURCE_DIR)/%.target
	~/bin/CG-DMR.sh $< $(SAMPLE_DIR) $(OUT_DIR) 2> $@
	~/bin/C-DMR.sh $< $(SAMPLE_DIR) $(OUT_DIR) 2> $@

clean :
	- rm -rf *.DMR.out

