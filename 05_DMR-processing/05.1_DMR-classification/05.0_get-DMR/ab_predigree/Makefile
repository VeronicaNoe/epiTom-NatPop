SOURCE_DIR = $(abspath ba_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.DMR.out, $(FILES))
SAMPLE_DIR = "/mnt/disk2/vibanez/04_methylome-comparison/ab_pedigree"
OUT_DIR    = $(CURDIR)
all : $(TARGETS)

$(CURDIR)/%.DMR.out : $(SOURCE_DIR)/%.target
	#~/bin/get-dmr-meth-sample_pedigree.sh $< 2> $@ # this is for sample dmr
	~/bin/CG-DMR.sh $< $(SAMPLE_DIR) $(OUT_DIR) 2> $@ # this do both CG and C-DMR

clean :
	- rm -rf *.out

