SOURCE_DIR = $(abspath $(wildcard *_targets))
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.DMR.out, $(FILES))
SAMPLE_DIR = "/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/ah_methylome-comparison"
OUT_DIR    = $(CURDIR)
all : $(TARGETS)

$(CURDIR)/%.DMR.out : $(SOURCE_DIR)/%.target
#	~/bin/get-dmr-meth-sample_denovo-pedigree_ONT.sh $< 2> $@ # this is for sample dmr -> use sampleTargets
	~/bin/CG-DMR_denovo-ONT.sh $< $(SAMPLE_DIR) $(OUT_DIR) 2> $@ # this do both CG and C-DMR -> use pair-wise-DMRs.tar>

clean :
	- rm -rf *.out


