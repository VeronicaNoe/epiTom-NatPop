
SOURCE_DIR = $(abspath da_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.out, $(FILES))
SIG     = /mnt/disk2/vibanez/10_data-analysis/Fig3/aa_GWAS-DMRs/sig
NONSIG  = /mnt/disk2/vibanez/10_data-analysis/Fig3/aa_GWAS-DMRs/nonSig-GIF
RESULTDIR = /mnt/disk2/vibanez/10_data-analysis/Fig3/ab_data-analysis/results/

all : $(TARGETS)

$(CURDIR)/%.out : $(SOURCE_DIR)/%.target
	~/bin/toRun_merge-GIF.sh $< 2> $@
	cat $(RESULTDIR)/00.0_list2move-from-sig2nonsig.tsv | parallel -j10 'ls -lh $(SIG)/{}.* && mv $(SIG)/{}.* $(NONSIG)/'
	#while read line; do \
	#	mv $(SIG)/$${line}".*" $(NONSIG)/; \
	#done < $(RESULTDIR)/00.0_list2move-from-sig2nonsig.tsv

clean :
	- rm -rf *.out
