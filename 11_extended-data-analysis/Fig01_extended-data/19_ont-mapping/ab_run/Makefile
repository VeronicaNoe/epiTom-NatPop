SOURCE_DIR = $(abspath $(wildcard *_targets))
FILES = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.mappingOut, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.mappingOut : $(SOURCE_DIR)/%.target
	/mnt/ssd123/vibanez/19_ont-mapping/ab_run/dorado-bam_barcoded-cervil.sh $< 2> $@
	#/mnt/ssd123/vibanez/19_ont-mapping/ab_run/dorado-bam_cervil.sh $< 2> $@

clean :
	- rm -rf *.mappingOut
