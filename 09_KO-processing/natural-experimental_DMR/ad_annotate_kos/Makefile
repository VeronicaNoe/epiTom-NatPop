SOURCE_DIR = $(abspath ba_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.out, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.out : $(SOURCE_DIR)/%.target
	~/bin/annotate_KOs.sh $< 2> $@
	#~/bin/descriptive_stat_annotate-cmt3-RdDM.sh $< 2> $@


clean :
	- rm -rf *.out

