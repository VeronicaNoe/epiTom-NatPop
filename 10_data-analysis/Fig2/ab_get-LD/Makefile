SOURCE_DIR = $(abspath ba_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.out, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.out : $(SOURCE_DIR)/%.target
	~/bin/get_mvcf-LD-decay.sh $< 2> $@
	#~/bin/get_LD-decay.sh $< 2> $@

clean :
	- rm -rf *.out
	- rm -rf *.bim
	- rm -rf *.bed
	- rm -rf *.fam
	- rm -rf *.log
	- rm -rf *.irem
	- rm -rf *.nosex


