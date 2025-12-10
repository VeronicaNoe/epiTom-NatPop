SOURCE_DIR = $(abspath ../bc_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.CG-DMR.out, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.CG-DMR.out : $(SOURCE_DIR)/%.target
	~/bin/ko_merge_CG-DMR.sh $< 2> $@

clean :
	- rm -rf *.CG-DMR.out

