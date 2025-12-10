SOURCE_DIR = $(abspath ba_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.out, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.out : $(SOURCE_DIR)/%.target
	~/bin/get_mSFS.sh $< 2> $@
	~/bin/get_SFS.sg $< 2> $@
	~/bin/get_SFS_group.sg $< 2> $@

clean :
	- rm -rf *.out

