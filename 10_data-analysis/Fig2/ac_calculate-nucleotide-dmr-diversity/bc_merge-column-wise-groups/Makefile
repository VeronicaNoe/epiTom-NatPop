SOURCE_DIR = $(abspath ca_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.out, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.out : $(SOURCE_DIR)/%.target
	~/bin/merge_epiPI-by-features.sh $< 2> $@
	~/bin/merge_PI-by-features.sh $< 2> $@

clean :
	- rm -rf *.out

