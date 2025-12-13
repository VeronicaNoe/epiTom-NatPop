SOURCE_DIR = $(abspath $(wildcard *_targets))
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.split, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.split : $(SOURCE_DIR)/%.target
	bin/split-vcf-group.sh $< 2> $@

clean :
	- rm -rf *.split

