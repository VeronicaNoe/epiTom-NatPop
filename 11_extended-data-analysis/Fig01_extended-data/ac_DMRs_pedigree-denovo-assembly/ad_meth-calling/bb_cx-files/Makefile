SOURCE_DIR = $(abspath $(wildcard *_targets))
FILES = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.mappingOut, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.mappingOut : $(SOURCE_DIR)/%.target
	~/work/bin/review_02_modkit_to_cx.sh  $< 2> $@

clean :
	- rm -rf *.mappingOut
