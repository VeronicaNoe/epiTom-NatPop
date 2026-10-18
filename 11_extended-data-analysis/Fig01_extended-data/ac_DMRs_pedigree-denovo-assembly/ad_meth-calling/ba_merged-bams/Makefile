SOURCE_DIR = $(abspath $(wildcard *_targets))
FILES = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.mappingOut, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.mappingOut : $(SOURCE_DIR)/%.target
	~/work/bin/review_01_prepare_final_modbam.sh $< 2> $@

clean :
	- rm -rf *.mappingOut
