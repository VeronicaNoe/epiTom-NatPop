SOURCE_DIR = $(abspath $(wildcard *_targets))
FILES = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.mappingOut, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.mappingOut : $(SOURCE_DIR)/%.target
	#~/work/bin/review_trim-reads-from-cervil.sh $< 2> $@
	~/work/bin/review_mapping-ONT-SL2.5.sh $< 2> $@
clean :
	- rm -rf *.mappingOut
