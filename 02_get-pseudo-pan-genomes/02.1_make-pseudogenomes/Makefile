SOURCE_DIR = $(abspath aa_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.genomePrepLog, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.genomePrepLog : $(SOURCE_DIR)/%.target
	~/bin/making_pseudoGenome.sh $^ > $@

clean :
	- rm -rf *.genomePrepLog
