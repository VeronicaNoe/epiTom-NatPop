SOURCE_DIR = $(abspath ba_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.out, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.out : $(SOURCE_DIR)/%.target
	~/bin/toRun_avg-meth-4-metaplots.sh $< 2> $@

clean :
	- rm -rf *.out

