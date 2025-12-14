SOURCE_DIR = $(abspath $(wildcard *_targets))
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.plinkformat, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.plinkformat : $(SOURCE_DIR)/%.target
	bin/toRun_plinkFormat.sh $< 2> $@

clean :
	- rm -rf *.plinkformat
