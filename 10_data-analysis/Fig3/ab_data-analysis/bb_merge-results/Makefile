SOURCE_DIR = $(abspath ca_targets)
FILES      = $(wildcard $(SOURCE_DIR)/*.target)
TARGETS    = $(patsubst $(SOURCE_DIR)/%.target, $(CURDIR)/%.out, $(FILES))

all : $(TARGETS)

$(CURDIR)/%.out : $(SOURCE_DIR)/%.target
	~/bin/merge_analysis-GWAS.sh $< 2> $@
	for file in *CG-DMR_sigSNPs.mQTL; do
          awk '{OFS="\t"} NR>1 {print $1,$4,$5,$7,$6}' "$file" | sed 's/:/\t/g' >> CG-DMR_sigSNPs.merged
        done
        for file in *C-DMR_sigSNPs.mQTL; do
          awk '{OFS="\t"} NR>1 {print $1,$4,$5,$7,$6}' "$file" | sed 's/:/\t/g' >> C-DMR_sigSNPs.merged
        done
        #chr startSNP endSNP pValues chrDMR  chrDMR_startDMR_DMR
        awk '{OFS="\t"} {printf "%02d\t%s\t%s\t%s\t%s\n", $1, $2, $2+1, $3, "ch"$4"_"$5"_"$6}' CG-DMR_sigSNPs.merged |\
	sortBed -i - |\
	mergeBed -i - -c 4,5,5 -o min,count_distinct,distinct > ../bf_annotation/CG-DMR_sigSNPs.bed

        awk '{OFS="\t"} {printf "%02d\t%s\t%s\t%s\t%s\n", $1, $2, $2+1, $3, "ch"$4"_"$5"_"$6}' C-DMR_sigSNPs.merged |\
	sortBed -i - |\
	mergeBed -i - -c 4,5,5 -o min,count_distinct,distinct > ../bc_annotation/C-DMR_sigSNPs.bed

clean :
	- rm -rf *.out

