#!/bin/bash
genome="/mnt/data6/vibanez/IGV/annotation/ITAG2.4_genomic.fasta"
output_xml="/mnt/data6/vibanez/IGV/session_file_RNAseq.xml"
sampleDir="/mnt/data6/vibanez/IGV/data/RNAseq"
sample_file="/mnt/data6/vibanez/IGV/RNAsampleNames.tsv"

# Read sample names from the file into an array
read -a samples <<< $(cat "$sample_file")

# Start building the XML content
xml_content="<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n<Session genome=\"$genome\" version=\"8\">\n\t<Resources>\n"
# Loop through samples and add resource entries
for sample in "${samples[@]}"; do
	echo ${sample}
	xml_content+="\t\t<Resource path=\"$sampleDir/${sample}_leaf_RNAseq.bw\" type=\"bw\"/>\n"
done
# Add annotation files
xml_content+="\t\t<Resource path=\"/mnt/data6/vibanez/IGV/annotation/SL2.50_REPET.sorted.gff3\" type=\"gff3\"/>\n"
xml_content+="\t\t<Resource path=\"/mnt/data6/vibanez/IGV/annotation/ITAG2.4_gene_models.gff3\" type=\"gff3\"/>\n"
xml_content+="\t\t<Resource path=\"/mnt/data6/vibanez/IGV/annotation/dmr_associated.sorted.bed\" type=\"bed\"/>\n"

# Complete the XML content
xml_content+="\t</Resources>\n\t<Panel height=\"4840\" name=\"DataPanel\" width=\"1901\">\n"
# Loop through samples and add track entries
for sample in "${samples[@]}"; do
  xml_content+="\t\t<Track alpha=\"0.4\" attributeKey=\"${sample}_leaf_RNAseq.bw\" autoScale=\"false\" clazz=\"org.broad.igv.track.DataSourceTrack\" color=\"153,153,153\" fontSize=\"8\" id=\"${sampleDir}/${sample}_leaf_RNAseq.bw\" name=\"${sample}\" renderer=\"BAR_CHART\" visible=\"true\" windowFunction=\"mean\">\n"
  xml_content+="\t\t\t<DataRange baseline=\"0.0\" drawBaseline=\"true\" flipAxis=\"false\" maximum=\"1000\" minimum=\"0.0\" type=\"LINEAR\"/>\n"
 # xml_content+="\t\t\t<Property name=\"height\">20</Property>\n"
  xml_content+="\t\t</Track>\n"
done

# Complete the XML content
xml_content+="\t</Panel>\n"
# Add the new tracks to the XML content
xml_content+="\t<Panel height=\"160\" name=\"FeaturePanel\" width=\"1901\">\n"
xml_content+="\t\t<Track attributeKey=\"Reference sequence\" clazz=\"org.broad.igv.track.SequenceTrack\" fontSize=\"10\" id=\"Reference sequence\" name=\"Reference sequence\" sequenceTranslationStrandValue=\"POSITIVE\" shouldShowTranslation=\"false\" visible=\"true\"/>\n"
xml_content+="\t\t<Track attributeKey=\"ITAG2.4_gene_models.gff3\" clazz=\"org.broad.igv.track.FeatureTrack\" fontSize=\"10\" groupByStrand=\"false\" id=\"/mnt/data6/vibanez/IGV/annotation/ITAG2.4_gene_models.gff3\" name=\"ITAG2.4_gene_models.gff3\" visible=\"true\"/>\n"
xml_content+="\t\t<Track attributeKey=\"SL2.50_REPET.sorted.gff3\" clazz=\"org.broad.igv.track.FeatureTrack\" color=\"255,153,0\" fontSize=\"10\" groupByStrand=\"false\" id=\"/mnt/data6/vibanez/IGV/annotation/SL2.50_REPET.sorted.gff3\" name=\"SL2.50_REPET.sorted.gff3\" visible=\"true\"/>\n"
xml_content+="\t\t<Track attributeKey=\"dmr_associated.sorted.bed\" clazz=\"org.broad.igv.track.FeatureTrack\" fontSize=\"10\" groupByStrand=\"false\" id=\"/mnt/data6/vibanez/IGV/annotation/dmr_associated.sorted.bed\" name=\"dmr_associated.sorted.bed\" visible=\"true\"/>\n"
xml_content+="\t</Panel>\n"

xml_content+="\t<PanelLayout dividerFractions=\"0.8221153846153846\"/>\n\t<HiddenAttributes>\n\t\t<Attribute name=\"DATA FILE\"/>\n\t\t<Attribute name=\"DATA TYPE\"/>\n\t\t<Attribute name=\"NAME\"/>\n\t</HiddenAttributes>\n</Session>"
# Save the XML content to the output file
echo -e "$xml_content" > "$output_xml"
