#!/bin/bash
genome="/mnt/data6/vibanez/IGV/annotation/ITAG2.4_genomic.fasta"
output_xml="/mnt/data6/vibanez/IGV/generated_session_file.xml"
sampleDir="/mnt/data6/vibanez/IGV/data"
sample_file="/mnt/data6/vibanez/IGV/samplesName.tsv"
context="CG CHG CHH"

# Read sample names from the file into an array
read -a samples <<< $(cat "$sample_file")

# Start building the XML content
xml_content="<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n<Session genome=\"$genome\" version=\"8\">\n<Resources>\n"

# Loop through samples and add resource entries
for sample in "${samples[@]}"; do
  echo ${sample}
  for context_sample in $context; do
    bw_file="$sampleDir/${sample}_leaf_Biseq"_${context_sample}".bw"
    xml_content+="\t<Resource path=\"$bw_file\" type=\"bw\"/>\n"
done
# Complete the XML content
xml_content+="</Resources>\n<Panel height=\"778\" name=\"DataPanel\" width=\"1901\">\n"
# Loop through samples and add track entries
for sample in "${samples[@]}"; do
  xml_content+="\t<Track alpha=\"0.24\" attributeKey=\"Overlay\" autoScale=\"false\" clazz=\"org.broad.igv.track.MergedTracks\" fontSize=\"10\" id=\"${sampleDir}/${sample}_leaf_Biseq\" name=\"${sample}\" renderer=\"BAR_CHART\" visible=\"true\">\n"
  ## add color per CG
  xml_content+="\t\t<Track altColor=\"204,204,0\" attributeKey=\"${sample}_leaf_Biseq.bw\" autoScale=\"false\" color=\"${color}\" fontSize=\"10\" id=\"${sampleDir}/${sample}_leaf_Biseq_CG.bw\" name=\"${sample}_leaf_Biseq_CG.bw\" renderer=\"BAR_CHART\" visible=\"true\" windowFunction=\"mean\"/>\n"
  ## add color per CHG
  xml_content+="\t\t<Track altColor=\"0,51,204\" attributeKey=\"${sample}_leaf_Biseq.bw\" autoScale=\"false\" color=\"${color}\" fontSize=\"10\" id=\"${sampleDir}/${sample}_leaf_Biseq_CHG.bw\" name=\"${sample}_leaf_Biseq_CHG.bw\" renderer=\"BAR_CHART\" visible=\"true\" windowFunction=\"mean\"/>\n"
  ## add color per CHH
  xml_content+="\t\t<Track altColor=\"204,0,102\" attributeKey=\"${sample}_leaf_Biseq.bw\" autoScale=\"false\" color=\"${color}\" fontSize=\"10\" id=\"${sampleDir}/${sample}_leaf_Biseq_CHH.bw\" name=\"${sample}_leaf_Biseq_CHH.bw\" renderer=\"BAR_CHART\" visible=\"true\" windowFunction=\"mean\"/>\n"
  xml_content+="\t\t<DataRange baseline=\"0.0\" drawBaseline=\"true\" flipAxis=\"false\" maximum=\"1.0\" minimum=\"0.0\" type=\"LINEAR\"/>\n"
  xml_content+="\t</Track>\n"
done
# Complete the XML content
xml_content+="</Panel>\n<Panel height=\"47\" name=\"FeaturePanel\" width=\"1901\">\n"
xml_content+="\t<Track attributeKey=\"Reference sequence\" clazz=\"org.broad.igv.track.SequenceTrack\" fontSize=\"10\" id=\"Reference sequence\" name=\"Reference sequence\" sequenceTranslationStrandValue=\"POSITIVE\" shouldShowTranslation=\"false\" visible=\"true\"/>\n"
xml_content+="</Panel>\n<PanelLayout dividerFractions=\"0.9375\"/>\n<HiddenAttributes>\n\t<Attribute name=\"DATA FILE\"/>\n\t<Attribute name=\"DATA TYPE\"/>\n\t<Attribute name=\"NAME\"/>\n</HiddenAttributes>\n</Session>"
# Save the XML content to the output file
echo -e "$xml_content" > "$output_xml"

# Launch IGV with the generated session file
igv -g $genome -b $output_xml -o
