#!/users/bioinfo/vibanez/anaconda3/bin/python
import os
import fnmatch  # Import the fnmatch module
import sys
args = sys.argv[1:]  # The first element in sys.argv is the script filename, so we skip it
# Now, args contains the command-line arguments provided when running the script
# Define the working directory and input arguments
#WD = "/mnt/data6/vibanez/SNPs/vcfiles/ah_GWAS/analysis_result/bc_snp-blocks/bc_result"
WD = "/mnt/disk2/vibanez/10_data-analysis/Fig3/aa_GWAS-DMRs/bd_results/sig"

OUTDIR = "/mnt/disk2/vibanez/10_data-analysis/Fig3/ab_data-analysis/bb_merge-results"
os.chdir(WD)
toProcess = args[0] + '_' + args[1]
print(toProcess)
# Get a list of sample files ending with '.mQTL'
samples = [f for f in os.listdir(WD) if fnmatch.fnmatch(f, f"{toProcess}*.mQTL")]
samples = [os.path.splitext(f)[0] for f in samples]
mQTL = []
for i in samples:
	#print(i)
	chrDMR = i.split('_')[0].replace('ch', '')
	DMR = i.split('_')[1]
	startDMR = i.split('_')[2]
	metadata = [chrDMR, startDMR, DMR]
	with open(os.path.join(WD, i + '.mQTL'), 'r') as f:
		lines = f.readlines()
		for line in lines:
			parts = line.strip().split('\t')
			qtlTmp = [parts[0], parts[1], parts[2],parts[3]]
			qtlTmp.extend(metadata)
			mQTL.append(qtlTmp)
		colnames = ['SNP', 'effect', 'effectSE','pValues','chrDMR','startDMR', 'DMR']
		#colnames = ['chrSNP', 'startSNP','endSNP','pValue', 'effect', 'effectSE', 'chrDMR', 'startDMR', 'DMR']


# Write the data to a '.tsv' file
with open(os.path.join(OUTDIR, f"{args[0]}_{args[1]}_{args[2]}.mQTL"), 'w') as output_file:
	output_file.write('\t'.join(colnames) + '\n')
	for row in mQTL:
		output_file.write('\t'.join(map(str, row)) + '\n')
