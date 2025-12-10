#!/users/bioinfo/vibanez/anaconda3/bin/python
import os
import fnmatch  # Import the fnmatch module
import sys
args = sys.argv[1:]  # The first element in sys.argv is the script filename, so we skip it
# Now, args contains the command-line arguments provided when running the script
# Define the working directory and input arguments
WD = "/mnt/disk2/vibanez/10_data-analysis/Fig3/aa_GWAS-DMRs/bd_results/nonSig"
OUTDIR = "/mnt/disk2/vibanez/10_data-analysis/Fig3/ab_data-analysis/bb_merge-results"
os.chdir(WD)
toProcess = args[0] + '_' + args[1]
print(toProcess)
# Get a list of sample files ending with '.mQTL'
samples = [f for f in os.listdir(WD) if fnmatch.fnmatch(f, f"{toProcess}*.mQTL")]
samples = [os.path.splitext(f)[0] for f in samples]
mQTL = []
colnames = [ 'logVariance', 'log', 'delta', 'sigmaGen', 'sigmaError', 'h2', 'chrDMR', 'startDMR', 'DMR']
for i in samples:
        print(i)
        chrDMR = i.split('_')[0]
        DMR = i.split('_')[1]
        startDMR = i.split('_')[2]
        metadata = [chrDMR, startDMR, DMR]
        with open(os.path.join(WD, i + '.reml'), 'r') as reml_file:
                reml_lines = reml_file.readlines()
                reml_values = [line.strip().split('\t')[0] for line in reml_lines]
                reml_values.extend(metadata)
                mQTL.append(reml_values)
                # Define the column names
#                colnames = [ 'logVariance', 'log', 'delta', 'sigmaGen', 'sigmaError', 'h2', 'chrDMR', 'startDMR', 'DMR']

# Write the data to a '.tsv' file
with open(os.path.join(OUTDIR, f"{args[0]}_{args[1]}_nonSigSNPs.reml"), 'w') as output_file:
        output_file.write('\t'.join(colnames) + '\n')
        for row in mQTL:
                output_file.write('\t'.join(map(str, row)) + '\n')
