#!/home/vibanez/anaconda3/bin/python
# read in gzip file
import gzip
import sys, argparse
import os
os.chdir("/home/vibanez/results")
outdir =  "/mnt/disk2/vibanez/GWAS/analysis-Results/sigSNPs/"
fileName = sys.argv[1]
prefix = sys.argv[2]
f = gzip.open(fileName, 'r')
f.readline()
#sig_snp = open(outdir + prefix + ".mQTL", "w")
output_filename = outdir + prefix + ".mQTL"
#write header
#sig_snp.write("chr:PosSNP\tsig\tdelta\tpValue\n")
#sig_snp.write()
result_dict = {'significant_lines': []}
counter = 0
for line in f:
			line = line.strip().split()
			values = line[3]
			significant_count = 0
			# check value in line
			for value in values:
				if float(values) <=0.00000005:
					#sig_snps.write(line)
					result_dict['significant_lines'].append(line)
					print("significant")

			counter += 1
			if (counter % 10000) == 0:
							print ("Read %d entires." % counter)

if result_dict['significant_lines']:
        with open(output_filename, 'w') as output_file:
                output_file.write("chr:PosSNP\tsig\tdelta\tpValue\n")

                for line_bytes in result_dict['significant_lines']:
                        val = [byte.decode('utf-8') for byte in line_bytes]
                        formatted_line = '\t'.join(val) + '\n'
                        #line_str = '\t'.join(map(str, line_list))
                        #output_file.write(line_str + '\n')
                        output_file.write(formatted_line)


#sig_snp.close()
f.close()
