for i in xrange(1, end, 1):

print_string = ""
sample_name = str(columns[i][0])
sample_name_re = re.sub('_*1*.fastq.gz', '', sample_name)
print sample_name
print_string = print_string + ">%s\n\n" % sample_name_re
variant_allele = ''.join(columns[i][1:])
print_string = print_string + str(variant_allele)
allele_variant_fasta = open("%s/%s.fa" % (args.filter2_only_snp_vcf_dir, sample_name_re), 'w+')
#print print_string
#print print_string.replace('_.*1.fastq.gz', '')
#print print_string_re
allele_variant_fasta.write(print_string)
allele_variant_fasta.close()
variant_allele_array = []
variant_allele_array.append(columns[i][1:])