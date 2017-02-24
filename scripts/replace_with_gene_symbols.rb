# USAGE:
# ruby scripts/replace_with_gene_symbols.rb RT-PCR/id2info.txt trend_data/RT-PCR_gene_expressions_edited.txt results/dose_response_curves/RT-PCR_gene_expressions_edited_w.genesymbols.txt

# Convert to gene symbols in PCR data.

require 'rubygems'
require 'csv'

ifile1 = ARGV[0]
ifile2 = ARGV[1]
ofile = ARGV[2]

# read the lists
genesymbols = Hash.new { |h,k| h[k] = [] }
CSV.foreach(ifile1, {:col_sep => "\t"}) do |row|
	if row[0] != ""
		genesymbols[row[0]] = row[1]
	end
end

genes = Hash.new { |h,k| h[k] = [] }
header = nil
CSV.foreach(ifile2, {:col_sep => "\t"}) do |row|
	if row[0] == "ID"
		header = row[0..20]
	end
	if row[0] != "" && row[0] != "ID" #&& genes.has_key?(row[0])
		genes[row[0]] = row[1..20]
	end
end

# output
CSV.open(ofile, "wb", {:col_sep => "\t"}) do |csv|
	csv << header
	genes.each do |gene, row|
		csv << [genesymbols[gene], row].flatten!(1)
	end
end
