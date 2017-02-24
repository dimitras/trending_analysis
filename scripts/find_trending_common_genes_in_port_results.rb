# USAGE:
# ruby scripts/find_trending_common_genes_in_port_results.rb results/venn_diagrams/19_common_genes2symbols.csv expression_data/FINAL_master_list_of_gene_counts_MIN.sense.radiation_renamed.sorted.txt results/venn_diagrams/common.genes.sense.0.1.PORTvalues.txt

# Find the 19 common genes found from the venn diagram, in the PORT results.

require 'rubygems'
require 'csv'

ifile1 = ARGV[0]
ifile2 = ARGV[1]
ofile = ARGV[2]

# read the lists
genes = Hash.new { |h,k| h[k] = [] }
CSV.foreach(ifile1) do |row|
	if row[0] != "" && row[0] != "ID"
		genes[row[0]] = row[1]
	end
end

retrieved_genes = Hash.new { |h,k| h[k] = [] }
header = nil
CSV.foreach(ifile2, {:col_sep => "\t"}) do |row|
	if row[0] == "id"
		header = row[0..20]
	end
	if row[0] != "" && row[0] != "id" && genes.has_key?(row[0])
		retrieved_genes[row[0]] = row[1..20]
	end
end

# output
CSV.open(ofile, "wb", {:col_sep => "\t"}) do |csv|
	csv << header
	retrieved_genes.each do |gene, row|
		csv << [genes[gene], row].flatten!(1)
	end
end
