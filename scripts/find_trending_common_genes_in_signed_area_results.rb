# USAGE:
# ruby scripts/find_trending_common_genes_in_signed_area_results.rb results/venn_diagrams/19_common_genes2symbols.csv results/signed_areas/RNA-Seq_signed_areas_reformatted.txt results/signed_areas/signed_areas_for_19_common_genes.txt


# Find the 19 common genes found from the venn diagram, in the signed areas results.

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
		header = row
	end
	if row[0] != "" && row[0] != "id" && genes.has_key?(row[0])
		retrieved_genes[row[0]] = row[1]
	end
end

# output
CSV.open(ofile, "wb", {:col_sep => "\t"}) do |csv|
	csv << header
	retrieved_genes.each do |gene, row|
		csv << [genes[gene], row]
	end
end
