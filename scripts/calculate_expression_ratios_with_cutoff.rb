# USAGE:
# ruby scripts/calculate_expression_ratios_with_cutoff.rb results/venn_diagrams/all_sense_proton_0.1.txt results/venn_diagrams/all_sense_gamma_0.1.txt results/expression_ratios_histograms/all_sense_0.1_ratios.txt

# ruby scripts/calculate_expression_ratios_with_cutoff.rb results/venn_diagrams/all_PCR_proton_0.1.txt results/venn_diagrams/all_PCR_gamma_0.1.txt results/expression_ratios_histograms/all_PCR_0.1_ratios.txt

# Calculate the expression ratios for every dose (n=2..10), for every gene with FDR<0.1, for each condition(radiation). expression at dose n/expression at baseline. Create 9 histograms of ratios over all genes per radiation.

require 'rubygems'
require 'csv'

ifile1 = ARGV[0]
ifile2 = ARGV[1]
ofile = ARGV[2]

# read the lists
allgenes1 = Hash.new { |h,k| h[k] = [] }
CSV.foreach(ifile1, {:col_sep => "\t"}) do |row|
	row_to_ratios = []
	if row[0] != "" && row[0] != "ID"
		for i in 6..14	
			row_to_ratios << ((row[i].to_f + 10)/(row[5].to_f + 10)).round(3)
		end
		mean = row_to_ratios[0..8].inject{|sum, x| sum + x}.to_f / row_to_ratios[0..8].size
		allgenes1[row[0]] = [row_to_ratios, row[1..3], mean, row[4]].flatten(2)
	end
end

allgenes2 = Hash.new { |h,k| h[k] = [] }
CSV.foreach(ifile2, {:col_sep => "\t"}) do |row|
	row_to_ratios = []
	if row[0] != "" && row[0] != "ID"
		for i in 6..14	
			row_to_ratios << ((row[i].to_f + 10)/(row[5].to_f + 10)).round(3)
		end
		mean = row_to_ratios[0..8].inject{|sum, x| sum + x}.to_f / row_to_ratios[0..8].size
		allgenes2[row[0]] = [row_to_ratios, row[1..3], mean, row[4]].flatten(2)
	end
end

all_genes_list = allgenes1.keys | allgenes2.keys # concatenate 2 genes lists and remove the duplicates

allgenes = Hash.new { |h,k| h[k] = [] }
all_genes_list.each do |gene|
	if !allgenes1.has_key?(gene)
		for i in 0..13
			allgenes1[gene][i] = "NA" 
		end
	end
	if !allgenes2.has_key?(gene)
		for i in 0..13
			allgenes2[gene][i] = "NA" 
		end
	end
	allgenes[gene] = [allgenes1[gene][0..8], allgenes2[gene][0..8], allgenes1[gene][9..13], allgenes2[gene][9..13]].flatten!(2)
end


protons_header = ["Proton5.S2/baseline", "Proton10.S3/baseline", "Proton25.S4/baseline", "Proton50.S5/baseline", "Proton75.S6/baseline", "Proton100.S7/baseline", "Proton125.S8/baseline", "Proton150.S9/baseline", "Proton200.S10/baseline"]
gammas_header = ["Gamma5.S12/baseline", "Gamma10.S13/baseline", "Gamma25.S14/baseline", "Gamma50.S15/baseline", "Gamma75.S16/baseline", "Gamma100.S17/baseline", "Gamma125.S18/baseline", "Gamma150.S19/baseline", "Gamma200.S20/baseline"]

# output
CSV.open(ofile, "wb", {:col_sep => "\t"}) do |csv|
	csv << ["ID", protons_header, gammas_header, "Proton.FDR", "Proton.direction", "Proton.type", "Proton.mean", "gene", "Gamma.FDR", "Gamma.direction", "Gamma.type", "Gamma.mean", "gene"].flatten!(2)
	allgenes.each do |gene, row|
		csv << [gene, row].flatten!(2)
	end
end

