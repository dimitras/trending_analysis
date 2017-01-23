# USAGE:
# ruby scripts/calculate_expression_ratios.rb expression_data/FINAL_master_list_of_gene_counts_MIN.sense.radiation_renamed.sorted.txt results/dose_response_curves/FINAL_master_list_of_gene_counts_MIN.sense.radiation_renamed.sorted_ratios_updated.txt
# ruby scripts/calculate_expression_ratios.rb expression_data/FINAL_master_list_of_gene_counts_MIN.antisense.radiation_renamed.sorted.txt results/dose_response_curves/FINAL_master_list_of_gene_counts_MIN.antisense.radiation_renamed.sorted_ratios_updated.txt

# Calculate the expression ratios for every dose (n=2..10), for every gene, for each condition(radiation). expression at dose n/expression at baseline. Create 9 histograms of ratios over all genes per radiation.

require 'rubygems'
require 'csv'

ifile = ARGV[0]
ofile = ARGV[1]

# read the list
header = nil
allgenes = Hash.new { |h,k| h[k] = [] }
CSV.foreach(ifile, {:col_sep => "\t"}) do |row|
	row_to_ratios = []
	if row[0] == "id"
		header = row
	end
	if row[0] != "" && row[0] != "id"
		for i in 2..10	
			row_to_ratios << ((row[i].to_f + 10)/(row[1].to_f + 10)).round(3)
		end
		for i in 12..20	
			row_to_ratios << ((row[i].to_f + 10)/(row[11].to_f + 10)).round(3)
		end
		c1_mean = row_to_ratios[0..8].inject{|sum, x| sum + x}.to_f / row_to_ratios[0..8].size
		c2_mean = row_to_ratios[9..17].inject{|sum, x| sum + x}.to_f / row_to_ratios[9..17].size
		allgenes[row[0]] = [row_to_ratios, row[21], row[22], c1_mean, c2_mean].flatten(1)
	end
end

protons_header = []
header[2..10].each do |col|
	protons_header << "#{col}/Proton0.S1"
end
gammas_header = []
header[12..20].each do |col|
	gammas_header << "#{col}/Gamma0.S11"
end

# output
CSV.open(ofile, "wb", {:col_sep => "\t"}) do |csv|
	csv << [header[0], protons_header, gammas_header, header[21..24],"Proton.mean", "Gamma.mean"].flatten!(2)
	allgenes.each do |gene, row|
		csv << [gene, row[0..21]].flatten!(1)
	end
end

