# USAGE:
# ruby scripts/calculate_expression_means_per_dose.rb expression_data/FINAL_master_list_of_gene_counts_MIN.sense.radiation_renamed.sorted.txt results/dose_response_curves/FINAL_master_list_of_gene_counts_MIN.sense.radiation_means_per_dose.txt 
# ruby scripts/calculate_expression_means_per_dose.rb expression_data/FINAL_master_list_of_gene_counts_MIN.antisense.radiation_renamed.sorted.txt results/dose_response_curves/FINAL_master_list_of_gene_counts_MIN.antisense.radiation_means_per_dose.txt 

# Calculate the expression mean for every dose, for the 2 conditions(radiations).

require 'rubygems'
require 'csv'

ifile = ARGV[0]
ofile = ARGV[1]

# read the list
header = nil
expressions_per_dose = Hash.new { |h,k| h[k] = [] }

CSV.foreach(ifile, {:col_sep => "\t"}) do |row|
	if row[0] == "id"
		header = row
	end
	if row[0] != "" && row[0] != "id"
		for i in 1..20
			expressions_per_dose[header[i]] << row[i].to_i
		end
	end
end

means_per_dose = Hash.new { |h,k| h[k] = [] }
expressions_per_dose.each do |dose, expressions|
	means_per_dose[dose] = (expressions.reduce(:+))/expressions.size
end

# output
CSV.open(ofile, "wb", {:col_sep => "\t"}) do |csv|
	means_per_dose.each do |dose, mean|
		csv << [dose, mean].flatten(1)
	end
end

