# USAGE:
# ruby scripts/calculate_expression_means.rb expression_data/FINAL_master_list_of_gene_counts_MIN.sense.radiation_renamed.sorted.txt results/dose_response_curves/FINAL_master_list_of_gene_counts_MIN.sense.radiation_means.txt 
# ruby scripts/calculate_expression_means.rb expression_data/FINAL_master_list_of_gene_counts_MIN.antisense.radiation_renamed.sorted.txt results/dose_response_curves/FINAL_master_list_of_gene_counts_MIN.antisense.radiation_means.txt 

# Calculate the expression mean for every gene, for the 2 conditions(radiations).

require 'rubygems'
require 'csv'

ifile = ARGV[0]
ofile = ARGV[1]

# read the list
header = nil
c1_sums = Hash.new { |h,k| h[k] = [] }
c2_sums = Hash.new { |h,k| h[k] = [] }

CSV.foreach(ifile, {:col_sep => "\t"}) do |row|
	c1_sum = 0
	c2_sum = 0
	if row[0] == "id"
		header = row
	end
	if row[0] != "" && row[0] != "id"
		for i in 1..10
			c1_sum += row[i].to_i
			c2_sum += row[i+10].to_i
		end
		c1_sums[row[0]] = c1_sum
		c2_sums[row[0]] = c2_sum
	end
end

c1_means = Hash.new { |h,k| h[k] = [] }
c1_sums.each do |gene, sum|
	# sum/number of doses(10) for each gene
	c1_means[gene] = (sum/10).to_f
end

c2_means = Hash.new { |h,k| h[k] = [] }
c2_sums.each do |gene, sum|
	# sum/number of doses(10) for each gene
	c2_means[gene] = (sum/10).to_f
end

# output
CSV.open(ofile, "wb", {:col_sep => "\t"}) do |csv|
	csv << ["id", "Proton.expression.mean", "Gamma.expression.mean"].flatten(1)
	Hash[c1_means.sort_by {|_, v| -v}].each do |gene, mean|
		csv << [gene, mean, c2_means[gene]].flatten(1)
	end
end

