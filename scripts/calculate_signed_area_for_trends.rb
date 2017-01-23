# USAGE:
# ruby scripts/calculate_signed_area_for_trends.rb expression_data/FINAL_master_list_of_gene_counts_MIN.sense.radiation_renamed.sorted.txt results/signed_areas/FINAL_master_list_of_gene_counts_MIN.sense.radiation_signed_areas.txt 
# ruby scripts/calculate_signed_area_for_trends.rb expression_data/FINAL_master_list_of_gene_counts_MIN.antisense.radiation_renamed.sorted.txt results/signed_areas/FINAL_master_list_of_gene_counts_MIN.antisense.radiation_signed_areas.txt 

# ruby scripts/calculate_signed_area_for_trends.rb trend_data/RT-PCR_gene_expressions.txt results/signed_areas/RT-PCR_gene_expressions_signed_areas.txt 

# Calculate the signed area for every gene, for the 2 conditions(radiations).

require 'rubygems'
require 'csv'

ifile = ARGV[0]
ofile = ARGV[1]

# read the list
header = nil
diffs = Hash.new { |h,k| h[k] = [] }
sign_changes = Hash.new { |h,k| h[k] = [] }

CSV.foreach(ifile, {:col_sep => "\t"}) do |row|
	if row[0].casecmp("id") == 0
		header = row
	end
	if row[0] != "" && row[0].casecmp("id") != 0
		sign_changed = nil
		sign = nil
		puts row[0]
		for i in 1..10
			# grab the proton - gamma expression differences per dose
			diffs[row[0]] << row[i].to_i - row[i+10].to_i
			# track sign changes
			dose = header[i].split(".")[0].split(/[\D+]/).join
			if sign_changed.nil? && sign.nil?
				sign_changed = true
				if row[i].to_i - row[i+10].to_i >= 0
					sign = "+"
				elsif row[i].to_i - row[i+10].to_i < 0
					sign = "-"
				end
				sign_changes[row[0]] << [sign, dose]
			elsif sign == "+"
				if row[i].to_i - row[i+10].to_i >= 0
					sign_changed = false
					sign = "+"
				elsif row[i].to_i - row[i+10].to_i < 0
					sign_changed = true
					sign = "-"
					sign_changes[row[0]] << [sign, dose]
				end
			elsif sign == "-"
				if row[i].to_i - row[i+10].to_i >= 0
					sign_changed = true
					sign = "+"
					sign_changes[row[0]] << [sign, dose]
				elsif row[i].to_i - row[i+10].to_i < 0
					sign_changed = false
					sign = "-"
				end
			end
		end
	end
end

max_diff = Hash.new { |h,k| h[k] = [] }
signed_areas = Hash.new { |h,k| h[k] = [] }
diffs.each do |gene, gene_diffs|
	# sum all the diffs 'proton-gamma per dose', per gene
	signed_areas[gene] = gene_diffs.reduce(:+)
	# grab the max diff proton-gamma for each gene
	max_diff[gene] = gene_diffs.map{|x| x.abs}.max
end

total_positives = total_negatives = total_zeros = 0
signed_areas.each do |gene, sum|
	# If sum>0, Proton affects expression more / dominates over Gamma
	# If sum<0, Gamma affects expression more / dominates over Proton
	if sum > 0
		total_positives += 1 
	elsif sum < 0 
		total_negatives += 1
	elsif sum == 0
		total_zeros += 1
	end
end

# output
CSV.open(ofile, "wb", {:col_sep => "\t"}) do |csv|
	csv << ["total positives", total_positives, "total negatives", total_negatives, "total zeros", total_zeros].flatten(1)
	csv << ["id", "signed_area(proton-gamma)", "max(abs(proton-gamma))", "sign changes(sign-dose pairs)"].flatten(1)
	Hash[signed_areas.sort_by {|_, v| -v}].each do |gene, signed_area|
		csv << [gene, signed_area, max_diff[gene], sign_changes[gene].flatten(1)].flatten(1)
	end
end

