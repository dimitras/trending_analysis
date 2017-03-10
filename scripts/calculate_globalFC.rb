# USAGE:
# ruby scripts/calculate_globalFC.rb trend_data/Gene_Sense_Gamma_trend_sorted.txt trend_data/Gene_Sense_Proton_trend_sorted.txt results/globalFCs.xlsx

# Filter by different FDR and calculate the global FC for each Radiation type.
require 'rubygems'
require 'csv'
require 'axlsx'
require "rinruby" 

ifile1 = ARGV[0]
ifile2 = ARGV[1]
ofile = ARGV[2]
fdr_cutoffs = [0.5, 0.4, 0.3, 0.2, 0.1, 0.05]

# read the lists

puts "******************* #{ifile1} *******************"

avg_fc_per_row = Hash.new { |h,k| h[k] = [] }
CSV.foreach(ifile1, {:col_sep => "\t"}) do |row|
	fdr_cutoffs.each do |fdr|
		if (row[0] != "ID") && (row[0] != "")
			if row[1].to_f <= fdr
				row_to_ratios = []
				for i in 6..14
					# puts "Counts per row: #{row[i]} over #{row[5]} turns to #{((row[i].to_f + 10)/(row[5].to_f + 10)).round(3)}"
					row_to_ratios << ((row[i].to_f + 10)/(row[5].to_f + 10)).round(3) # calculate the FCs for each dose over baseline
					# puts "FCs per row: #{row_to_ratios.join(" + ")}"
				end
				avg_fc_per_row[fdr] << row_to_ratios[0..8].inject{|sum, x| sum + x}.to_f / row_to_ratios[0..8].size # push the avg FC per gene
				# puts "AVG FC per row: #{avg_fc_per_row[fdr].last}"
			end
		end
	end
end

puts "******************* #{ifile2} *******************"

avg_fc_per_row_c2 = Hash.new { |h,k| h[k] = [] }
CSV.foreach(ifile2, {:col_sep => "\t"}) do |row|
	fdr_cutoffs.each do |fdr|
		if (row[0] != "ID") && (row[0] != "")
			if row[1].to_f <= fdr
				row_to_ratios = []
				for i in 6..14
					# puts "Counts per row: #{row[i]} over #{row[5]} turns to #{((row[i].to_f + 10)/(row[5].to_f + 10)).round(3)}"
					row_to_ratios << ((row[i].to_f + 10)/(row[5].to_f + 10)).round(3) # calculate the FCs for each dose over baseline
					# puts "FCs per #{row[0]}: #{row_to_ratios.join(" + ")}"
				end
				avg_fc_per_row_c2[fdr] << row_to_ratios[0..8].inject{|sum, x| sum + x}.to_f / row_to_ratios[0..8].size
				# puts "AVG FC per #{row[0]} at FDR<=#{fdr}: #{avg_fc_per_row_c2[fdr].last}"
			end
		end
	end
end

puts "******************* #{ofile} *******************"

# output
results_xlsx = Axlsx::Package.new
results_wb = results_xlsx.workbook

# fdr_cutoffs.each do |fdr|
# 	puts "Global FC for cond2 at FDR<#{fdr}: #{avg_fc_per_row_c2[fdr].join(" + ")} =\n #{avg_fc_per_row_c2[fdr].inject{|sum, x| sum + x}}"
# end
# global FCs
results_wb.add_worksheet(:name => "global FCs") do |sheet|
	sheet.add_row ["FDR","GammaFC","ProtonFC"].flatten(1)
	fdr_cutoffs.each do |fdr|
		sheet.add_row [fdr, avg_fc_per_row[fdr].inject{|sum, x| sum + x}.to_f / avg_fc_per_row[fdr].size, avg_fc_per_row_c2[fdr].inject{|sum, x| sum + x}.to_f / avg_fc_per_row_c2[fdr].size].flatten(1)
	end
end

# write xlsx file
results_xlsx.serialize(ofile)