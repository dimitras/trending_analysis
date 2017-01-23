# USAGE:
# ruby scripts/double_venn.rb trend_data/Gene_Sense_Gamma_trend.txt trend_data/Gene_Sense_Proton_trend.txt results/venn_diagrams/common_genes_from_sense_trend_analysis_0.1.xlsx 0.1
# ruby scripts/double_venn.rb trend_data/Gene_Antisense_Gamma_trend.txt trend_data/Gene_Antisense_Proton_trend.txt results/venn_diagrams/common_genes_from_antisense_trend_analysis_0.1.xlsx 0.1

# ruby scripts/double_venn.rb trend_data/RT-PCR_Gamma_trending.txt trend_data/RT-PCR_Proton_trending.txt results/venn_diagrams/common_genes_from_RT-PCR_trend_analysis_0.1.xlsx 0.1

# Create genes lists with the common genes between 2 DE comparisons.

require 'rubygems'
require 'csv'
require 'axlsx'
require "rinruby" 

strand = "Sense"
# strand = "Antisense"
ifile1 = ARGV[0]
ifile2 = ARGV[1]
ofile = ARGV[2]
cutoff = ARGV[3]

# read the list
counts12_list = Hash.new { |h,k| h[k] = [] }
header = ""
genes_cond1 = Hash.new { |h,k| h[k] = [] }
genes_cond1_list = Hash.new { |h,k| h[k] = [] }
CSV.foreach(ifile1, {:col_sep => "\t"}) do |row|
	if row[0] == "ID"
		header = row
	end
	if row[0] != "ID" && row[0] != ""
		if row[1].to_f <= cutoff.to_f
			genes_cond1[row[0]] << 1
			genes_cond1_list[row[0]] = row
		end
	end
end

genes_cond2 = Hash.new { |h,k| h[k] = [] }
genes_cond2_list = Hash.new { |h,k| h[k] = [] }
CSV.foreach(ifile2, {:col_sep => "\t"}) do |row|
	if row[0] != "ID" && row[0] != "" 
		if row[1].to_f <= cutoff.to_f
			genes_cond2[row[0]] << 1
			genes_cond2_list[row[0]] = row
		end
	end
end

common_genes_list = Hash.new { |h,k| h[k] = [] }
genes_cond1_list.each do |gene, row|
	if genes_cond2_list.has_key?(gene)
		common_genes_list[gene] << [row, genes_cond2_list[gene]].flatten(1)
	end
end

genes_cond1.each do |gene, count|
	if genes_cond2.has_key?(gene)
		counts12_list[gene] << 1
	end
end

total_1 = 0
genes_cond1.each do |gene, count|
	total_1 += genes_cond1[gene].inject(0){|sum, i| sum + i}
end
total_2 = 0
genes_cond2.each do |gene, count|
	total_2 += genes_cond2[gene].inject(0){|sum, i| sum + i}
end

# output
results_xlsx = Axlsx::Package.new
results_wb = results_xlsx.workbook
total_12 = 0

# create genes lists
results_wb.add_worksheet(:name => "#{strand} Commons FDR<#{cutoff}") do |sheet|
	sheet.add_row ["Gamma","","","","","","","","","","","","","","","Proton","","","","","","","","","","","","","",""].flatten(1)
	sheet.add_row [header, "","","","","","","","","", header].flatten(1)
	common_genes_list.each do |gene, row|
		sheet.add_row row.flatten(1)
	end
end
results_wb.add_worksheet(:name => "All #{strand} Gamma FDR<#{cutoff}") do |sheet|
	sheet.add_row header
	genes_cond1_list.each do |gene, row|
		sheet.add_row row
	end
end
results_wb.add_worksheet(:name => "All #{strand} Proton FDR<.#{cutoff}") do |sheet|
	sheet.add_row header
	genes_cond2_list.each do |gene, row|
		sheet.add_row row
	end
end

results_wb.add_worksheet(:name => "#{strand} Gamma-Only FDR<#{cutoff}") do |sheet|
	sheet.add_row header
	genes_cond1_list.each do |gene, row|
		if !common_genes_list.has_key?(gene)
			sheet.add_row row
		end
	end
end
results_wb.add_worksheet(:name => "All #{strand} Proton-Only FDR<.#{cutoff}") do |sheet|
	sheet.add_row header
	genes_cond2_list.each do |gene, row|
		if !common_genes_list.has_key?(gene)
			sheet.add_row row
		end
	end
end


# create list with common peptides between 1,2 conditions
results_wb.add_worksheet(:name => "#{strand} Gamma-Proton") do |sheet|
	sheet.add_row ["gene", "commons", "Gamma count", "Proton count"]
	counts12_list.each do |gene, count|
		sheet.add_row [gene, count.inject(0){|sum,i| sum + i}, genes_cond1[gene].inject(0){|sum,i| sum + i}, genes_cond2[gene].inject(0){|sum,i| sum + i}]
		total_12 += count.inject(0){|sum,i| sum + i}
	end
	sheet.add_row ["totals", total_12, total_1, total_2]
end

# create summary table
results_wb.add_worksheet(:name => "summary") do |sheet|
	sheet.add_row ["total 12", "total 1", "total 2"]
	sheet.add_row [total_12, total_1, total_2]
end

# write xlsx file
results_xlsx.serialize(ofile)

# create venn diagram
R.eval <<EOF
library(VennDiagram)
library(gridExtra)
g = draw.pairwise.venn(area1 = #{total_1}, area2 = #{total_2}, cross.area = #{total_12}, category = c('Gamma', 'Proton'), lty = 'blank', fill = c('skyblue', 'mediumorchid'), fontfamily = 'sans', cex = 4.5, cat.cex = 1, cat.fontfamily = 'sans', cat.pos = c(3,35), cat.dist = c(0.02, 0.03))
pdf('#{ofile}.venn.pdf')
grid.arrange(gTree(children=g), top='genes with FDR<#{cutoff}')
dev.off()
EOF