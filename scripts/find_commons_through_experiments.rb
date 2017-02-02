# USAGE:
# ruby scripts/find_commons_through_experiments.rb trend_data/Gene_Sense_Gamma_trend_sorted_filtered.txt trend_data/Gene_Sense_Proton_trend_sorted_filtered.txt trend_data/RT-PCR_Gamma_trending_sorted_filtered.txt trend_data/RT-PCR_Proton_trending_sorted_filtered.txt results/expression_ratios_histograms/common_genes_across_experiments_0.1.txt

# Find the common genes across multiple experiments.

require 'rubygems'
require 'csv'

ifile1 = ARGV[0]
ifile2 = ARGV[1]
ifile3 = ARGV[2]
ifile4 = ARGV[3]
ofile = ARGV[4]

# read the lists
genes1 = Hash.new { |h,k| h[k] = [] }
CSV.foreach(ifile1, {:col_sep => "\t"}) do |row|
	if row[0] != "" && row[0] != "ID" && !row[5].nil?
		rsum = 0
		for i in 5..14
			rsum += row[i].to_f
		end
		genes1[row[0]] = [row[0..1], (rsum/10).round(3), row[5..14], row[4]].flatten(2)
	end
end

genes2 = Hash.new { |h,k| h[k] = [] }
CSV.foreach(ifile2, {:col_sep => "\t"}) do |row|
	if row[0] != "" && row[0] != "ID" && !row[5].nil?
		rsum = 0
		for i in 5..14	
			rsum += row[i].to_f
		end
		genes2[row[0]] = [row[0..1], (rsum/10).round(3), row[5..14]].flatten(2)
	end
end

genes3 = Hash.new { |h,k| h[k] = [] }
CSV.foreach(ifile3, {:col_sep => "\t"}) do |row|
	if row[0] != "" && row[0] != "ID" && !row[5].nil? && (!row[0].include? "RNA-Seq")
		rsum = 0
		for i in 5..14	
			rsum += row[i].to_f
		end
		genes3[row[0].split("_")[0]] = [row[0..1], (rsum/10).round(3), row[5..14]].flatten(2)
	end
end

genes4 = Hash.new { |h,k| h[k] = [] }
CSV.foreach(ifile4, {:col_sep => "\t"}) do |row|
	if row[0] != "" && row[0] != "ID" && !row[5].nil? && (!row[0].include? "RNA-Seq")
		rsum = 0
		for i in 5..14	
			rsum += row[i].to_f
		end
		genes4[row[0].split("_")[0]] = [row[0..1], (rsum/10).round(3), row[5..14]].flatten(2)
	end
end



common_genes = Hash.new { |h,k| h[k] = [] }
genes1.each do |gene, row|
	if genes2.has_key?(gene) && genes3.has_key?(gene) && genes4.has_key?(gene) 
		common_genes[gene] = [gene, row[1], genes2[gene][1], genes3[gene][1], genes4[gene][1], row[2], genes2[gene][2], genes3[gene][2], genes4[gene][2], row[3..12], genes2[gene][3..12], genes3[gene][3..12], genes4[gene][3..12], row[13]].flatten!(2)
	end
end

header = ["RNA-Seq.Gamma.FDR", "RNA-Seq.Proton.FDR", "PCR.Gamma.FDR", "PCR.Proton.FDR", "RNA-Seq.Gamma.mean", "RNA-Seq.Proton.mean", "PCR.Gamma.mean", "PCR.Proton.mean"]
protons_header = ["Proton0.S1", "Proton5.S2", "Proton10.S3", "Proton25.S4", "Proton50.S5", "Proton75.S6", "Proton100.S7", "Proton125.S8", "Proton150.S9", "Proton200.S10"]
gammas_header = ["Gamma0.S11", "Gamma5.S12", "Gamma10.S13", "Gamma25.S14", "Gamma50.S15", "Gamma75.S16", "Gamma100.S17", "Gamma125.S18", "Gamma150.S19", "Gamma200.S20"]

a = "RNA."
b = "PCR."
# output
CSV.open(ofile, "wb", {:col_sep => "\t"}) do |csv|
	csv << ["ID", header, gammas_header.map{|x| a+x}, protons_header.map{|x| a+x}, gammas_header.map{|x| b+x}, protons_header.map{|x| b+x}, "gene"].flatten!(2)
	common_genes.each do |gene, row|
		csv << [row].flatten!(2)
	end
end

