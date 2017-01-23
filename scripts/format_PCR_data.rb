# USAGE:
# ruby scripts/format_PCR_data.rb RT-PCR/expression_tables/ trend_data/RT-PCR_Proton_gene_expressions.txt trend_data/RT-PCR_Gamma_gene_expressions.txt

# Format PCR data to a list.

require 'rubygems'
require 'csv'

idir = ARGV[0]
ofile1 = ARGV[1]
ofile2 = ARGV[2]


all_proton = Hash.new { |h,k| h[k] = [] }
all_gamma = Hash.new { |h,k| h[k] = [] }

Dir.foreach(idir) do |ifile|
	next if ifile == "." or ifile == ".."
	
	genesymbol = ifile.split("_")[0]
	tissue = ifile.split("_")[1].split(".")[0]

	raw_table1 = []
	raw_table2 = []
	CSV.foreach(idir+ifile, {:col_sep => "\t"}) do |row|
		if row[0] != "" && row[0] != "Dose"
			raw_table1 << [row[1]]
			raw_table2 << [row[2]]
		end
	end

	transposed_table1 = raw_table1.transpose
	transposed_table2 = raw_table2.transpose
	all_proton[genesymbol+"_"+tissue] = transposed_table1.flatten!(1)
	all_gamma[genesymbol+"_"+tissue] = transposed_table2.flatten!(1)
end

# output
CSV.open(ofile1, "wb", {:col_sep => "\t"}) do |csv|
	csv << ["ID", "Proton0.S1", "Proton5.S2", "Proton10.S3", "Proton25.S4", "Proton50.S5", "Proton75.S6", "Proton100.S7", "Proton125.S8", "Proton150.S9", "Proton200.S10"]
	all_proton.each do |gene, row|
		csv << [gene, row].flatten!(1)
	end
end

# output
CSV.open(ofile2, "wb", {:col_sep => "\t"}) do |csv|
	csv << ["ID", "Gamma0.S11", "Gamma5.S12", "Gamma10.S13", "Gamma25.S14", "Gamma50.S15", "Gamma75.S16", "Gamma100.S17", "Gamma125.S18", "Gamma150.S19", "Gamma200.S20"]
	all_gamma.each do |gene, row|
		csv << [gene, row].flatten!(1)
	end
end