#!/usr/bin/env ruby
$genes = {} # Hash mapping genes to GO terms
$goterms_genes = {} # Hash mapping GO terms to genes for background
$compgoterms_genes = {} # Hash mapping GO terms to genes for test set
$totalref = 0 # Total number of genes in ref set
$freqs = {} # Hash mapping reference GO terms to counts
$totalcomp = 0 # Total number of genes in in comp set with GO terms
$totaldegs = 0 # Total number of DEGs in comp set
$compfreqs = {} # Hash mapping comparison GO terms to counts
$goterms = {} # Hash of GoTerms
#--------------------------------------------------------------
class GoTerm
	attr_accessor	:id, :name, :namespace, :definition
	def initialize(id)
		@id = id
		@name = @namspace = @definition = ""
	end
end
#--------------------------------------------------------------
def read_go
	currentgo = nil
	alt_ids = [] # Array of alt_id
	File.open("go.obo") do |f1|
		while line = f1.gets
			if line[0..2] == "id:"
				$goterms[currentgo.id] = currentgo unless currentgo.nil?
				for alt_id in alt_ids
					$goterms[alt_id] = currentgo
				end
				alt_ids.clear
				currentgo = GoTerm.new(line[4...-1])
			elsif line[0..4] == "name:"
				currentgo.name = line[6...-1]
			elsif line[0..9] == "namespace:"
				currentgo.namespace = line[11...-1]
			elsif line[0..3] == "def:"
				currentgo.definition = line[5...-1]
			elsif line[0..6] == "alt_id:"
				alt_ids.push(line[8...-1])
			end
		end
	end
	$goterms[currentgo.id] = currentgo
	for alt_id in alt_ids
		$goterms[alt_id] = currentgo
	end
end
#--------------------------------------------------------------
def read_ref_data(ref_data)
	File.open(ref_data) do |f1|
		while line = f1.gets
			$totalref += 1
			line_arr = line[0...-1].split("\t")
			go_arr = line_arr[1].split(",")
			$genes[line_arr[0]] = go_arr
			for go in go_arr
				if $freqs.keys.include?(go)
					$freqs[go] += 1
				else
					$freqs[go] = 1
				end
				if $goterms_genes.keys.include?(go) # Map GO Terms back to genes
					$goterms_genes[go].push(line_arr[0])
				else
					$goterms_genes[go] = [line_arr[0]]
				end
			end
		end
	end					
end
#--------------------------------------------------------------
def read_comp_data(comp_data)
	File.open(comp_data) do |f2|
		while line = f2.gets
			gene_arr = line[0...-1].split(",")
			go_arr = $genes[gene_arr[0]]
			$totaldegs += 1 if (gene_arr[6].to_f < 0.10 && gene_arr[6] != "NA" && gene_arr[6] != "padj")
			unless go_arr.nil?
				if (gene_arr[6].to_f < 0.10 && gene_arr[6] != "NA")
					$totalcomp += 1
					for go in go_arr
						if $compfreqs.keys.include?(go)
							$compfreqs[go] += 1
						else
							$compfreqs[go] = 1
						end
						if $compgoterms_genes.keys.include?(go) # Map GO Terms back to genes
							$compgoterms_genes[go].push(gene_arr[0])
						else
							$compgoterms_genes[go] = [gene_arr[0]]
						end
					end
				end
			end
		end
	end
end
#--------------------------------------------------------------
def factorial(val)
	total = 1
	while val > 1
		total *= val
		val -= 1
	end
	return total
end
#--------------------------------------------------------------
def fisher(go, fcomp, fref)
	comp_go = $compfreqs[go]
	ref_go = $freqs[go]
	not_comp_go = $totalcomp - comp_go
	not_ref_go = $totalref - ref_go
	denom = factorial(ref_go)*factorial(comp_go)*factorial(not_ref_go)*factorial(not_comp_go)
	numer = factorial(ref_go+comp_go)*factorial(not_comp_go+not_ref_go)*fcomp*fref
	pval = numer.to_f/denom.to_f
	puts go + "\t" + pval.to_s
end
#--------------------------------------------------------------
def chi_cum_prob(test) # This calculates chi distribution cumulative probability distribution under special case of df = 1
	return Math.erf(Math.sqrt(test/2.0))
end
#--------------------------------------------------------------
def chi2(go)
	comp_go = $compfreqs[go]
	ref_go = $freqs[go]
	not_comp_go = $totalcomp - comp_go
	not_ref_go = $totalref - ref_go
	numer = ((ref_go*not_comp_go - comp_go*not_ref_go)**2)*(ref_go+not_comp_go+comp_go+not_ref_go)
	denom = (ref_go+comp_go)*(not_ref_go+not_comp_go)*(comp_go+not_comp_go)*(ref_go+not_ref_go)
	chi2 = numer.to_f/denom.to_f
	return 1.0 - chi_cum_prob(chi2)
end
#--------------------------------------------------------------
def enrichment(go)
	enrichment = ($compfreqs[go].to_f * $totalref.to_f)/($freqs[go].to_f *  $totalcomp.to_f) # Division of fractions
	return enrichment
end
#--------------------------------------------------------------
def benjamini_hochberg(pvalues)
	rank = maxrank = 1.0
	for gene in  pvalues
		crit = (rank/$totalcomp.to_f)*ARGV[3].to_f
		maxrank = rank if gene[1] <= crit
		rank += 1.0
		gene.push(crit)
	end
	rank = 1.0
	for gene in pvalues
		rank <= maxrank ? gene.push("Yes") : gene.push("No")
		rank += 1.0
	end
	return pvalues
end
#--------------------------------------------------------------
read_go
read_ref_data(ARGV[0])
read_comp_data(ARGV[1])
fcomp = factorial($totalcomp)
fref = factorial($totalref)
$chi2_p = {} # Hash of chi2 p-values
for go in $compfreqs.keys
	#fisher(go, fcomp, fref) # Values too large and computer rounds factorials to Infinity
	$chi2_p[go] = chi2(go)
end
$chi2_p = $chi2_p.sort_by { |key, value| value }
$chi2_p = benjamini_hochberg($chi2_p)
for go in $chi2_p
	enrich = enrichment(go[0])
	go.push(enrich)
end
puts "SummaryStatistics"
puts "Total DEGs: " + $totaldegs.to_s
puts "DEGs with GO Terms: " + $totalcomp.to_s
puts "Total background with GO Terms: " + $totalref.to_s
puts "\nGOTerm\tchi2_p-value\tRawSignificant(alpha=#{ARGV[2]})\tBonferroniAdjusted_p-value\tBonferroniSignificant\tBenjamini-Hochberg_CriticalValue\tBenjamini-HochbergSignificant(FDR=#{ARGV[3]})\tFoldEnrichment\tGO_Name\tGO_Namespace\tGO_Definition\tAllGO-Associated_Genes\tTestGO-Associated_Genes"
for gene in $chi2_p
	out = gene[0] + "\t" + gene[1].to_s + "\t"
	gene[1] <= ARGV[2].to_f ? out << "Yes" : out << "No"
	cpval = gene[1]*$totalcomp.to_f # Bonferroni corrected pvalue
	cpval = 1.0 if cpval > 1.0
	out << "\t" + cpval.to_s + "\t"
	cpval <= ARGV[2].to_f ? out << "Yes" : out << "No"
	out << "\t" + gene[2].to_s + "\t" + gene[3] + "\t" + gene[4].to_s + "\t" + $goterms[gene[0]].name + "\t" + $goterms[gene[0]].namespace + "\t" + $goterms[gene[0]].definition + "\t" + $goterms_genes[gene[0]].join(",") + "\t" + $compgoterms_genes[gene[0]].join(",")
	puts out
end
