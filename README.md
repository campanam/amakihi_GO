# amakihi_GO  
Script to perform gene ontology (GO) enrichment tests for avian malaria challenge experiments in Hawai'i 'amakihi  

Michael G. Campana, 2019  
Smithsonian's National Zoo and Conservation Biology Institute  

## License  
The software is made available under the Smithsonian Institution [terms of use](https://www.si.edu/termsofuse).  

## Citation  
Please cite:  
Paxton, K.L., Cassin-Sackett, L., Atkinson, C.T., Videvall, E., Campana, M.G., Fleischer, R.C. 2023. Gene expression reveals immune response strategies of naïve Hawaiian honeycreepers experimentally infected with introduced avian malaria. *Journal of Heredity*. 114(4): 326–340. [doi: 10.1093/jhered/esad017](https://dx.doi.org/10.1093/jhered/esad017).  

## Installation  
The script requires [Ruby](https://www.ruby-lang.org/). After installation of the language, you can download and install the script using the following commands:  
`git clone https://github.com/campanam/amakihi_GO`  
`mv amakihi_GO/amakihi_GO.rb /some/directory/`  
`cd /some/directory`  

## Usage  
1. Place the gene ontology file (must be named `go.obo`) in the directory where you will execute the `amakihi_GO.rb` script. Gene ontology files are obtainable from [here](http://geneontology.org/docs/download-ontology/) [1-2].  

2. Generate a tab-separated headerless table associating each gene with its GO terms. The first column gives the gene name. The second lists the GO terms, separated by commas. For example:  
`Hmun_k71_000001.g101	GO:0003743,GO:0006413`  
`Hmun_k71_000001.g109	GO:0003735,GO:0005622,GO:0005840,GO:0006412`  
`Hmun_k71_000001.g112	GO:0004930,GO:0007186,GO:0016021`  
`Hmun_k71_000001.g114	GO:0004930,GO:0007186,GO:0016021`  

3. Generate a comma-separated value (CSV) file comparing differentially expressed genes (DEGs) between two treatments. The script expects the gene name to be in the first column and the adjusted P-value to be in the seventh. Other columns are ignored. Any gene with an adjusted P-value less than 0.10 is considered differentially expressed. If a header is used, the script looks for the seventh column to be named "padj". Other values (besides "NA") will be treated as valid data. For example:  
`Amakihi.ID,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj`  
`Hmun_k71_000278.g19718,1527.065763,6.626286975,0.852146564,7.775994479,7.49E-15,1.01E-10`  
`Hmun_k71_000004.g887,155.0433743,-4.485323613,0.588367358,-7.623338638,2.47E-14,1.67E-10`  
`Hmun_k71_000012.g2532,88.95984331,6.220675188,0.838190667,7.421551483,1.16E-13,5.20E-10`  

4. Choose an alpha and an FDR value to evaluate GO-term enrichment significance.  

5. Execute the script using:  
`ruby amakihi_GO.rb <gene_GO_table.tsv> <treatment_comparison.csv> <alpha> <FDR> > <output.tsv>`  

## Output  
The script will generate an output tab-separated values (TSV) file giving the following summary statistics information:  
* TotalDEGs: Total number of DEGs observed in the treatment comparison  
* DEGs with GO Terms: Total number of DEGs that have assigned GO terms  
* Total background with GO Terms: Total number of genes in the database that have assigned GO terms  

It then provides a table of individual GO Terms with the following columns:  
* GOTerm: GO identification number  
* chi2_p-value: P-value for chi-squared test for GO enrichment  
* RawSignificant(alpha=\<VALUE\>): Whether the GO term is significantly enriched based on the chosen alpha  
* BonferroniAdjusted_p-value: P-value after Bonferroni correction  
* BonferroniSignificant: Whether the GO term is significantly enriched after Bonferroni correction  
* Benjamini-Hochberg_CriticalValue: Critical value for enrichment using the Benjamini-Hochberg FDR  
* Benjamini-HochbergSignificant(FDR=\<VALUE\>): Whether the GO term is significantly enriched based on the chosen FDR  
* FoldEnrichment: The GO term's enrichment factor between the treatment and the background  
* GO_Name: Text name of the GO term  
* GO_Namespace: General type (namespace) of the GO term  
* GO_Definition: Text definition of the GO term  
* AllGO-Associated_Genes: List of genes in total database associated with GO term  
* TestGO-Associated_Genes: List of genes in treatment dataset associated with GO term  

## References  
1. Ashburner *et al.* 2000. Gene ontology: tool for the unification of biology. *Nat. Genet.* 25: 25-29.  
2. Gene Ontology Consortium. 2021. The Gene Ontology resource: enriching a GOld mine. *Nucleic Acids Res.* 49: D325-D334.  
