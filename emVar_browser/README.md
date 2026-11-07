# Setup for browsing bovine aortic endothelial cell emVars from SuRE

The .json file is a session file for IGV. It will load the data files from the Zenodo record for the cattle SuRE paper.

To use it, go to [web based IGV](https://igv.org/app/), go to the 'Session' dropdown menu and either find your locally saved copy or use the URL option to point it to the [file on GitHub](https://raw.githubusercontent.com/evotools/SuRE-gg-variants/refs/heads/main/emVar_browser/emVars.json). This will load all the tracks from the Zenodo record.

The image file IGV_screenshot.png shows the view that will appear on loading. This shows a region on chromosome 29 but any location can be navigated to using the location boxes at the top left.

At the top is a feature track giving details of all the cattle variants that were assayed with SuRE. Below are pairs of tracks for each data set where the cattle variants were tested (cattle aortic endothelial cells, human aortic endothelial cells, and the model for an interaction between them. The upper track ('p-value') gives the -log10-transformed nominal p-value in the allele-specific expression model; the lower track ('details') gives further model results, including Benjamini-Hochberg corrected p-value (FDR) and estimated effect size.

Full list of track fields:
SuRE variants track:
Name		Variant ID as given by cactus and used in some data files for the paper
VarType		Variant type (SNP, indel etc)
RefLength	Length in nt of reference allele
AltLength   Length in nt of alternative allele
Gt_Bt		Genotype in Bos taurus SuRE donor animal
Gt_Bi		Genotype in Bos indicus SuRE donor animal
RefCount	Number of fragments assayed that carried the reference allele
AltCount	Number of fragments assayed that carried the alternative allele
AltBases	Sequences of alternative allele (reference sequence can be obtained from the IGV genome)
Location	Variant chromosome and base position

p-value tracks:
Value	-log10-transformed model p-value

details tracks:
Name		Placeholder, all are .
Score		Placeholder, all are 1000
IfSig		Whether variant is an emVar
FDR			Benjamini-Hochberg corrected p-value
Effect		Estimated proportional effect on expression of going from reference to alternative allele
RefScore	Mean normalised SuRE score (see paper Methods) for fragments carrying the reference allele
AltScore	Mean normalised SuRE score for fragments carrying the alternative allele
Skew		Skew of normalised scores for all fragments (if skew>20 were excluded from analysis)
Location	Variant chromosome and base position
