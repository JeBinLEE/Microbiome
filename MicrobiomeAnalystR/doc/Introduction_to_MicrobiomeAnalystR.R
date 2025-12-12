## ----eval=FALSE---------------------------------------------------------------
# # Load MicrobiomeAnalystR
# library(MicrobiomeAnalystR)

## ----eval=FALSE---------------------------------------------------------------
# # Initiate the mbSetObj
# mbSet<-Init.mbSetObj()
# 
# # Set the analysis type
# mbSet<-SetModuleType(mbSet, "mdp")
# 
# # Read in the abundance file, file format is .txt, taxa labels not included in the abundance file,
# # taxonomy labels are Greengenes, metadata file is included separately, and the selected module is 16S.
# mbSet<-Read16SAbundData(mbSet, dataName="ibd_asv_table.txt", format="text", taxa_type="Greengenes", ismetafile="T");
# 
# # Read in the metadata file
# mbSet<-ReadSampleTable(mbSet, dataName="ibd_meta.csv");
# 
# # Read in the taxonomy table
# mbSet<-Read16STaxaTable(mbSet, dataName="ibd_taxa.txt");
# 
# # Read in the phylogenetic tree file
# mbSet<-ReadTreeFile(mbSet, dataName="ibd_tree.tre");
# 
# # Data integrity check
# mbSet<-SanityCheckData(mbSet, filetype="text");
# 
# # Create phyloseq object
# mbSet<-CreatePhyloseqObj(mbSet, type="text", taxa_type="Greengenes", taxalabel="F")

## ----eval=FALSE---------------------------------------------------------------
# # Low abundance filtering, based on the prevalence of features. In this case, by
# # setting the sample percentage to 0.2, only features with >4 counts in at least 20% of all samples
# # will be retained.
# mbSet<-ApplyAbundanceFilter(mbSet, filt.opt="prevalence", count=4, smpl.perc=0.2);
# 
# # Low variance filtering, based on the inter-quantile range. Here, variance is calculated using
# # IQR and 10% of the lowest variance features will be removed.
# mbSet<-ApplyVarianceFilter(mbSet, filtopt="iqr", filtPerct=0.1);

## ----eval=FALSE---------------------------------------------------------------
# # Option 1: No rarefying, total sum scaling and no transformation
# mbSet<-PerformNormalization(mbSet, rare.opt="none", scale.opt="colsum", transform.opt="none");
# 
# # Option 2: Rarefying to the minimum library size + plotting the rarefraction curves
# mbSet<-PerformNormalization(mbSet, rare.opt="rarewi", scale.opt="colsum", transform.opt="none");
# mbSet<-PlotRareCurve(mbSet, graphName="rarefraction_curve.png", variable="X")

## ----eval=FALSE---------------------------------------------------------------
# # Create Biomarker Sweave report
# PreparePDFReport(mSet, "User Name")

