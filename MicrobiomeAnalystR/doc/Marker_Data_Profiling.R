## ----eval=FALSE---------------------------------------------------------------
# # Stacked bar plot, percentage abundance, genus level, viridis color palette
# mbSet<-PlotTaxaAundanceBar(mbSet, "taxa_alpha_1","Genus","Group", "none", "barraw",10, "set3","sum",10, "bottom", "F", "png");
# 
# # Stacked area plot, class level, viridis color palette
# mbSet<-PlotTaxaAbundanceArea(mbSet, "taxa_alpha_3","overview","Genus","Group", 10, "viridis","sum",10,"bottom", "F", "png");

## ----eval=FALSE---------------------------------------------------------------
# # Overall pie chart
# mbSet<-PlotOverallPieGraph(mbSet, "Phylum", 10,"sum", 10, "bottom");
# GetSeriesColors()
# mbSet<-SavePiechartImg(mbSet, "Phylum","primary_piechart_0","png");

## ----eval=FALSE---------------------------------------------------------------
# # Heat tree at the genus level comparing pediatric CD vs healthy controls. Only
# # significant taxon (p < 0.05) are labelled.
# mbSet<-PrepareHeatTreePlot(mbSet, meta="Class", taxalvl="Genus", color="plasma", layoutOpt="reda", comparison="CD_vs_Control",
#                           wilcox.cutoff = 0.05, imgName="heat_tree_0", format="png", dpi=300)

## ----eval=FALSE---------------------------------------------------------------
# # Create a dot plot of alpha diversity measures. Samples are grouped by the "Class" experimental factor and colored
# # using the viridis color palette.
# mbSet<-PlotAlphaData(mbSet, "filt","alpha_diver_0","Chao1","Class","OTU", "default", "png");
# 
# # Create summary box plot of alpha diversity measures using the viridis color palette.
# mbSet<-PlotAlphaBoxData(mbSet, "alpha_diverbox_0","Chao1","Class","default", "png");
# 
# # Calculate alpha-diversity significance using the parametric t-test (two groups).
# mbSet<-PerformAlphaDiversityComp(mbSet, "tt","Class");

## ----eval=FALSE---------------------------------------------------------------
# # Create a PCoA score plot using the viridis custom color palette, data points are colored
# # by their group label.
# mbSet<-PlotBetaDiversity(mbSet, plotNm="beta_diver_0", ordmeth="PCoA", distName="bray", colopt="expfac", metadata="Class",
#                          showlabel="none", taxrank="OTU", taxa="null", alphaopt="Chao1", ellopt="yes", format="png",
#                          dpi=72, custom_col="viridis");
# 
# # Creates a json file for 3D PCoA, use the web to view the interactive 3D PCoA plot
# 	mbSet<-PCoA3D.Anal(mbSet, "PCoA","bray","OTU","expfac","Group","","Chao1","beta_diver3d_0.json")
# 
# # Calculate the beta-diversity significance
# mbSet<-PerformCategoryComp(mbSet, "OTU", "adonis","bray","Group");

## ----eval=FALSE---------------------------------------------------------------
# # Create a core microbiome plot using the viridis color palette
# mbSet<-CoreMicrobeAnalysis(mbSet, "core_micro_0",0.2,0.01,"OTU","bwm","overview", "all_samples", "Group", "CD", "png");

## ----eval=FALSE---------------------------------------------------------------
# # Create a heatmap at the OTU level using the plasma color palette
# mbSet<-PlotHeatmap(mbSet, "heatmap_0","euclidean","ward.D","bwm","Group","OTU","overview","F", "png","T","F");

## ----eval=FALSE---------------------------------------------------------------
# # Create a dendogram of the microbiome data at the OTU level.
# mbSet<-PlotTreeGraph(mbSet, "plot_tree_0","bray","ward.D","Group","OTU", "default", "png");

## ----eval=FALSE---------------------------------------------------------------
# 

## ----eval=FALSE---------------------------------------------------------------
# # Identify and plot the pattern
# mbSet<-Match.Pattern(mbSet, "pearson", "1-2", "Genus", "Group")
# mbSet<-PlotCorr(mbSet, "ptn_1", "png", width=NA)

## ----eval=FALSE---------------------------------------------------------------
# # Classical non-parametric univariate analysis: Mann-Whitney/Kruskall Wallis
# mbSet<-PerformUnivarTest(mbSet, "Group",0.05,"NA","OTU","nonpar");
# 
# # Classical parametric univariate analysis: T-test/ANOVA
# mbSet<-PerformUnivarTest(mbSet, "Group",0.05,"NA","OTU","tt");

## ----eval=FALSE---------------------------------------------------------------
# # edgeR at the OTU level
# mbSet<-PerformRNAseqDE(mbSet, "EdgeR",0.05,"Group","NA","OTU");
# 
# # deseq2 at OTU level
# mbSet<-PerformRNAseqDE(mbSet, "EdgeR",0.05,"Group","NA","OTU");

## ----eval=FALSE---------------------------------------------------------------
# # metegenomeSeq using the zero-inflated Gaussian fit model
# mbSet<-PerformMetagenomeSeqAnal(mbSet, "Group",0.05,"NA","OTU","zigfit");

## ----eval=FALSE---------------------------------------------------------------
# # First perform LefSe analysis
# mbSet<-PerformLefseAnal(mbSet, 0.1, "fdr", 2.0, "Group","F","NA","OTU");
# 
# # Plot LefSe results
# mbSet<-PlotLEfSeSummary(mbSet, 15, "dot", "bar_graph_0","png");

## ----eval=FALSE---------------------------------------------------------------
# # First perform RF analysis
# mbSet<-RF.Anal(mbSet, 500,7,1,"Group","OTU")
# 
# # Plot the RF classification results
# mbSet<-PlotRF.Classify(mbSet, 15, "rf_cls_0","png", width=NA)
# 
# # Plot the RF Variable Importance Plot
# mbSet<-PlotRF.VIP(mbSet, 15, "rf_imp_0","png", width=NA)

