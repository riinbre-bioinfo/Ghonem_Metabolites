# Ghonem_Metabolites

R scripts and commands used to analyze metabolite data in (manuscript)

Data input files available upon request

Metabolite data was analyzed using linear mixed models and ANOVA as described in (manuscript).  Briefly, Treatment (Healthy, Ursodiol(mono) and Combined (Adjunct)) and Cohort (PSC and PBC) were combined into a single factor with five levels (Healthy, UrsodiolPBC, UrsodiolPSC, CombinedPBC, CombinedPSC) and a linear model was constructed using TreatCohort as the fixed effect and Subject as the random effect to account for samples taken from the same subject.  The model was run in limma with contrasts to compare Healthy vs. Ursodiol, Healthy vs. Combined and Ursodiol vs. Combined (all cohorts combined).  ANOVA was run using Treatment as the fixed effect and Subject as the blocking factor with multiple comparisons made using the Tukey test using the multcomp package.  Heatmaps were generated using the heatmap.2 function of gplots package and all other figures were made using ggplot2 package.  PCA results were generated using prcomp package and plotted using ggbiplot package.
