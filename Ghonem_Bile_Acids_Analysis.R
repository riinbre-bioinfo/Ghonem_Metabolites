# Linear mixed model w/ Treatment and Cohort as a single variable (TreatCohort) w/ five levels and Subject as random effect + contrasts

bileAcids_full.df <- read.csv("bileAcids6.csv", stringsAsFactors = TRUE)
bileAcids_full.df$TreatCohort <- factor(bileAcids_full.df$TreatCohort, levels = c("Healthy", "MonoPBC", "MonoPSC", "AdjunctPBC", "AdjunctPSC"))
bileAcids_full.df$Treatment <- factor(bileAcids_full.df$Treatment, levels = c("Healthy", "Mono","Adjunct"))
bileAcids_meta.df <- bileAcids_full.df[, 1:6]
bileAcids.df <- bileAcids_full.df[-c(1:6)]
bileAcids_log.df <- log(bileAcids.df)
bileAcids_log.df[16, "C4"] <- mean(bileAcids_log.df$C4, na.rm = TRUE) # replace single missing value in C4 column (sample 21) with C4 column mean
bileAcids_design.mat <- model.matrix(~0+TreatCohort, data = bileAcids_meta.df)
colnames(bileAcids_design.mat) <- c("Healthy", "MonoPBC", "MonoPSC", "AdjunctPBC", "AdjunctPSC")
bileAcids_corfit <- duplicateCorrelation(t(bileAcids_log.df), bileAcids_design.mat, block=bileAcids_meta.df$Subject)
bileAcids_contrasts.mat <- makeContrasts(MonoVsHealthy=((MonoPBC+MonoPSC)/2)-Healthy, AdjunctVsHealthy=((AdjunctPBC+AdjunctPSC)/2)-Healthy, PSCVsHealthy=((MonoPSC+AdjunctPSC)/2)-Healthy, PBCVsHealthy=((MonoPBC+AdjunctPBC)/2)-Healthy, MonoPBCVsMonoPSC=MonoPBC-MonoPSC, AdjunctPBCVsAdjunctPSC=AdjunctPBC-AdjunctPSC, AdjunctVsMono=(AdjunctPSC+AdjunctPBC)-(MonoPBC+MonoPSC), PBCVsPSC=(MonoPBC+AdjunctPBC)-(MonoPSC+AdjunctPSC), levels = bileAcids_design.mat)
rownames(bileAcids_contrasts.mat) <- c("Healthy", "MonoPBC", "MonoPSC", "AdjunctPBC", "AdjunctPSC")
bileAcids_SimpleFit <- lmFit(t(bileAcids_log.df), bileAcids_design.mat, block=bileAcids_meta.df$Subject, correlation=bileAcids_corfit$consensus)
bileAcids_contrastsFit <- contrasts.fit(bileAcids_SimpleFit, bileAcids_contrasts.mat)
bileAcids_contrastsFit <- eBayes(bileAcids_contrastsFit)

#PCA plot

bileAcids.pca <- prcomp(bileAcids_log.df, center = TRUE, scale. = TRUE)
ggbiplot(bileAcids.pca, ellipse = TRUE, groups = bileAcids_meta.df$TreatCohort, varname.size = 3, varname.adjust = 3) +
  scale_color_manual(name = "TreatCohort", values = c("darkolivegreen", "magenta", "orange", "blue", "cyan")) +
  scale_shape_manual(name = "TreatCohort", values = c(15, 1, 16, 2, 17)) +
  geom_point(aes(colour=bileAcids_meta.df$TreatCohort, shape = bileAcids_meta.df$TreatCohort), size = 3) +
  theme_bw() +
  theme(legend.direction = "horizontal", legend.position = "top")

# Heat Maps

heatmap.2(t(as.matrix(bileAcids_log.df)), labCol = bileAcids_meta.df$Subject, ColSideColors = c(rep("magenta", 20), rep("blue", 12), rep("orange", 11), rep("cyan", 11), rep("darkolivegreen", 26)), trace = "none", cexCol = 0.75, col = rev(brewer.pal(7, "RdBu")), scale = "row", dendrogram="row", Colv = FALSE)
par(lend = 1)
legend("topright", legend = c("Healthy", "MonoPBC", "AdjunctPBC", "MonoPSC", "AdjunctPSC"), col = c("darkolivegreen", "magenta", "blue", "orange","cyan"), lty = 1, lwd = 10)

# strip chart input

lapply(
  names(bileAcids_full_log.df[c(7:31)]), 
  plotMetabolites, 
  metdata = bileAcids_full_log.df, 
  xvar = "Treatment", 
  shape = "Treatment",
  color = "TreatCohort",
  shapes = c(16, 17,15),
  colors = c("#556A2F", "#FF00FD", "#FFA300", "#0000FC", "#00FDFD"),
  brewer = FALSE,
  grouped = FALSE)

