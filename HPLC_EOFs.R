# Sasha Kramer | HPLC Pigment Clustering — R
# Translation of Kramer_cluster.m | Kramer and Siegel, 2019
# https://github.com/sashajane19/HPLC_EOF

# Load required library
library(R.matlab)

# --- Column names ---
pigment_cols <- c(
  "Tchla","Tchlb","Tchlc","ABcaro","ButFuco",
  "HexFuco","Allo","Diadino","Diato","Fuco",
  "Perid","Zea","MVChla","DVchla","Chllide",
  "MVChlb","DVchlb","Chlc12","Chlc3","Lut",
  "Neo","Viola","Phytin","Phide","Pras")

# --- Load data ---
mat_data     <- readMat("Global_HPLC_all.mat")
Global_RHPLC <- as.data.frame(mat_data$Global.RHPLC)
Global_Rlat <- as.vector(mat_data$Global.Rlat)
Global_Rlon <- as.vector(mat_data$Global.Rlon)
colnames(Global_RHPLC) <- pigment_cols

# --- QC ---
detection_limits <- list(
  Phide=0.003, Perid=0.003, MVChlb=0.003, Phytin=0.003,
  ButFuco=0.002, Fuco=0.002, Neo=0.002,
  Pras=0.002, HexFuco=0.002, DVchla=0.002)

for (pigment in names(detection_limits)) {
  limit <- detection_limits[[pigment]]
  mask  <- Global_RHPLC[[pigment]] != 0 & Global_RHPLC[[pigment]] <= limit
  Global_RHPLC[mask, pigment] <- 0
  cat(sprintf("  %8s: %3d value(s) <= %.3f set to 0\n", pigment, sum(mask), limit))
}

# --- Check % below detection ---
n_samples <- nrow(Global_RHPLC)
below_det <- sapply(pigment_cols, function(col)
  100 * sum(Global_RHPLC[[col]] <= 0.001, na.rm=TRUE) / n_samples)
below_det_df <- data.frame(Pigment=names(below_det),
                           Pct=round(below_det,1), row.names=NULL)
print(below_det_df[order(-below_det_df$Pct), ])

# --- Remove degradation products ---
Rpigcluster1 <- Global_RHPLC[, !colnames(Global_RHPLC) %in%
                               c("Chllide","Phytin","Phide")]

# --- Remove redundant / below-detection pigments ---
redundant_pigments <- c("Tchlb","Tchlc","ABcaro","Diadino","Diato",
                        "MVChla","DVchlb","Lut","Pras")
Rpigcluster2 <- Rpigcluster1[, !colnames(Rpigcluster1) %in% redundant_pigments]
label2 <- colnames(Rpigcluster2)

# --- Cluster pigments (absolute) ---
D2 <- as.dist(1 - cor(Rpigcluster2, use="pairwise.complete.obs"))
Z2 <- hclust(D2, method="ward.D2")
plot(Z2, labels=label2, main="Dendrogram — Absolute Pigment Values",
     xlab="", ylab="Linkage Distance", sub="", hang=-1)

# --- Normalize to Tchla ---
normlabel  <- label2[label2 != "Tchla"]
normchl_df <- Rpigcluster2[, normlabel] / Global_RHPLC[["Tchla"]]

# --- Define colors and reorder pigments ---
colors <- rainbow(5)

# Define pigment ratios in order of cluster results
# (These indices are relative to the normchl matrix)
# Original: normchl columns are: ButFuco(1), HexFuco(2), Allo(3), Fuco(4),
#           Perid(5), Zea(6), DVchla(7), MVchlb(8), Chlc12(9), Chlc3(10), Neo(11), Viola(12)
normchlS <- normchl_df[, c(7, 12, 11, 3, 1, 2, 4, 8, 9, 6, 10, 5)]
pigEOF <- normchlS

# --- Mean center pigments before PCA ---
# Calculate means and standard deviations (ignoring NAs)
pigmeans <- colMeans(pigEOF, na.rm = TRUE)
pigstd <- apply(pigEOF, 2, sd, na.rm = TRUE)

# Center and standardize
pigs_center <- sweep(sweep(pigEOF, 2, pigmeans), 2, pigstd, "/")

# --- Calculate EOFs using PCA ---
# prcomp: EOFs_pig = loadings (eigenvectors), AFs_pig = scores (amplitude functions)
# scale=FALSE because we already standardized manually
pca_result <- prcomp(pigs_center, center = FALSE, scale. = FALSE)

EOFs_pig <- pca_result$rotation
AFs_pig <- pca_result$x
eigvalues_pig <- (pca_result$sdev)^2

# --- Calculate variance explained
var_explained <- eigvalues_pig / sum(eigvalues_pig)
var_cutoff <- 0.031  # keep top 6 modes

# Filter modes above cutoff
idx_keep <- which(var_explained > var_cutoff)
expvar <- var_explained[idx_keep]
expvar_percent <- expvar * 100

cat("\nVariance explained by each mode (%):\n")
print(expvar_percent)
cat("\nCumulative variance explained (%):\n")
print(cumsum(expvar_percent))

# Keep top 6 modes (or however many exceed var_cutoff)
EOFs_pigC <- EOFs_pig[, idx_keep]
AFs_pigC <- AFs_pig[, idx_keep]

# Create amplitude function and lat/lon dataframe for later use
latlon <- cbind(Global_Rlat, Global_Rlon)
AFs_pig_latlon <- cbind(latlon, AFs_pigC)

# Extract individual amplitude functions
AF1p <- AFs_pigC[, 1]
AF2p <- AFs_pigC[, 2]
AF3p <- if (ncol(AFs_pigC) >= 3) AFs_pigC[, 3] else NA
AF4p <- if (ncol(AFs_pigC) >= 4) AFs_pigC[, 4] else NA
AF5p <- if (ncol(AFs_pigC) >= 5) AFs_pigC[, 5] else NA
AF6p <- if (ncol(AFs_pigC) >= 6) AFs_pigC[, 6] else NA

# --- Correlate amplitude functions with pigments
# Correlation with lat/lon
corrmat1 <- cbind(AF1p, AF2p, AF3p, AF4p, AF5p, AF6p, pigEOF, latlon)
R_p1 <- cor(corrmat1, use = "pairwise.complete.obs")
# P_p1 <- cor.test would require loop for all pairs; skipped for simplicity

# Correlation of AFs with pigments only
corrmat <- cbind(AF1p, AF2p, AF3p, AF4p, AF5p, AF6p, pigEOF)
R_all <- cor(corrmat, use = "pairwise.complete.obs")

# Trim to AF × pigment correlations
n_modes <- ncol(AFs_pigC)
n_pigs <- ncol(pigEOF)
R_pig_trim <- R_all[1:n_modes, (n_modes + 1):(n_modes + n_pigs)]

cat("\nAF-Pigment correlations (first 3 modes):\n")
print(R_pig_trim[1:min(3, nrow(R_pig_trim)), ])

# --- Prepare plotting labels
# Rearrange based on order of pigments
plotlabel <- c("DVchla","Viola","Neo","Allo","ButFuco","HexFuco","Fuco","MVchlb","Chlc12","Zea","Chlc3","Perid")
plotlabel <- plotlabel[c(10,1,4,2,8,3,5,6,12,11,7,9)]
  
# Rearrange loadings
EOFs_rearr <- EOFs_pigC[c(10,1,4,2,8,3,5,6,12,11,7,9), ]
R_pigR <- t(R_pig_trim)
R_rearr <- R_pigR[c(10,1,4,2,8,3,5,6,12,11,7,9), ]

# --- Plot loadings for Mode 1 (example)
mode_to_plot <- 1

# Create figure with simple bar plot
png("EOF_loadings_mode1.png", width = 1000, height = 600)
par(mar = c(10, 5, 3, 2), cex = 1.3)

# Prepare bar positions and colors (mapping cluster groups)
bar_vals <- EOFs_rearr[1:length(plotlabel), mode_to_plot]
bar_colors <- c(colors[1], colors[1],                               # Zea, DVchla
                colors[3], colors[3], colors[3], colors[3],         # MVchlb, Neo, Viola, Allo
                colors[2], colors[2],                               # ButFuco, HexFuco
                colors[5],                                          # Perid
                colors[4], colors[4], colors[4]                    # Fuco, Chlc12, Chlc3
                )                               

bp <- barplot(bar_vals, names.arg = plotlabel, ylim = c(-1, 1), 
              col = bar_colors, main = paste("EOF Mode", mode_to_plot, "-", 
                                             round(expvar_percent[mode_to_plot], 1), "% variance"),
              ylab = "Loading", las = 2, grid = NA)

abline(h = 0, col = "black", lwd = 1)
grid(NA, NULL, lwd = 1)

# Add correlation values below each bar (if desired, optional)
for (i in 1:length(bar_vals)) {
  text(bp[i], -0.95, paste(round(R_rearr[i, mode_to_plot] * 100, 0)), 
       cex = 1.0, font = 2)
}

dev.off()
cat("\nPlot saved: EOF_loadings_mode1.png\n")

# --- Plot amplitude functions in space 
# For mode 1
# Make sure the 'maps' package is installed and loaded
if (!requireNamespace("maps", quietly = TRUE)) install.packages("maps")
library(maps)

png("AF1_spatial.png", width = 1000, height = 600)
par(mar = c(5, 5, 3, 3))

# Simple scatter plot of AF1 to set up the axes and limits
plot(Global_Rlon, Global_Rlat, col = "white", main = "Amplitude Function Mode 1",
     xlab = "Longitude", ylab = "Latitude", cex = 0.1)

# Color scale
zlim <- c(-2, 2)
col_range <- colorRamp(c("blue", "white", "red"))
n_colors <- 100
colors_scale <- rgb(col_range(seq(0, 1, length.out = n_colors)), maxColorValue = 255)

# Map AF values to colors
AF1_scaled <- pmax(pmin(AF1p, zlim[2]), zlim[1])  # clip to zlim
AF1_norm <- (AF1_scaled - zlim[1]) / (zlim[2] - zlim[1])
AF1_colors <- colors_scale[round(AF1_norm * (n_colors - 1)) + 1]

# Plot points first
points(Global_Rlon, Global_Rlat, col = AF1_colors, pch = 16, cex = 1.5)

# Add global map outline on top of the points
map("world", add = TRUE, col = "black", lwd = 1) 

# Add colorbar (simple text version)
mtext("AF1 magnitude [-2 to 2]", side = 4, line = 1)

dev.off()
cat("Plot saved: AF1_spatial.png\n")