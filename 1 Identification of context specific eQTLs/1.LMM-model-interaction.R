
#-------------------------------------------------------------------------------
# Get Command-Line Arguments
#-------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript script_name.R <SRS_data_file> <cell_counts_file>")
}

# Assign input files from command-line arguments
srs_data_file <- args[1]
cell_counts_file <- args[2]

cat("Files used in the analysis:\n")
cat("1. SRS data:", srs_data_file, "\n")
cat("2. Cell counts data:", cell_counts_file, "\n\n")


# Load required libraries
library(ggplot2)
library(data.table)
library(lmerTest)
library(interactions)
library(reshape2)

# Load SRS data
SRS <- fread(srs_data_file, stringsAsFactors = FALSE)
SRS[, donor := gsub("_.*", "", V1)]  # Extract donor ID
SRS[, SRS := as.factor(SRS)]         # Convert SRS to factor

# Function to fit LMM and extract summary
fit_lmm <- function(formula, data, reml = FALSE) {
  model <- lmer(formula, data = data, REML = reml)
  return(summary(model))
}

# Function to extract p-value and beta coefficient
extract_stats <- function(model_summary, row = 2) {
  p_value <- signif(model_summary$coefficients[row, 5], digits = 3)
  beta <- signif(model_summary$coefficients[row, 1], digits = 3)
  return(list(p_value = p_value, beta = beta))
}

# Model 1: AdjustedExprn ~ SNP + SRS + (1|donor)
form1 <- AdjustedExprn ~ SNP + SRS + (1 | donor)
lm1_summary <- fit_lmm(form1, SRS)
stats1 <- extract_stats(lm1_summary)
cat("Model 1 - P-value:", stats1$p_value, "Beta:", stats1$beta, "\n")

# Model 2: AdjustedExprn ~ SNP + SRS + SNP*SRS + (1|donor)
form2 <- AdjustedExprn ~ SNP + SRS + SNP * SRS + (1 | donor)
lm2_summary <- fit_lmm(form2, SRS, reml = FALSE)
stats2 <- extract_stats(lm2_summary)
p_interaction <- signif(lm2_summary$coefficients[4, 5], digits = 2)

# Interaction plot for SNP and SRS
plot1 <- interact_plot(lmer(form2, SRS, REML = FALSE), 
                       pred = SNP, modx = SRS, 
                       plot.points = TRUE, jitter = c(0.1, 0), 
                       colors = c("darkorange", "cyan4"), interval = TRUE)
plot1 + theme_bw() + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1, decimal.mark = ',')) +
  annotate("text", x = 1.0, y = 4.1, label = paste0("Pinteraction = ", p_interaction))
ggsave("fig.1d.pdf", width = 3, height = 2.5)

# Load and preprocess cell counts data
Cell.counts <- fread(cell_counts_file, stringsAsFactors = FALSE)
Cell.counts <- Cell.counts[, .(
  `Subject Id`,
  `Day 1 Polymorphonuclear leukocytes`, `Day 3 Polymorphonuclear leukocytes`, `Day 5 Polymorphonuclear leukocytes`,
  `Day 1 Mononuclear leukocytes`, `Day 3 Mononuclear leukocytes`, `Day 5 Mononuclear leukocytes`,
  `Day 1 Lymphocytes`, `Day 3 Lymphocytes`, `Day 5 Lymphocytes`
)]

# Reshape cell counts data using reshape2::melt
aql <- melt(Cell.counts, id.vars = "Subject Id", variable.name = "Day_Cell", value.name = "Count")

# Use colsplit to split Day_Cell into Day, cell_1, and cell_2
aql$Day_Cell = gsub("Day ", "", aql$Day_Cell)
split_cols <- colsplit(aql$Day_Cell, " ", names = c("Day", "cell_1", "cell_2"))

# Combine the split columns with the original data
aql <- cbind(aql, split_cols)

# Create the ID column for merging
aql$ID <- paste(aql$`Subject Id`, aql$Day, sep = "_")

# Combine cell_1 and cell_2 to create Cell.type
aql$Cell.type <- paste(aql$cell_1, aql$cell_2, sep = " ")

# Reshape back to wide format
aql <- dcast(aql, ID ~ Cell.type, value.var = "Count")

# Calculate N.to.L ratio
aql$N.to.L = log2(aql$`Polymorphonuclear leukocytes` / aql$Lymphocytes + 1)

# Merge SRS data with cell counts
data <- merge(SRS, aql, by.x = "V1", by.y = "ID")
data <- data[!is.na(N.to.L) & !is.infinite(N.to.L)]

# Model 3: AdjustedExprn ~ SNP + N.to.L + (1|donor)
form3 <- AdjustedExprn ~ SNP + N.to.L + (1 | donor)
lm3_summary <- fit_lmm(form3, data)
stats3 <- extract_stats(lm3_summary)
cat("Model 3 - P-value:", stats3$p_value, "Beta:", stats3$beta, "\n")

# Model 4: AdjustedExprn ~ SNP + N.to.L + SNP*N.to.L + (1|donor)
form4 <- AdjustedExprn ~ SNP + N.to.L + SNP * N.to.L + (1 | donor)
lm4_summary <- fit_lmm(form4, data, reml = FALSE)
stats4 <- extract_stats(lm4_summary)
p_interaction_n_to_l <- signif(lm4_summary$coefficients[4, 5], digits = 2)

# Interaction plot for SNP and N.to.L
plot2 <- interact_plot(lmer(form4, data, REML = FALSE), 
                       pred = SNP, modx = N.to.L, 
                       plot.points = TRUE, point.alpha = 0.5, jitter = c(0.1, 0), 
                       colors = c("blue", "green", "red"), interval = TRUE)
plot2 + theme_bw() + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1, decimal.mark = ',')) +
  annotate("text", x = 1.0, y = 4.1, label = paste0("Pinteraction = ", p_interaction_n_to_l))
ggsave("fig.1e.pdf", width = 3.4, height = 2.5)