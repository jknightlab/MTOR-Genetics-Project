#------------------------------------------------------------------------------#
# Survival and Logistic Regression
#------------------------------------------------------------------------------#
library(data.table)
library(survival)
library(survminer)
library(ggplot2)
library(ggrepel)

#------------------------------------------------------------------------------#
# File path
#------------------------------------------------------------------------------#
data_file <- "~/MTOR.project/UKB/Sepsis_Any_Any_J15.tab"

#------------------------------------------------------------------------------#
# LOAD DATA
#------------------------------------------------------------------------------#
dat <- fread(data_file)
dat[, days_between := as.numeric(as.character(days_between))]

# Ethnicity filtering
dat2 <- dat[genetic_ethnic_grouping_f22006_0_0 == "Caucasian" & !is.na(`chr1:11246222:G:C_C`), ]

#------------------------------------------------------------------------------#
# Genotype recoding
#------------------------------------------------------------------------------#
dat2$genotype <- ifelse(dat2$`chr1:11246222:G:C_C` %in% "2", 0, 
                        ifelse(dat2$`chr1:11246222:G:C_C` %in% "1", 1,  
                               ifelse(dat2$`chr1:11246222:G:C_C` %in% "0", 2, dat2$`chr1:11246222:G:C_C`)))

#------------------------------------------------------------------------------#
# Survival analysis (Cox PH model)
#------------------------------------------------------------------------------#
# Define event and censoring
dat2[, `:=`(Event_time = pmin(days_between, 28),
            Death_flag_28day = as.integer(days_between <= 28))]

# Fit Cox-PH model
cox <- coxph(Surv(Event_time, Death_flag_28day) ~ genotype + Age_0 + sex_f31_0_0, data = dat2)
summary(cox)

# Plot Forest Plot
ggforest(cox, data = as.data.frame(dat2), main = "28-day Mortality", noDigits = 2)

#------------------------------------------------------------------------------#
# KM curve
#------------------------------------------------------------------------------#
fit <- survfit(Surv(Event_time, Death_flag_28day) ~ genotype, data = dat2)
#
survivalplot <- ggsurvplot(fit,
                           xlab = "Time in days",
                           ylab = "Survival Probability (%)",
                           conf.int = TRUE,
                           surv.scale = "default",
                           pval = signif(summary(cox)$coef[1,5], digits = 4),
                           risk.table = TRUE,
                           break.time.by = 4,
                           xlim = c(0, 28),
                           palette = c("darkorange", "lightblue", "cyan4"),
                           ggtheme = theme_bw(base_size = 13),
                           fontsize = 3.5)

# Adjust y-axis range
survivalplot$plot <- survivalplot$plot +
  coord_cartesian(ylim = c(0.7, 1)) +
  scale_y_continuous(breaks = seq(0.7, 1, by = 0.1), labels = c(70, 80, 90, 100)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Remove elements from risk table
survivalplot$table <- survivalplot$table +
  theme(axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())

# Display the survival plot
survivalplot

#------------------------------------------------------------------------------#
# logistic regression
#------------------------------------------------------------------------------#
survival2 <- dat2[IS == "0", ]

model <- glm(Death_flag_28day ~ genotype + Age_0 + sex_f31_0_0 +
               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7, 
             data = survival2, family = binomial)
summary(model)

# Extract regression results
logistic_results <- data.frame(
  Sample.size = nrow(model.frame(model)),
  beta = coef(summary(model))[2, 1],
  se = coef(summary(model))[2, 2],
  z.value = coef(summary(model))[2, 3],
  CI_beta.low = confint(model, "genotype", level = 0.95)[1],
  CI_beta.hi = confint(model, "genotype", level = 0.95)[2],
  OR = exp(coef(summary(model))[2, 1]),
  P = coef(summary(model))[2, 4]
)

logistic_results



