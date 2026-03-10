# ============================================================
# Trophic complexity and ecosystem stability in marine mesocosms
# Combined Analysis Script
# ============================================================
#
# Authors: [Author names]
# Manuscript: [Journal, Year]
#
# Description:
#   This script reproduces all statistical analyses and figures
#   in the manuscript. Data are read from ../data/. Figures are
#   written to ../figures/.
#
# Experiment overview:
#   Fully-crossed mesocosm factorial experiment (n=144 buckets)
#   testing effects of four biotic factors on ecosystem structure
#   and stability: Crab (predator), Amphipods (grazers),
#   Ulva (macroalgae), and Kelp (structure/productivity).
#   Response variables include benthic chlorophyll, oxygen
#   metabolism (GPP, NCP, Respiration), nutrient concentrations,
#   organism biomass, and mesofauna community composition on
#   artificial pot-scrubber substrates. Temporal stability is
#   measured as the variance of detrended time series per bucket.
#
# Script sections:
#   1.  Setup (libraries, paths)
#   2.  Helper functions (temporal variance)
#   3.  Data loading
#   4.  Data wrangling and preparation
#   5.  Biomass calculations
#   6.  Temporal variance calculations
#   7.  Statistical analyses
#       7a. Benthic chlorophyll (BenthoTorch)
#       7b. Algae cover (Ectocarpus)
#       7c. Organism biomass (crab, amphipods, Ulva, kelp, detritus)
#       7d. Nutrients (nitrate, ammonium, phosphate)
#       7e. Oxygen metabolism (NCP, GPP, Respiration)
#       7f. Mesofauna community (pot scrubber)
#   8.  Figures
#       8a. Figure 1: Benthic chlorophyll & metabolism vs. treatments
#       8b. Figure 2: Organism biomass vs. treatments
#       8c. Figure 3: Mesofauna community vs. treatments
#
# To run: set your working directory to the R/ folder, or use
#   the RStudio 'Source' button, before running.
# ============================================================


# ============================================================
# SECTION 1: SETUP
# ============================================================

rm(list = ls())

# --- Libraries ---
library(vegan)
library(car)
library(reshape)
library(psych)
library(lme4)
library(MASS)
library(plyr)
library(lmerTest)
library(mvabund)

# --- Paths (relative to this script's location) ---
data_dir <- file.path("..", "data")
fig_dir  <- file.path("..", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)


# ============================================================
# SECTION 2: HELPER FUNCTIONS
# ============================================================

#' Calculate detrended temporal variance (single response variable)
#'
#' For each mesocosm bucket (1-144), fits a linear regression of
#' log10(Value + 1) ~ Sample (time point) and returns the variance
#' of the residuals. This detrends the time series and quantifies
#' temporal instability independent of any directional trend.
#'
#' @param dd  Data frame with columns: Overall.bucket (integer 1-144),
#'            Sample (time point index), Value (response variable).
#' @return    Data frame with one column Temp_variances (length 144).
calc_temporal_variance <- function(dd) {
  Temp_variances <- numeric(144)
  for (i in 1:144) {
    ind              <- dd$Overall.bucket == i
    model            <- lm(log10(dd$Value[ind] + 1) ~ dd$Sample[ind])
    Temp_variances[i]<- var(resid(model))
  }
  data.frame(Temp_variances = Temp_variances)
}

#' Calculate detrended temporal variance (multivariate: 10 mesofauna species)
#'
#' Same approach as calc_temporal_variance() but applied to each of 10
#' mesofauna taxa simultaneously (columns 9:18 of dd).
#'
#' @param dd  Data frame with columns: bucket, Time, and species abundances
#'            in columns 9:18 (Ostracod, Copepod, Gastropod, Foraminifera,
#'            Isopod, Mite, Brittle.star, Amphipod, Larvae, Nematode).
#' @return    Data frame (144 rows x 10 species) of temporal variances.
calc_temporal_variance_multi <- function(dd) {
  sp.ids <- c("Ostracod", "Copepod", "Gastropod", "Foraminifera", "Isopod",
              "Mite", "Brittle.star", "Amphipod", "Larvae", "Nematode")
  species        <- dd[, 9:18]
  Temp_variances <- matrix(0, 144, 10)
  for (j in 1:10) {
    for (i in 1:144) {
      ind                  <- dd$bucket == i
      model                <- lm(log10(species[ind, j] + 1) ~ dd$Time[ind])
      Temp_variances[i, j] <- var(resid(model))
    }
  }
  colnames(Temp_variances) <- sp.ids
  data.frame(Temp_variances)
}


# ============================================================
# SECTION 3: DATA LOADING
# ============================================================

# --- Experimental design and treatment structure ---
dat    <- read.table(file.path(data_dir, "experiment_design.csv"),
                     header=T, sep=',')

# --- Organism weights and counts ---
crab   <- read.table(file.path(data_dir, "crab_weights.csv"),
                     header=T, sep=',')
pod    <- read.table(file.path(data_dir, "amphipod_counts.csv"),
                     header=T, sep=',')
dry    <- read.table(file.path(data_dir, "dry_weights.csv"),
                     header=T, sep=',')

# --- Benthic chlorophyll (BenthoTorch fluorometry, 4 sampling periods) ---
# Sample 1 is embedded in experiment_design.csv (columns 16-19, 24)
ben2   <- read.table(file.path(data_dir, "benthic_chl_sample2.csv"),
                     header=T, sep=',')
ben3   <- read.table(file.path(data_dir, "benthic_chl_sample3.csv"),
                     header=T, sep=',')
ben4   <- read.table(file.path(data_dir, "benthic_chl_sample4.csv"),
                     header=T, sep=',')

# --- Algae ---
algae  <- read.table(file.path(data_dir, "algae_cover.csv"),
                     header=T, sep=',')

# --- Nutrients ---
nut1   <- read.table(file.path(data_dir, "nutrients.csv"),
                     header=T, sep=',')
nutfin <- read.table(file.path(data_dir, "nitrate_initial_final.csv"),
                     header=T, sep=',')

# --- Oxygen metabolism (raw incubation data) ---
oxy    <- read.table(file.path(data_dir, "oxygen_metabolism.csv"),
                     header=T, sep=',')

# --- Mesofauna community (pot scrubber substrates) ---
ps         <- read.table(file.path(data_dir, "mesofauna_community_counts.csv"),
                         header=T, sep=',')
pslengths  <- read.table(file.path(data_dir, "mesofauna_scrubber_lengths.csv"),
                         header=T, sep=',')
psfaunasize<- read.table(file.path(data_dir, "mesofauna_size_measurements.csv"),
                         header=T, sep=',')
psfaunaconv<- read.table(file.path(data_dir, "length_weight_conversions.csv"),
                         header=T, sep=',')


# ============================================================
# SECTION 4: DATA WRANGLING AND PREPARATION
# ============================================================

# --- Lookup table: bucket positions within mesocosm tables ---
# 8 tables x 3 columns x 1 row = 24 positions; replicated across 6 tables = 144 buckets
tabloc          <- matrix(0, 24, 3)
tabloc[, 1]     <- seq(1, 24)
tabloc[, 2]     <- rep(1:8, each=3)
tabloc[, 3]     <- rep(1:3, times=8)
tabloc          <- data.frame(tabloc)
names(tabloc)   <- c("Table.bucket", "Row", "Column")

# --- Main design data: extract treatments and flow rate ---
dat1            <- dat[, 1:7]
dat1$flowpermin <- 2 * rowMeans(dat[, 12:14])   # mean of 3 flow meter readings x2
dat1$ecto1sam   <- dat$Ecto.est.1st.sam          # Ectocarpus at first sampling
dat1$Kelp       <- factor(dat1$Kelp, levels=c("No", "One piece", "Ground up"))
dat1            <- join(dat1, tabloc, type="left", by="Table.bucket")
Flow            <- dat1[, c(3, 8)]               # used to add flow to other datasets
dat1$Kelp       <- factor(dat1$Kelp, levels=c("No", "One piece", "Ground up"))

# --- Oxygen metabolism: fix data entry errors, reshape, calculate rates ---
oxy[306, 17] <- 15.13       # corrected DO value
oxy[567, 10] <- "12:26"     # corrected time
oxy[567, 13] <- 1.92        # corrected value
# Reshape: two incubation periods (light/dark) into long format
oxy1       <- oxy[, c(1:8, 10, 11, 9, 12, 14, 13)]; oxy1$inc_num <- 1
oxy2       <- oxy[, c(1:8, 12, 14:18)];              oxy2$inc_num <- 2
oxy2       <- setNames(oxy2, names(oxy1))
oxy3       <- rbind(oxy1, oxy2)
# Carbon flux per hour (mgC h-1) from DO change
oxy3$Cperh <- ((oxy3$do2 - oxy3$do1) * 10 * 0.375) / oxy3$inc1hour
# Calculate GPP = NCP_light - NCP_dark (assumes Resp constant in light/dark)
oxyl       <- oxy3[oxy3$incub == "light", ]
oxyd       <- oxy3[oxy3$incub == "dark",  ]
oxysimp    <- merge(oxyl, oxyd, by=c("Sample", "Table", "Overall.bucket"))
oxysimp$GPP<- oxysimp$Cperh.x - oxysimp$Cperh.y
oxysimp    <- oxysimp[, c(1:17, 22:30)]
names(oxysimp)[c(16, 25, 5, 6, 7, 8)] <- c("NCP", "Resp", "Ulva", "Pods", "Crab", "Kelp")
oxysimp$Kelp <- factor(oxysimp$Kelp, levels=c("No", "One piece", "Ground up"))
oxysimp      <- join(oxysimp, dat1, type="left", by="Overall.bucket")

# --- Benthic chlorophyll: combine 4 sampling periods into one data frame ---
ben        <- dat[, c(1:7, 16:19, 24)]; ben$sample <- 1
ben2       <- cbind(dat[, 1:7], ben2[, 4:8]);       ben2$sample <- 2
ben3       <- cbind(dat[, 1:7], ben3[, 3:7]);       ben3$sample <- 3
ben4       <- cbind(dat[, 1:7], ben4[, 3:7]);       ben4$sample <- 4
benth      <- rbind(ben, ben2, ben3, ben4)
benth$Cyano       <- as.numeric(benth$Cyano)
benth$Green.Algae <- as.numeric(benth$Green.Algae)
benth$Diatoms     <- as.numeric(benth$Diatoms)
benth$Total.Conc. <- as.numeric(benth$Total.Conc.)
benth      <- join(benth, tabloc, type="left", by="Table.bucket")
benth      <- join(benth, Flow,   type="left", by="Overall.bucket")

# --- Dry weights: split by organism type ---
names(dry)[2] <- "Overall.bucket"
# Ulva
ud       <- dry[dry$Type == "ulva", ]
ud       <- ud[order(ud$Overall.bucket), ]; rownames(ud) <- NULL
ulvadry  <- join(dat1, ud, type="left", by="Overall.bucket")
ulvadry  <- ulvadry[ulvadry$Ulva == "Yes", ]
ulvadry[ulvadry$Overall.bucket == 6, 17] <- 0   # dead Ulva, set to zero
ulvadry[ulvadry$Overall.bucket == 9, 17] <- 0
# Amphipods
ad       <- dry[dry$Type == "pod", ]
ad       <- ad[order(ad$Overall.bucket), ]; rownames(ad) <- NULL
amphdry  <- join(dat1, ad, type="left", by="Overall.bucket")
amphdry  <- amphdry[amphdry$Pods == "Yes", ]
# Kelp
kp       <- dry[dry$Type == "kelp", ]
kp       <- kp[order(kp$Overall.bucket), ]; rownames(kp) <- NULL
kelpdry  <- join(dat1, kp, type="left", by="Overall.bucket")
kelpdry  <- kelpdry[kelpdry$Kelp != "No", ]
kelpdry$biomass[kelpdry$Overall.bucket %in% c(33, 90, 113)] <- 0  # missing, set to zero
kelpdry  <- kelpdry[-71, ]
# Detritus (0.5 mm sieve fraction)
det      <- dry[dry$Type == "mm0.5", ]
det      <- det[order(det$Overall.bucket), ]; rownames(det) <- NULL
detdry   <- join(dat1, det, type="left", by="Overall.bucket")
# Crab dry weight
cdry     <- dry[dry$Type == "crab", ]
cdry     <- cdry[order(cdry$Overall.bucket), ]; rownames(cdry) <- NULL
cd1      <- cdry[, c(2, 7)]; names(cd1)[2] <- "biomass_dry"

# --- Nitrate (initial vs. final) ---
names(nutfin)[1] <- "Overall.bucket"
nutfin   <- join(dat1, nutfin, type="left", by="Overall.bucket")


# ============================================================
# SECTION 5: BIOMASS CALCULATIONS
# ============================================================

# --- Crab biomass change (wet weight initial to final) ---
crabstay          <- crab[is.na(crab$weight2), ]        # keep crabs not replaced
crabstay          <- crabstay[-c(3, 16, 19, 47), ]       # remove incomplete records
c1                <- merge(crabstay, dat1, by.x="bucket",
                           by.y="Overall.bucket", all.x=T, sort=F)
c1$changebiomass  <- c1$weightfin - c1$weight
c1$pbio           <- c1$weightfin - c1$weight / c1$weight  # proportional change
names(c1)[1]      <- "Overall.bucket"

# --- Dry:wet conversion factors ---
# Crab: mean dry:wet ratio = 0.22 (measured subset)
crabconv       <- join(c1, cd1, type="left", by="Overall.bucket")
crabconv$drywet<- crabconv$biomass_dry / crabconv$weightfin
names(crab)[1] <- "Overall.bucket"
crab$weightfin_dry <- crab$weightfin * 0.22

# Amphipods: mean dry:wet ratio = 0.21 (measured subset)
names(pod)[1]  <- "Overall.bucket"
names(pod)[3]  <- "wet_biomass"
amphconv       <- join(pod[, c(1, 3)], amphdry[, c(3, 17)],
                       type="left", by="Overall.bucket")
amphconv$per_conv <- amphconv$biomass / amphconv$wet_biomass

# Amphipod data with experimental design
pod1   <- merge(pod, dat1, by="Overall.bucket", all.y=T, sort=F)
podi   <- pod1[pod1$Pods == "Yes", ]
podi[62:72, 2:3] <- 0   # zero counts for buckets with no individuals at end

# --- Final biomass table: all organism groups ---
wt           <- join(dat1[, 1:9],       crab[, c(1, 17)],    type="left", by="Overall.bucket")
wt           <- join(wt,                cd1,                  type="left", by="Overall.bucket")
names(wt)[10:11] <- c("Crab_initial", "Crab_final")
wt           <- join(wt,                amphdry[, c(3, 17)],  type="left", by="Overall.bucket")
names(wt)[12]<- "Amphipod_final"
wt           <- join(wt,                ulvadry[, c(3, 17)],  type="left", by="Overall.bucket")
names(wt)[13]<- "Ulva_final"
wt           <- join(wt,                kelpdry[, c(3, 17)],  type="left", by="Overall.bucket")
names(wt)[14]<- "Kelp_final"
write.table(wt, file.path(data_dir, "final_biomass.csv"), row.names=F, sep=",")

# --- Trophic interaction strength: organism-level (log ratio / biomass / time) ---
w      <- wt[, c(6, 5, 4, 7, 11:14)]
levels(w$Kelp)[levels(w$Kelp) %in% c("One piece", "Ground up")] <- "Yes"
t_exp  <- 35   # experiment duration (days)
tab    <- as.data.frame(matrix(0, 4, 4))
names(tab)     <- c("Crab", "Pods", "Ulva", "Kelp")
row.names(tab) <- c("Crab", "Pods", "Ulva", "Kelp")
for (i in 1:4) {
  for (j in 5:8) {
    y          <- mean(w[w[, i] == "Yes", j], na.rm=T)
    n          <- mean(w[w[, i] == "No",  j], na.rm=T)
    a          <- mean(w[, i + 4], na.rm=T)
    tab[j-4, i]<- log(y / n) / (a * t_exp)
  }
}
tab[1, 1] <- NaN   # no crab-free control for crab biomass response

# --- Mesofauna biomass: calculate from individual length measurements ---
lc             <- 1.19   # eye lens conversion factor: 1 unit = 1.19 mm
names(psfaunaconv)[1] <- "Species"
ps1            <- melt(ps, id=c(1:8))
names(ps1)[9:10]<- c("Species", "Abundance")
# Standardise species names across data frames
levels(psfaunasize$Species)[levels(psfaunasize$Species)=="Amphipods"]    <- "Amphipod"
levels(psfaunasize$Species)[levels(psfaunasize$Species)=="Isopods"]      <- "Isopod"
levels(psfaunasize$Species)[levels(psfaunasize$Species)=="Mites"]        <- "Mite"
levels(psfaunasize$Species)[levels(psfaunasize$Species)=="Brittle stars"]<- "Brittle.star"
levels(psfaunasize$Species)[levels(psfaunasize$Species)=="Nematodes"]    <- "Nematode"
ps2    <- join(ps1, psfaunasize, type="left")
ps3    <- cbind(ps2[, 1:10], ps2[, 11:43] * lc)    # convert pixel units to mm
ps3.1  <- join(ps3, psfaunaconv[, 1:3], type="left")
ps3.2  <- (ps3.1$m * ps3.1[, 11:43]) + ps3.1$b     # apply L-W regression
ps3.3  <- cbind(ps3.1[1:10], ps3.2)
ps3.3$mean.biomass       <- rowMeans(ps3.3[11:43], na.rm=T)
ps3.3$biomass            <- ps3.3$Abundance * ps3.3$mean.biomass
ps3.4  <- join(ps3.3, pslengths, type="left")
ps3.4$biomassperpslength <- ps3.4$biomass / ps3.4$Scrubber.length
psf    <- ps3.4[, -c(31:43)]    # remove intermediate length columns
write.table(psf, file.path(data_dir, "mesofauna_biomass_calculated.csv"),
            row.names=F, sep=",")

# --- Trophic interaction strength: mesofauna community ---
ps4    <- ps3.4[ps3.4$Time == 4, c(3, 7, 8, 6, 5, 9, 47)]
names(ps4)[2:5] <- c("Crab", "Pods", "Ulva", "Kelp")
ps5    <- cast(ps4, Crab+Pods+Ulva+Kelp+bucket ~ Species,
               mean, value='biomassperpslength')
ps6    <- cbind(ps5[, 1:4], wt[, 11:14], ps5[, c(13,12,7,9,8,10,14,11,15,6)])
levels(ps6$Kelp)[levels(ps6$Kelp) %in% c("Solid", "Ground")] <- "Yes"
tabf   <- as.data.frame(matrix(0, 10, 4))
names(tabf)     <- c("Crab", "Pods", "Ulva", "Kelp")
row.names(tabf) <- c("Amphipod","Brittle.star","Copepod","Foraminifera",
                     "Gastropod","Isopod","Larvae","Mite","Nematode","Ostracod")
for (i in 1:4) {
  for (j in 9:18) {
    y           <- mean(ps6[ps6[, i] == "Yes", j], na.rm=T)
    n           <- mean(ps6[ps6[, i] == "No",  j], na.rm=T)
    a           <- mean(ps6[, i + 4], na.rm=T)
    tabf[j-8, i]<- log(y / n) / (a * t_exp)
  }
}
write.table(ps6,  file.path(data_dir, "community_biomass_interactions.csv"),
            row.names=F, sep=",")
write.table(tabf, file.path(data_dir, "community_interactions.csv"),
            row.names=F, sep=",")


# ============================================================
# SECTION 6: TEMPORAL VARIANCE CALCULATIONS
# ============================================================

# --- Benthic chlorophyll ---
dd_benth         <- benth
dd_benth$Sample  <- dd_benth$sample
dd_benth$Value   <- dd_benth$Total.Conc.
benthvar         <- cbind(dat1, calc_temporal_variance(dd_benth))
benthvar$Kelp    <- factor(benthvar$Kelp, levels=c("No","One piece","Ground up"))
names(benthvar)[12] <- "bentvar"

# --- Oxygen metabolism (NCP and GPP) ---
dd_oxy        <- oxysimp
dd_oxy$Value  <- dd_oxy$NCP
ncpvar        <- cbind(dat1, calc_temporal_variance(dd_oxy))
names(ncpvar)[12] <- "ncpvar"

dd_oxy$Value  <- dd_oxy$GPP
gppvar        <- cbind(dat1, calc_temporal_variance(dd_oxy))
names(gppvar)[12] <- "gppvar"

# --- Mesofauna community on pot scrubbers ---
# Merge community counts with experimental design
ps_main        <- merge(ps, dat1, by.x="bucket", by.y="Overall.bucket",
                        all.x=T, sort=F)
ps_main2       <- ps_main[, -c(19:24)]
ps_main2$total_abund <- rowSums(ps_main2[9:18])

# Aggregate biomass by bucket and time point
psbio2    <- aggregate(cbind(biomass, Scrubber.length) ~ bucket + Time,
                       FUN=sum, data=psf)
psbio3    <- merge(ps_main2[, c(1:7, 19:22)], psbio2, all.x=T, sort=F)
ps_full   <- merge(ps_main2, psbio2, all.x=T, sort=F)
ps_full$rich      <- rowSums(ps_full[, 9:18] > 0, na.rm=TRUE)
ps_full$diversity <- apply(ps_full[, 9:18], 1,
                           function(x) diversity(x, "simpson"))

# Temporal variance: abundance
dd_ps                <- na.omit(ps_full)
dd_ps$Sample         <- dd_ps$Time
dd_ps$Value          <- dd_ps$total_abund
dd_ps$Overall.bucket <- dd_ps$bucket
psvar     <- cbind(dat1, calc_temporal_variance(dd_ps))
psvar$Kelp<- factor(psvar$Kelp, levels=c("No","One piece","Ground up"))

# Temporal variance: biomass
dd_bio                <- na.omit(psbio3)
dd_bio$Sample         <- dd_bio$Time
dd_bio$Value          <- dd_bio$biomass
dd_bio$Overall.bucket <- dd_bio$bucket
psbiovar  <- cbind(dat1, calc_temporal_variance(dd_bio))

# Temporal variance: richness
dd_rich                <- na.omit(ps_full)
dd_rich$Sample         <- dd_rich$Time
dd_rich$Value          <- dd_rich$rich
dd_rich$Overall.bucket <- dd_rich$bucket
psrichvar <- cbind(dat1, calc_temporal_variance(dd_rich))

# Temporal variance: diversity (Simpson)
dd_div                <- na.omit(ps_full)
dd_div$Sample         <- dd_div$Time
dd_div$Value          <- dd_div$diversity
dd_div$Overall.bucket <- dd_div$bucket
psdivvar  <- cbind(dat1, calc_temporal_variance(dd_div))

# Combined variance table (used for Figure 3)
varg          <- cbind(psvar, psbiovar[, 12], psrichvar[, 12], psdivvar[, 12])
colnames(varg)[12:15] <- c("abund", "biomass", "rich", "div")
write.table(varg, file.path(data_dir, "mesofauna_variance.csv"),
            col.names=T, sep=",")


# ============================================================
# SECTION 7: STATISTICAL ANALYSES
# ============================================================

# -------------------------------------------------------
# 7a. Benthic chlorophyll (BenthoTorch: Total Concentration)
# -------------------------------------------------------

benth_clean <- na.omit(benth)

# -- All sampling periods: find Box-Cox transformation ---
mod_benth <- lm(Total.Conc. + 1 ~ Crab*Pods*Ulva*Kelp, data=benth_clean)
bc_benth  <- boxcox(mod_benth, plotit=FALSE)
mm_benth  <- bc_benth$x[which.max(bc_benth$y)]

# Best model: main effects + flow + random slope for sample period over table
m_benth <- lmer((Total.Conc. + 1)^mm_benth ~ Crab + Pods + Ulva + Kelp +
                  flowpermin + Crab:flowpermin + (sample | Table),
                dat=benth_clean)
summary(m_benth)
Anova(m_benth)

# -- Final sampling period only ---
b4      <- benth_clean[benth_clean$sample == 4, ]
mod_b4  <- lm(Total.Conc. ~ Crab*Pods*Ulva*Kelp, data=b4)
bc_b4   <- boxcox(mod_b4, plotit=FALSE)
mm_b4   <- bc_b4$x[which.max(bc_b4$y)]
m_benth4<- lmer((Total.Conc.)^mm_b4 ~ Crab + Pods + Ulva + Kelp +
                  flowpermin + (1 | Table),
                dat=b4)
summary(m_benth4)
Anova(m_benth4)

# -- Benthic chlorophyll temporal variance ---
mod_bv  <- lm(bentvar ~ Crab*Pods*Ulva*Kelp, data=benthvar)
bc_bv   <- boxcox(mod_bv, plotit=FALSE)
mm_bv   <- bc_bv$x[which.max(bc_bv$y)]
m_bv    <- lmer((bentvar)^mm_bv ~ Crab + Pods + Ulva + Kelp + flowpermin +
                  (1 | Table),
                dat=benthvar)
summary(m_bv)
Anova(m_bv)


# -------------------------------------------------------
# 7b. Algae cover (Ectocarpus)
# -------------------------------------------------------

alg1      <- join(algae, Flow, type="left", by="Overall.bucket")
alg1$resp <- alg1$Ectocarpus.cover
mod_alg   <- lm(resp + 1 ~ Crab*Pods*Ulva*Kelp, data=alg1)
bc_alg    <- boxcox(mod_alg, plotit=FALSE)
mm_alg    <- bc_alg$x[which.max(bc_alg$y)]
m_alg     <- lmer(sqrt(resp) ~ Crab + Pods + Ulva + Kelp + flowpermin +
                    (1 | Table),
                  dat=alg1)
summary(m_alg)
Anova(m_alg)


# -------------------------------------------------------
# 7c. Organism biomass
# -------------------------------------------------------

# -- Crab proportional biomass change ---
c1$resp <- c1$pbio
mod_c   <- lm(resp + 1 ~ Pods*Ulva*Kelp, data=c1)
bc_c    <- boxcox(mod_c, plotit=FALSE)
mm_c    <- bc_c$x[which.max(bc_c$y)]
m_crab  <- lmer((resp)^mm_c ~ Pods + Ulva + Kelp + flowpermin +
                  (1 | Table),
                REML=F, dat=c1)
summary(m_crab)
Anova(m_crab)

# -- Ulva dry biomass ---
ulvadry$resp <- ulvadry$biomass
mod_u        <- lm(resp + 1 ~ Crab*Pods*Kelp, data=ulvadry)
bc_u         <- boxcox(mod_u, plotit=FALSE)
mm_u         <- bc_u$x[which.max(bc_u$y)]
m_ulva       <- lmer((resp)^mm_u ~ Crab + Pods + Kelp + flowpermin +
                       (1 | Table),
                     REML=F, dat=ulvadry)
summary(m_ulva)
Anova(m_ulva)

# -- Amphipod dry biomass ---
amphdry$resp <- amphdry$biomass
mod_a        <- lm(resp ~ Crab*Ulva*Kelp, data=amphdry)
bc_a         <- boxcox(mod_a, plotit=FALSE)
mm_a         <- bc_a$x[which.max(bc_a$y)]
m_amph       <- lmer((resp)^mm_a ~ Crab + Ulva + Kelp + flowpermin +
                       (1 | Table),
                     REML=F, dat=amphdry)
summary(m_amph)
Anova(m_amph)

# -- Kelp dry biomass ---
kelpdry$resp <- kelpdry$biomass + 1
mod_k        <- lm(resp ~ Crab*Pods*Ulva, data=kelpdry)
bc_k         <- boxcox(mod_k, plotit=FALSE)
mm_k         <- bc_k$x[which.max(bc_k$y)]
m_kelp       <- lmer((resp)^mm_k ~ Crab + Pods + Ulva + flowpermin +
                       (1 | Table),
                     dat=kelpdry)
summary(m_kelp)
Anova(m_kelp)

# -- Detritus (0.5 mm fraction) dry biomass ---
detdry$resp  <- detdry$biomass + 1
mod_d        <- lm(resp ~ Crab*Pods*Ulva*Kelp, data=detdry)
bc_d         <- boxcox(mod_d, plotit=FALSE)
mm_d         <- bc_d$x[which.max(bc_d$y)]
m_det        <- lmer((resp)^mm_d ~ Crab + Pods + Ulva + Kelp + flowpermin +
                       (1 | Table),
                     dat=detdry)
summary(m_det)
Anova(m_det)


# -------------------------------------------------------
# 7d. Nutrients
# -------------------------------------------------------

nut1     <- nut1[1:144, 1:11]
nut1$Kelp<- factor(nut1$Kelp, levels=c("No","One piece","Ground up"))
nut1$Crab<- factor(nut1$Crab)
nut1$Pods<- factor(nut1$Pods)
nut1$Ulva<- factor(nut1$Ulva)

# Nitrate (log-normal)
m_no3  <- lm(log(Nitrate) ~ Crab + Pods + Ulva + Kelp, dat=nut1)
summary(m_no3); anova(m_no3)

# Ammonium
m_nh4  <- lm(log(Ammonium + 1) ~ Crab + Pods + Ulva + Kelp, dat=nut1)
summary(m_nh4); anova(m_nh4)

# Phosphate (log-normal)
m_po4  <- lm(log(Phosphate) ~ Crab + Pods + Ulva + Kelp, dat=nut1)
summary(m_po4); anova(m_po4)

# Nitrate final concentration
nutfin$resp <- nutfin$Final
mod_nf      <- lm(resp ~ Crab*Pods*Ulva*Kelp, data=nutfin)
bc_nf       <- boxcox(mod_nf, plotit=FALSE)
mm_nf       <- bc_nf$x[which.max(bc_nf$y)]
m_nitfin    <- lmer((resp)^mm_nf ~ Crab + Pods + Ulva + Kelp + flowpermin +
                      (1 | Table),
                    dat=nutfin)
summary(m_nitfin)
Anova(m_nitfin)


# -------------------------------------------------------
# 7e. Oxygen metabolism (NCP, GPP, Respiration)
# -------------------------------------------------------

# Final sampling period
oo      <- oxysimp[oxysimp$Sample == 4, ]

# -- Net Community Production (NCP) ---
oo$resp <- oo$NCP
m_ncp   <- lmer(resp ~ Crab + Pods + Ulva + Kelp + flowpermin +
                  (1 | Table),
                dat=oo, REML=F)
summary(m_ncp)
Anova(m_ncp)

# -- Gross Primary Production (GPP) ---
oo$resp <- oo$GPP
m_gpp   <- lmer(resp ~ Crab + Pods + Ulva + Kelp + flowpermin +
                  (1 | Table / Column),
                dat=oo, REML=F)
summary(m_gpp)
Anova(m_gpp)

# -- Respiration ---
oo$resp <- oo$Resp
m_resp  <- lmer(resp ~ Crab + Pods + Ulva + Kelp + flowpermin +
                  (1 | Table),
                dat=oo, REML=F)
summary(m_resp)
Anova(m_resp)

# -- NCP temporal variance ---
oxvar_ncp      <- ncpvar
oxvar_ncp$resp <- oxvar_ncp$ncpvar
mod_ov         <- lm(resp ~ Crab*Pods*Ulva*Kelp, data=oxvar_ncp)
bc_ov          <- boxcox(mod_ov, plotit=FALSE)
mm_ov          <- bc_ov$x[which.max(bc_ov$y)]
m_ncpvar       <- lmer((resp)^mm_ov ~ Crab + Pods + Ulva + Kelp + flowpermin +
                         (1 | Table),
                       dat=oxvar_ncp)
summary(m_ncpvar)
Anova(m_ncpvar)

# -- GPP temporal variance ---
oxvar_gpp      <- gppvar
oxvar_gpp$resp <- oxvar_gpp$gppvar
mod_gv         <- lm(resp ~ Crab*Pods*Ulva*Kelp, data=oxvar_gpp)
bc_gv          <- boxcox(mod_gv, plotit=FALSE)
mm_gv          <- bc_gv$x[which.max(bc_gv$y)]
m_gppvar       <- lmer((resp)^mm_gv ~ Crab + Pods + Ulva + Kelp + flowpermin +
                         (1 | Table),
                       dat=oxvar_gpp)
summary(m_gppvar)
Anova(m_gppvar)


# -------------------------------------------------------
# 7f. Mesofauna community (pot scrubber)
# -------------------------------------------------------

# Subset to final time point for cross-sectional analyses
t4   <- na.omit(ps_full[ps_full$Time == 4, ])
cc4  <- t4[, 9:18]
var4 <- t4[, c(4, 20:26)]

# -- Multivariate community composition: PERMANOVA (adonis) ---
ad4  <- adonis(log(cc4 + 1) ~ Crab + Pods + Ulva + Kelp + flowpermin,
               Strata=t4$Table.x, var4)
print(ad4)

# -- Multivariate community composition: RDA ---
mod_rda <- rda(log(cc4 + 1) ~ Crab + Pods + Ulva + Kelp + flowpermin, var4)
vif.cca(mod_rda)
summary(mod_rda)

# -- Univariate: total abundance at final time point ---
rr      <- ps_full[ps_full$Time == 4, ]
rr$resp <- rr$total_abund
mod_ps  <- lm(resp ~ Crab*Pods*Ulva*Kelp, data=rr)
bc_ps   <- boxcox(mod_ps, plotit=FALSE)
mm_ps   <- bc_ps$x[which.max(bc_ps$y)]
m_psab  <- lmer((resp)^mm_ps ~ Crab + Pods + Ulva + Kelp + flowpermin +
                  (1 | Table.x),
                REML=F, dat=rr)
summary(m_psab)
Anova(m_psab)

# -- Univariate: richness at final time point ---
rr$resp  <- rr$rich
m_psrich <- lmer((resp)^mm_ps ~ Crab + Pods + Ulva + Kelp + flowpermin +
                   (1 | Table.x),
                 REML=F, dat=rr)
summary(m_psrich)
Anova(m_psrich)

# -- Univariate: diversity at final time point ---
rr$resp  <- rr$diversity
m_psdiv  <- lmer((resp)^mm_ps ~ Crab + Pods + Ulva + Kelp + flowpermin +
                   (1 | Table.x),
                 REML=F, dat=rr)
summary(m_psdiv)
Anova(m_psdiv)

# -- Temporal variance: abundance, richness, diversity ---
rr_var      <- psdivvar
rr_var$resp <- rr_var$Temp_variances
mod_pv      <- lm(resp ~ Crab*Pods*Ulva*Kelp, data=rr_var)
bc_pv       <- boxcox(mod_pv, plotit=FALSE)
mm_pv       <- bc_pv$x[which.max(bc_pv$y)]
m_psdivvar  <- lmer((resp)^mm_pv ~ Crab + Pods + Ulva + Kelp + flowpermin +
                      (1 | Table),
                    REML=F, dat=rr_var)
summary(m_psdivvar)
Anova(m_psdivvar)

# -- Multivariate GLM (mvabund): species-level effects ---
mv_fom   <- mvabund::mvabund(cc4)
mod_mv   <- mvabund::manyglm(mv_fom ~ (Crab + Pods + Ulva + Kelp + flowpermin),
                              data=var4, family="negative.binomial",
                              cor.type="I", composition=FALSE)
mod_sum  <- anova(mod_mv)
mod_full <- anova(mod_mv, p.uni="adjusted")
print(mod_sum)


# ============================================================
# SECTION 8: FIGURES
# ============================================================

# -------------------------------------------------------
# 8a. Figure 1: Benthic chlorophyll and metabolism vs. treatments
# -------------------------------------------------------
# Combines BenthoTorch Chl, Chl variance, NCP, NCP variance,
# GPP, GPP variance as rows; Crab, Amphipods, Ulva, Kelp, Water
# flow as columns. Open circles = means; * = p < 0.05.

# Assemble plotting data frame
benth4   <- benth[benth$sample == 4, ]
oxy4     <- oxysimp[oxysimp$Sample == 4, ]
all_fig1 <- cbind(dat1,
                  nut1[1:144, c(8, 11, 10)],
                  benth4$Total.Conc.,
                  benthvar$bentvar,
                  oxy4$NCP,
                  ncpvar$ncpvar,
                  oxy4$GPP,
                  gppvar$gppvar)

pdf(file.path(fig_dir, "Figure1_chl_metabolism.pdf"), width=10, height=12)
par(mfrow=c(6, 5), mar=c(.3, .5, .2, .5), oma=c(3, 5.5, .5, 0))
for (i in 15:20) {
  for (j in c(6, 5, 4, 7, 8)) {
    y    <- c(0, max(all_fig1[, i], na.rm=T))
    xll  <- c(.5, 2.5); space <- c(1:2); cl <- c("white", "gray45")
    xl   <- c("No", "Yes")
    if (j == 7) { xl <- c("No","Yes","Grnd"); xll <- c(.5,3.5)
                  space <- c(1:3); cl <- c("white","gray80","gray45") }
    ysig <- max(y) * 0.95
    if (i == 16) { y <- c(0, 0.032); ysig <- 0.027 }
    if (i == 17) { y <- c(0, 20) }
    if (i == 20) { y <- c(0, 0.051); ysig <- 0.045 }
    if (j != 8) {
      resp <- all_fig1[, i]; fac1 <- all_fig1[, j]
      boxplot(resp ~ fac1, outline=F, col=cl, axes=F, at=space,
              xlim=xll, ylim=y, bty='l')
      if (j == 6) axis(2, las=2) else axis(2, labels=F)
      if (i == 20) { axis(1, las=1, at=space, labels=xl, tck=0, padj=-1)
                     abline(-0.002, 0) } else axis(1, labels=F, tck=0)
      q <- tapply(resp, list(fac1), mean, na.rm=T)
      points(space, as.vector(q), pch=1, col='black', cex=1.2)
      # Significance indicators (p < 0.05 from Section 7)
      if (j==5 & i==15) points(1.5, ysig, pch=8, col='black', cex=1.2) # amphipod chl
      if (j==6 & i==16) points(1.5, ysig, pch=8, col='black', cex=1.2) # crab chl-var
      if (j==6 & i==17) points(1.5, ysig, pch=8, col='black', cex=1.2) # crab NCP
      if (j==4 & i==17) points(1.5, ysig, pch=8, col='black', cex=1.2) # ulva NCP
      if (j==6 & i==19) points(1.5, ysig, pch=8, col='black', cex=1.2) # crab GPP
      if (j==4 & i==19) points(1.5, ysig, pch=8, col='black', cex=1.2) # ulva GPP
      if (j==4 & i==20) points(1.5, ysig, pch=8, col='black', cex=1.2) # ulva GPP-var
    } else {
      plot(all_fig1[,j], all_fig1[,i], pch=1, cex=.6, ylab="", xlab="",
           xaxt="n", yaxt="n", bty='l', xlim=c(150,900), ylim=y)
      axis(2, labels=F)
      if (i == 20) axis(1, las=1, tick=F, padj=-1)
      if (i %in% c(15,16,17,19)) abline(lm(all_fig1[,i] ~ all_fig1[,j]),
                                         lty=1, lwd=2)
    }
  }
}
cc <- .8; ll <- 4; lll <- 2.5
mtext(expression(paste("Chl (",mu,"g*cm"^{-2},")")),      2, las=0, line=ll, cex=cc, at=.92, outer=T)
mtext(expression(paste("Chl variance")),                   2, las=0, line=ll, cex=cc, at=.75, outer=T)
mtext(expression(paste("NCP (mgC*h"^{-1},")")),            2, las=0, line=ll, cex=cc, at=.58, outer=T)
mtext(expression(paste("NCP variance")),                   2, las=0, line=ll, cex=cc, at=.43, outer=T)
mtext(expression(paste("GPP (mgC*h"^{-1},")")),            2, las=0, line=ll, cex=cc, at=.25, outer=T)
mtext(expression(paste("GPP variance")),                   2, las=0, line=ll, cex=cc, at=.08, outer=T)
mtext("Crab",       1, las=0, line=1, cex=.8, at=.1,  outer=T)
mtext("Amphipods",  1, las=0, line=1, cex=.8, at=.3,  outer=T)
mtext("Ulva",       1, las=0, line=1, cex=.8, at=.5,  outer=T)
mtext("Kelp",       1, las=0, line=1, cex=.8, at=.7,  outer=T)
mtext(expression(paste("Water input (ml*min"^{-1},")")), 1, las=0, line=1.2, cex=.8, at=.89, outer=T)
dev.off()


# -------------------------------------------------------
# 8b. Figure 2: Organism biomass vs. treatments
# -------------------------------------------------------
# Rows: crab growth (%), amphipod biomass, Ulva biomass,
# kelp biomass, detritus biomass. Columns: treatments + flow.

ul     <- ulvadry[, c(3, 17)]
alg4   <- algae[algae$Sample == 1, ]
c2     <- c1[, c(1, 28)]; names(c2)[1] <- "Overall.bucket"
pod4   <- amphdry[, c(3, 17)];  names(pod4)[2] <- "podbio"
kelp4  <- kelpdry[, c(3, 17)];  names(kelp4)[2] <- "kelpbio"
det4   <- detdry[, c(3, 17)];   names(det4)[2]  <- "detbio"

all_fig2 <- cbind(dat1, nut1[1:144, c(8, 11, 10)])
all_fig2 <- join(all_fig2, c2,    type="left", by="Overall.bucket")
all_fig2 <- join(all_fig2, pod4,  type="left", by="Overall.bucket")
all_fig2 <- join(all_fig2, ul,    type="left", by="Overall.bucket")
all_fig2 <- join(all_fig2, kelp4, type="left", by="Overall.bucket")
all_fig2 <- join(all_fig2, det4,  type="left", by="Overall.bucket")
all_fig2 <- cbind(all_fig2, alg4$Ectocarpus.cover)
names(all_fig2)[20] <- "Ecto"
all_fig2$Ecto <- all_fig2$Ecto / 1000
# Set biomass to zero where organism was absent (was NA)
all_fig2[all_fig2$Crab == "No", 15] <- 0
all_fig2[all_fig2$Pods == "No", 16] <- 0
all_fig2[all_fig2$Ulva == "No", 17] <- 0

pdf(file.path(fig_dir, "Figure2_biomass_factors.pdf"), width=10, height=10)
par(mfrow=c(5, 5), mar=c(.3, .5, .2, .5), oma=c(3, 5.5, .5, 0))
for (i in 15:19) {
  for (j in c(6, 5, 4, 7, 8)) {
    y    <- c(0, max(all_fig2[, i], na.rm=T))
    xll  <- c(.5, 2.5); space <- c(1:2); cl <- c("white", "gray45")
    xl   <- c("No", "Yes")
    if (j == 7) { xl <- c("No","Yes","Grnd"); xll <- c(.5,3.5)
                  space <- c(1:3); cl <- c("white","gray80","gray45") }
    if (i == 15) { y <- c(0, 10);  ysig <- 9 }
    if (i == 16) { y <- c(0, .5);  ysig <- .45 }
    if (i == 17) { y <- c(0, 5);   ysig <- 4.5 }
    if (i == 19) { y <- c(0, 2);   ysig <- 1.7 }
    # Skip self-reference panels (organism plotted against its own treatment)
    self_ref <- (i==15 & j==6) | (i==16 & j==5) | (i==17 & j==4) | (i==18 & j==7)
    if (j != 8) {
      resp <- all_fig2[, i]; fac1 <- all_fig2[, j]
      if (self_ref) {
        plot(NA, xlim=xll, ylim=y, axes=F, xlab="", ylab="")
        if (j == 6) axis(2, las=2) else axis(2, labels=F)
        if (i == 19) axis(1, las=1, at=space, labels=xl, tck=0, padj=-1)
        else axis(1, labels=F, tck=0)
      } else {
        boxplot(resp ~ fac1, outline=F, col=cl, axes=F, at=space,
                xlim=xll, ylim=y, bty='l')
        if (j == 6) axis(2, las=2) else axis(2, labels=F)
        if (i == 19) { axis(1, las=1, at=space, labels=xl, tck=0, padj=-1)
                       abline(-1, 0) } else axis(1, labels=F, tck=0)
        q <- tapply(resp, list(fac1), mean, na.rm=T)
        points(space, as.vector(q), pch=1, col='black', cex=1.2)
        if (j==7 & i==15) text(c(1,2,3), ysig, labels=c("a","b","a"))
        if (j==6 & i==16) points(1.5, ysig, pch=8, col='black', cex=1.2)
        if (j==5 & i==19) points(1.5, ysig, pch=8, col='black', cex=1.2)
      }
    } else {
      plot(all_fig2[,j], all_fig2[,i], pch=1, cex=.6, ylab="", xlab="",
           xaxt="n", yaxt="n", bty='l', xlim=c(150,900), ylim=y)
      axis(2, labels=F)
      if (i == 20) axis(1, las=1, tick=F, padj=-1)
      if (j==8 & i==20) abline(lm(all_fig2[,i] ~ all_fig2[,j]), lty=1, lwd=2)
    }
  }
}
mtext(expression(paste("% crab growth (g)")), 2, las=0, line=2.5, cex=.8, at=.91, outer=T)
mtext(expression(paste("Amphipods (g)")),      2, las=0, line=2.5, cex=.8, at=.7,  outer=T)
mtext(expression(paste("Ulva (g)")),           2, las=0, line=2.5, cex=.8, at=.5,  outer=T)
mtext(expression(paste("Kelp (g)")),           2, las=0, line=2.5, cex=.8, at=.3,  outer=T)
mtext(expression(paste("Detritus (g)")),       2, las=0, line=2.5, cex=.8, at=.1,  outer=T)
mtext("Crab",       1, las=0, line=1, cex=.8, at=.1,  outer=T)
mtext("Amphipods",  1, las=0, line=1, cex=.8, at=.3,  outer=T)
mtext("Ulva",       1, las=0, line=1, cex=.8, at=.5,  outer=T)
mtext("Kelp",       1, las=0, line=1, cex=.8, at=.7,  outer=T)
mtext(expression(paste("Water input (ml*min"^{-1},")")), 1, las=0, line=1.2, cex=.8, at=.89, outer=T)
dev.off()


# -------------------------------------------------------
# 8c. Figure 3: Mesofauna community vs. treatments
# -------------------------------------------------------
# Rows: abundance, abundance variance, richness, richness
# variance, diversity, diversity variance.
# Columns: Crab, Amphipods, Ulva, Kelp, Water flow.

# Assemble plotting data frame (final time point + variances)
dat_ps4 <- ps_full[ps_full$Time == 4, ]
dat_ps4 <- dat_ps4[order(dat_ps4$bucket), ]
dat3    <- cbind(varg[, 1:11],
                 dat_ps4[, c(23, 24, 26, 27)],
                 varg[, 12:15])
dat3$Kelp <- factor(dat3$Kelp, levels=c("No","One piece","Ground up"))
# Reorder columns for plotting
dat3 <- dat3[, c(1:12, 15, 13, 16, 14, 17)]
all_fig3 <- dat3

pdf(file.path(fig_dir, "Figure3_mesofauna_community.pdf"), width=10, height=12)
par(mfrow=c(6, 5), mar=c(.3, .5, .2, .5), oma=c(3, 5.5, .5, 0))
for (i in 12:17) {
  for (j in c(6, 5, 4, 7, 8)) {
    y    <- c(0, max(all_fig3[, i], na.rm=T))
    xll  <- c(.5, 2.5); space <- c(1:2); cl <- c("white", "gray45")
    xl   <- c("No", "Yes")
    if (j == 7) { xl <- c("No","Yes","Grnd"); xll <- c(.5,3.5)
                  space <- c(1:3); cl <- c("white","gray80","gray45") }
    if (i == 12) { y <- c(0, 500);   ysig <- 480 }
    if (i == 13) { y <- c(0, .4);    ysig <- 0.39 }
    if (i == 14) { y <- c(0, 10);    ysig <- 9.5 }
    if (i == 15) { y <- c(0, 0.025); ysig <- 0.023 }
    if (i == 16) { y <- c(0, 1);     ysig <- 0.95 }
    if (i == 17) { y <- c(0, 0.012); ysig <- 0.011 }
    if (j != 8) {
      resp <- all_fig3[, i]; fac1 <- all_fig3[, j]
      keep <- !is.na(resp)
      n_lev     <- nlevels(droplevels(factor(fac1[keep])))
      space_use <- space[1:n_lev]; cl_use <- cl[1:n_lev]
      if (n_lev < length(space)) {
        plot.new()
      } else {
        boxplot(resp ~ fac1, outline=F, col=cl_use, axes=F, at=space_use,
                xlim=xll, ylim=y, bty='l')
        if (j == 6) axis(2, las=2) else axis(2, labels=F)
        axis(1, labels=F, tck=0)
        if (i == 17) axis(1, las=1, at=space, labels=xl, tck=0, padj=-1, lwd=0)
        q <- tapply(resp, list(fac1), mean, na.rm=T)
        points(space, as.vector(q), pch=1, col='black', cex=1.2)
        if (j==5 & i==12) points(1.5, ysig, pch=8, col='black', cex=1.2) # pod abund
        if (j==7 & i==13) text(c(1,2,3), ysig, labels=c("a","b","ab"))   # kelp abund-var
        if (j==5 & i==14) points(1.5, ysig, pch=8, col='black', cex=1.2) # pod richness
      }
    } else {
      plot(all_fig3[,j], all_fig3[,i], pch=1, cex=.6, ylab="", xlab="",
           xaxt="n", yaxt="n", bty='l', xlim=c(150,900), ylim=y)
      axis(2, labels=F)
      if (i == 17) axis(1, las=1, tick=F, padj=-1)
      if (j==8 & i==12) abline(lm(all_fig3[,i] ~ all_fig3[,j]), lty=1, lwd=2)
      if (j==8 & i==16) abline(lm(all_fig3[,i] ~ all_fig3[,j]), lty=1, lwd=2)
    }
  }
}
mtext("Abundance",           2, las=0, line=2.8, cex=.8, at=.92, outer=T)
mtext("Abundance variance",  2, las=0, line=2.8, cex=.8, at=.75, outer=T)
mtext("Richness",            2, las=0, line=2.8, cex=.8, at=.58, outer=T)
mtext("Richness variance",   2, las=0, line=2.8, cex=.8, at=.42, outer=T)
mtext("Diversity",           2, las=0, line=2.8, cex=.8, at=.25, outer=T)
mtext("Diversity variance",  2, las=0, line=2.8, cex=.8, at=.08, outer=T)
mtext("Crab",       1, las=0, line=1, cex=.8, at=.1,  outer=T)
mtext("Amphipods",  1, las=0, line=1, cex=.8, at=.3,  outer=T)
mtext("Ulva",       1, las=0, line=1, cex=.8, at=.5,  outer=T)
mtext("Kelp",       1, las=0, line=1, cex=.8, at=.7,  outer=T)
mtext(expression(paste("Water input (ml*min"^{-1},")")), 1, las=0, line=1.2, cex=.8, at=.89, outer=T)
dev.off()


# ============================================================
# END OF SCRIPT
# ============================================================
message("Analysis complete. Figures saved to: ", normalizePath(fig_dir))
