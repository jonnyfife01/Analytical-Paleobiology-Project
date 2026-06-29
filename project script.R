library(fossilx)
library(chronosphere)
library(divDyn)

# Set working directory
setwd("C:/Users/jonny/Documents/R projects/Analytical Paleobiology/Project")

# Pulling the data from pbdb
refined <- fossilx::pbdb_extend(x = "latest",
                             tax.level = "genus",
                             tax.combine = "clgen",
                             include = list("class:Mammalia"),
                             env.categories = "divDyn",
                             strat=c("stg_1.0", "ten_1.0", "stb")
                             )

# Rows containing marine fauna
index <- c(
  
  grep("marine", refined$taxon_environment),
  grep("oceanic", refined$taxon_environment),
  grep("coastal", refined$taxon_environment)
  
)

marmam <- refined[index, ]

# Removing nonmarine taxa from the data, update as needed
# Currently Mustelids are the only definitely nonmarine group
nonmarine <- c("Mustelidae")

marmam <- marmam[marmam$family != nonmarine, ]

# Get stages abd keys from divdyn
data("stages")
data("keys")

# Stage entries
stgMin <- categorize(marmam[, "early_interval"], keys$stgInt)
stgMax <- categorize(marmam[, "late_interval"], keys$stgInt)

# convert to numeric
stgMin <- as.numeric(stgMin)
stgMax <- as.numeric(stgMax)

# empty container
marmam$stg <- rep(NA, nrow(marmam))

# Select entries, where
stgCondition <- c(
  # The early and late interval fields indicate the same stg
  which(stgMax == stgMin),
  # or the late_interval field is empty
  which(stgMax == -1))

# in these entires, use the stg indicated by the early_interval
marmam$stg[stgCondition] <- stgMin[stgCondition]


# find number of occurrences and collections
sampStg <- binstat(marmam, tax = "genus", bin = "stg", col = "collection_no", duplicates = FALSE)

# plot diversity
tsplot(stages, boxes = c("series", "sys"),
       shading = "stg",
       xlim = 84:95,
       ylim = c(0,2000))
lines(stages$mid, sampStg$occs,
      lwd = 2)

################
# Strontium Data

strontium <- read.csv("data/Strontium.csv", skip = 6)

#Remove empty columns
strontium <- Filter(function(x)!all(is.na(x)), strontium)

# Get stage for each row
stgMin <- categorize(strontium$AGE...STAGE, keys$stgInt)
stgMax <- categorize(strontium$AGE...STAGE, keys$stgInt)


# convert to numeric
stgMin <- as.numeric(stgMin)
stgMax <- as.numeric(stgMax)

# empty container
strontium$stg <- rep(NA, nrow(strontium))

# Select entries where
stgCondition <- c(
  # The early and late interval fields indicate the same stg
  which(stgMax == stgMin),
  # or the late_interval field is empty
  which(stgMax == -1))

# use the stg indicated
strontium$stg[stgCondition] <- stgMin[stgCondition]

# find mean 87/86 Sr for each stage
SrRatio <- c()

for(i in 1:95){
  SrRatio[i] <- mean(strontium$X87Sr.86Sr.mean[strontium$stg == i], na.rm = TRUE)
}

# plot strontium
tsplot(stages, boxes = c("series", "sys"),
       shading = "stg",
       xlim = 84:95,
       ylim = c(.705,.71))
lines(stages$mid, SrRatio,
      lwd = 2, col = "brown")

# Correlation between the two
corr <- cor.test(SrRatio, sampStg$occs)

corr

stats <- lm(sampStg$occs ~ SrRatio)

summary(stats)

###########################################
# Check for autocorrelations
# in marine mammals
mam.div <- diff(sampStg$occs,)
acf(mam.div[84:94])

# autocorrelation in strontium
sr.div <- diff(SrRatio)
acf(sr.div[69:94])

# Plot them together
plot(SrRatio, sampStg$occs,
     xlim = c(0.7075, 0.7095),
     pch = 16,
     col = "darkgreen",
     xlab = expression({}^{87}*Sr/{}^{86}*Sr),
     ylab = "Marine Mammal Fossil Occurrences")
abline(stats)
