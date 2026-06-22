library(fossilx)
library(chronosphere)
library(divDyn)

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
       xlim = 84:93,
       ylim = c(0,2000))
lines(stages$mid, sampStg$occs)

