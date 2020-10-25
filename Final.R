dataDir  = "../data/"
fileName = "wvs_data.rds"

data = readRDS(paste0(dataDir, fileName))

#Loading packages.
library(mice) # For missing data descriptives
library(MASS) # For robust stats
library("dplyr")

dataV <-data %>% select(V7,V8,V10,V11,V23,V24,V45,V46,V48,V53,V55,V59,V60,V81,V96,V97,V121,V123,V139,V143,V181,V239,V248,V242,V240)

summary(dataV)
str(dataV)
names(dataV)

# Removing all the negative values.
dataV[dataV<0] <- NA

# Univariate outlier treatment.

bpOutliers <- function(x) {
  ## Compute inner and outer fences:
  iFen <- boxplot.stats(x, coef = 1.5)$stats[c(1, 5)]
  oFen <- boxplot.stats(x, coef = 3.0)$stats[c(1, 5)]
  
  ## Return the row indices of flagged 'possible' and 'probable' outliers:
  list(possible = which(x < iFen[1] | x > iFen[2]),
       probable = which(x < oFen[1] | x > oFen[2])
  )
}


# deleting all the probable univariate outliers for each variables.
uni_outliers <- lapply(dataV, FUN = bpOutliers)


out_clean <- function(y){
  indx=uni_outliers$y$probable
  dataV[indx, y] <- NA  
}

for (i in names(dataV)){
  out_clean(i)
}

cm <- colSums(is.na(dataV))
cm
pm <- colMeans(is.na(dataV))
pm

## Summarize PM:
range(pm)
mean(pm)
median(pm)

## Find variables with PM greater than 10%:
pm[pm > 0.1]

## Compute covariance coverage:
cc <- md.pairs(dataV)$rr / nrow(dataV)
cc

## Summarize coverages:
range(cc)

#multile imputation
## Define our own method vector:
meth        <- rep("norm", ncol(dataV))
names(meth) <- colnames(dataV)

dataV$V24 <- as.factor(dataV$V24)
dataV$V240 <- as.factor(dataV$V240)
dataV$V248 <- as.factor(dataV$V248)
dataV$V60 <- as.factor(dataV$V60)
dataV$V81 <- as.factor(dataV$V81)



meth["V24"]    <- "logreg"
meth["V240"]    <- "logreg"
meth["V248"] <- "polyreg"
meth["V60"] <- "polyreg"
meth["V81"] <- "polyreg"
meth["V242"]<- "norm"

## Use mice::quickpred to generate a predictor matrix:
predMat <- quickpred(dataV, mincor = 0.2, include = "V240")
#predMat <- quickpred(bfi, mincor = 0.2)

## Impute missing using the predictor matrix from above:
miceOut <- mice(data            = dataV,
                m               = 20,
                maxit           = 10,
                method          = meth,
                predictorMatrix = predMat,
                seed            = 235711)

summary(miceOut)


miceOut2 <-
  subset_datlist(datlist = miceOut,
                 select  = setdiff(colnames(dataV), c("V24", "V240","V248","V60","V81")),
                 toclass = "mids")

## Create list of multiply imputed datasets:
impList <- complete(miceOut2, "all")

mdOutliers <-
  function(data, critProb, statType = "mcd", ratio = 0.75, seed = NULL)
  {
    ## Set a seed, if one is provided:
    if(!is.null(seed)) set.seed(seed)
    
    ## Compute (robust) estimates of the mean and covariance matrix:
    stats <- cov.rob(x             = data,
                     quantile.used = floor(ratio * nrow(data)),
                     method        = statType)
    
    ## Compute robust squared Mahalanobis distances
    md <- mahalanobis(x = data, center = stats$center, cov = stats$cov)
    
    ## Find the cutoff value:
    crit <- qchisq(critProb, df = ncol(data))
    
    ## Return row indices of flagged observations:
    which(md > crit)
  }

olList <- lapply(impList, mdOutliers, critProb = 0.99)

## Count the number of times each observation is flagged as an outlier:
olCounts <- table(unlist(olList))

## Define the threshold for voting (will be 10 in this case):
thresh <- ceiling(miceOut$m / 2)

## Define a vector of row indices for outliers:
outs <- as.numeric(names(olCounts[olCounts >= thresh]))

## Exclude outlying observations from mids object:
miceOut3 <- subset_datlist(datlist = miceOut, # We're using the original imputations
                           subset  = setdiff(1 : nrow(dataV), outs),
                           toclass = "mids")


## Sanity check the imputations by plotting observed vs. imputed densities:
densityplot(miceOut3)

impList2 <- complete(miceOut3, "all")

save.image(file = "my_work_space.RData")
