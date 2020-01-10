# This header makes parallel execution (on 5 cores) on a server cluster possible.
# If you don't need it, just ignore it and it will run single-threaded.
k = as.integer(Sys.getenv("PBS_ARRAYID"))
if(is.na(k)){
  k <- 1:5
}
cat("ID=", k, "\n")

# Choose cases for this script:
cases <- expand.grid(nItems = c(100, 500, 1000, 5000, 10000), 
                     nCosts = c(3, 5, 10))

comp <- matrix(c(1, 6, 15,
                 11, 2, 10,
                 7, 12, 5,
                 3, 8, 14,
                 13, 4, 9), byrow = TRUE, ncol = 3)
comp <- comp[k, ]

# -------

library(microbenchmark)


# simulateQuality - generates a quality time series by superimposing sinus waves
# Input:
#  foilLength - desired total length of the foil
#  nData - desired number of measurements along the foil
#  nWave - number of sinus waves to superimpose
# Output:
#  dataframe with two columns:
#   x - the positions of the measurements
#   q - the quality values at the positions
simulateQuality <- function(foilLength = 100, nData = 1000, nWave = 10){
  per <- runif(nWave, 5, 50)
  amp <- runif(nWave, 0, 0.3)
  xsh <- runif(nWave, -pi, pi)
  ysh <- runif(nWave, -0.25, 0.25)
  rnd <- rnorm(nData, 0, 0.1)
  x <- seq(from = 0, to = foilLength, length.out = nData + 1)[-1]
  res <- sapply(x, function(x) sum(ysh + amp * sin(x / per * 2 * pi + xsh))) + rnd
  return(data.frame(x = x, q = res))
}



# computeRelOut - takes a quality time series and transforms it into its 
#                 relative amount of outliers for a foil that ends at each of 
#                 the measurement positions
# Input:
#  measure - dataframe as created by simulateQuality
#  foilLength - total length of the foil
#  sheetLength - desired length of each sheet to be cut out
#  maxRelOutlier - maximum allowed amount of outliers per sheet before it is 
#                  regarded n.i.o.
# Output:
#  dataframe with two columns:
#   x - the positions of the measurements
#   r - the relative amount of outliers for a sheet that ends at position x
computeRelOut <- function(measure, foilLength = 100, sheetLength = 5, maxRelOutlier = 0.2){
  nData <- dim(measure)[1]
  nSheetPoints <- sheetLength / foilLength * nData
  outlier <- measure$q < -1 | measure$q > 1
  res <- data.frame(x = measure$x, r = NA)
  # In the beginning, no sheet can be selected
  res$r[res$x < sheetLength] <- Inf
  res$r[res$x >= sheetLength] <- filter(outlier, rep(1/nSheetPoints, nSheetPoints), sides = 1)[-(seq(nSheetPoints - 1))]
  res$r[res$r >= maxRelOutlier] <- Inf
  return(res)
}


# calcSheetCount - uses the naive algorithm to compute the maximum possible
#                  number of sheets in a relative outliers time series
# Input:
#  cost - matrix with one computeRelOutliers timeseries per row
#  foilLength - total length of the foil
#  sheetLength - length of each sheet
# Output:
#  integer, number of sheets that can be cut out of the foil
calcSheetCount <- function(cost, foilLength = 100, sheetLength = 5){
  nData <- dim(cost)[2]
  nSheetPoints <- sheetLength / foilLength * nData
  pos <- nSheetPoints
  sheetCount <- 0
  while(pos < nData){
    if(all(cost[, pos] < Inf)){
      sheetCount <- sheetCount + 1
      pos <- pos + nSheetPoints
    } else {
      pos <- pos + 1
    }
  }
  return(sheetCount)
}


# simulate - generates an optimization problem, solves it and returns the 
#            time needed to solve it
# Input:
#  nCosts - number of quality functions to regard simultaniously
#  foilLength - total length of the simulated foil
#  sheetLength - desired length of each sheet to be cut out
#  nData - number of simulated measurements along the foil
#  maxRelOutlier - maximum relative amount of outliers before a sheet is
#                  considered invalid
# Output:
#  numeric, the time needed to solve the optimization in milliseconds. 
#  Can return NA if no sheet could be found in the given simulated foil.
simulate <- function(nCosts = 10, foilLength = 100, sheetLength = 5, nData = 1000, maxRelOutlier = 0.2){
  moKnapsack <- function(cost, nItems = nData, nCosts = nCosts, nChoose = nChoose){
    for(c in seq(nCosts)){
      ownCost <- matrix(Inf, ncol = nItems, nrow = nChoose)
      for(row in seq(nChoose)){
        for(col in seq(nItems)){
          if(any(pred[, col, row])){
            if(row == 1){
              ownCost[row, col] <- cost[c, col]
            } else {
              predCost <- min(ownCost[row - 1, pred[, col, row]])
              pred[, col, row] <<- pred[, col, row] & 
                (ownCost[row - 1, ] == predCost)
              ownCost[row, col] <- cost[c, col] + predCost
            }
          }
        }
      }
      isLowest <- ownCost[nChoose, ] == min(ownCost[nChoose, ])
      pred[, !isLowest, nChoose] <<- FALSE
      for(s in rev(seq(nChoose)[-nChoose])){
        for(j in 1:nItems){
          # delete all its preds if it has no succeder
          if(!any(pred[j, ,s + 1])){
            pred[, j, s] <<- FALSE
          }
        }
      }
    }
    # Walk through pred to give out the optimal solution
    # start with the one that has the lowest owncosts
    res <- integer(nChoose)
    res[nChoose] <- which.min(ownCost[nChoose, ])
    for(i in rev(seq(nChoose)[-1])){
      res[i - 1] <- max(which(pred[, res[i], i]))
    }
    return(res)
  }
  
  # Generate cost matrix
  cost <- matrix(NA, nrow = nCosts, ncol = nData)
  for(i in seq(nCosts)){
    cost[i, ] <- computeRelOut(simulateQuality(foilLength = foilLength, nData = nData), 
                               foilLength = foilLength, sheetLength = sheetLength,
                               maxRelOutlier = maxRelOutlier)$r
  }
  
  # Find number of possible sheets:
  S <- calcSheetCount(cost, foilLength = foilLength, sheetLength = sheetLength)
  if(S > 0){
    # generate pred array
    nSheetPoints <- sheetLength / foilLength * nData
    pred <- array(FALSE, dim = c(nData, nData, S))
    pred[, nSheetPoints:nData, 1] <- TRUE
    for(i in seq(S)[-1]){
      pred[nSheetPoints:nData, nSheetPoints:nData, i] <- 
        upper.tri(matrix(TRUE, ncol = nData - nSheetPoints + 1, nrow = nData - nSheetPoints + 1), 
                  diag = FALSE)
    }
    
    # find optimal sheets:
    return(summary(microbenchmark(moKnapsack(cost, nChoose = S, nCosts = nCosts, nItems = nData), times = 1), unit = "ms")$mean)
  } else {
    return(NA)
  }
}

# main - a wrapper that executes simulate() several times
# Input:
#  nCosts - number of quality functions to regard simultaniously
#  foilLength - total length of the simulated foil
#  sheetLength - desired length of each sheet to be cut out
#  nData - number of simulated measurements along the foil
#  reps - number of repititions of the simulation
# Output:
#  numeric of length reps, includes time in milliseconds needed to solve each 
#  problem
main <- function(nCosts = 10, foilLength = 100, sheetLength = 5, nData = 1000, reps = 100){
  set.seed(123)
  res <- numeric(reps)
  i <- 1
  while(i <= reps){
    cat(i, "\n")
    res[i] <- simulate(nCosts = nCosts, foilLength = foilLength,
                       sheetLength = sheetLength, nData = nData, maxRelOutlier = 0.2)
    if(!is.na(res[i]))
      i <- i + 1
  }
  return(res)
}



# This is the loop that executes the simulation for all configurations given
# in the cases object
for(c in comp){
  res <- main(nData = cases[c, "nItems"], nCosts = cases[c, "nCosts"])
  save(res, file = paste0(c, "_simKnapsackRealistic.RData"))
}


# In the following, there is some testing code to see how the simulated quality
# looks and how many sheets can be found in it.
#
#
# test <- main()
# 
# for(i in 1:1000){
#   Sys.sleep(1)
#   par(mfrow = c(2, 1), mar = c(1, 3, 1, 1))
#   sim <- simulateQuality()
#   rel <- computeRelOut(sim)
#   rel$r[rel$r == Inf] <- NA
#   plot(sim, ylim = c(-2, 2))
#   abline(h = c(-1, 1))
#   plot(rel)
#   abline(h = 0.2)
# }
# 
# simulateOnlyS <- function(nCosts = 10, foilLength = 100, sheetLength = 5, nData = 1000){
#   # Generate cost matrix
#   cost <- matrix(NA, nrow = nCosts, ncol = nData)
#   for(i in seq(nCosts)){
#     cost[i, ] <- computeRelOut(simulateQuality(foilLength = foilLength, nData = nData),
#                                foilLength = foilLength, sheetLength = sheetLength)$r
#   }
#   S <- calcSheetCount(cost, foilLength = foilLength, sheetLength = sheetLength)
#   return(S)
# }
# mainOnlyS <- function(nCosts = 10, foilLength = 100, sheetLength = 5, nData = 1000, reps = 100){
#   cat(nCosts, nData, "\n")
#   set.seed(123)
#   res <- numeric(reps)
#   i <- 1
#   while(i <= reps){
#     res[i] <- simulateOnlyS(nCosts = nCosts, foilLength = foilLength,
#                             sheetLength = sheetLength, nData = nData)
#     if(res[i] > 0)
#       i <- i + 1
#   }
#   return(res)
# }
# Sdata <- apply(cases, 1, function(x) mainOnlyS(nCosts = x[2], nData = x[1]))
# apply(Sdata, 2, summary)
# # Combine over all nDatas, only differ between nCosts
# apply(matrix(Sdata, ncol = 3), 2, summary)
