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

# ----------

library(microbenchmark)


# sim - generates a simulation problem, solves it and tracks the runtime needed
#       for that
# Input:
#  nItems - integer, number of items / measurements to choose from 
#           (= J, in paper notation)
#  nCosts - integer, numer of costs to regard (= I, in paper notation)
#  reps - integer, number of times to run the simulation
# Output:
#  numeric, giving the runtime in milliseconds needed for solving the problems
sim <- function(nItems = 1000, nCosts = 5, reps = 100){
  moKnapsack <- function(cost, nItems = nItems, nCosts = nCosts, nChoose = nChoose){
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
      for(s in seq(nChoose - 1, 1, by = -1)){
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
    for(i in seq(nChoose, 2, by = -1)){
      res[i - 1] <- max(which(pred[, res[i], i]))
    }
    return(res)
  }
  
  
  nChoose <- 10
  cpuTime <- numeric(reps)
  for(i in 1:reps){
    # generate a problem:
    set.seed(i)
    cost <- matrix(0, ncol = nItems, nrow = nCosts)
    for(costInd in seq(nCosts)){
      cost[costInd, ] <- cumsum(rbinom(nItems, 20, 0.5) - 10)
    }
    
    # create predecessor structure (as array, for simplicity)
    # per column: which is a valid predecessor in an optimal solution
    pred <- array(FALSE, dim = c(nItems, nItems, nChoose))
    # third dimension: how many Items are before it
    # seq from 2, because if 1 is upper tri matrix, then the first element never
    # can be used
    pred[,,1] <- TRUE
    for(j in seq(2, nChoose)){
      # disallow each sheet to use the 5% of positions before it
      pred[1:(0.95*nItems),,j] <- upper.tri(matrix(TRUE, ncol = nItems, nrow = nItems))[-c(1:(0.05*nItems)),]
    }
    
    cpuTime[i] <- summary(microbenchmark(moKnapsack(cost, nChoose = 10, nCosts = nCosts, nItems = nItems), times = 1, unit = "ms"))$mean
    cat(i, "\n")
  }
  return(cpuTime)
}


# This is the loop that executes the simulation for all configurations given
# in the cases object
for(c in comp){
  res <- sim(nItems = cases[c, "nItems"], nCosts = cases[c, "nCosts"])
  save(res, file = paste0(c, "_simKnapsackAutoCor.RData"))
}