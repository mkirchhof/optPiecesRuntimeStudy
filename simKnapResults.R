# Results of the runtime simulations
# Reads in all result files and recombines them to the runtime tables


# ----- simKnap -----
cases <- expand.grid(nItems = c(100, 500, 1000, 5000, 10000), 
                     nCosts = c(3, 5, 10),
                     nChoose = c(10, 20, 50))
cases$Min <- NA
cases$Q1 <- NA
cases$Median <- NA
cases$Mean <- NA
cases$Q3 <- NA
cases$Max <- NA
for(i in seq(dim(cases)[1])){
  filename <- paste0("./results/", i, "_simKnapsack.RData", collapse = "")
  if(file.exists(filename)){
    load(filename)
    cases[i, 4:9] <- round(summary(res)) / 1000
  } else {
    cat(i, "not found\n")
  }
}
View(cases[,c(1,2,3,4,6,9)])



# ----- simKnapAutoCor -----
cases <- expand.grid(nItems = c(100, 500, 1000, 5000, 10000), 
                     nCosts = c(3, 5, 10))
cases$Min <- NA
cases$Q1 <- NA
cases$Median <- NA
cases$Mean <- NA
cases$Q3 <- NA
cases$Max <- NA
for(i in seq(dim(cases)[1])){
  filename <- paste0("./results/", i, "_simKnapsackAutoCor.RData", collapse = "")
  if(file.exists(filename)){
    load(filename)
    cases[i, 3:8] <- round(summary(res)) / 1000
  }
}
View(cases[,c(1,2,3,5,8)])



# ----- simKnapRealistic -----
cases <- expand.grid(nItems = c(100, 500, 1000, 5000, 10000), 
                     nCosts = c(3, 5, 10))
cases$Min <- NA
cases$Q1 <- NA
cases$Median <- NA
cases$Mean <- NA
cases$Q3 <- NA
cases$Max <- NA
for(i in seq(dim(cases)[1])){
  filename <- paste0("./results/", i, "_simKnapsackRealistic.RData", collapse = "")
  if(file.exists(filename)){
    load(filename)
    cases[i, 3:8] <- summary(res)
  }
}
View(cases[,c(1,2,3,5,8)])
