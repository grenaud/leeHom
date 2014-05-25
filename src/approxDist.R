#!/usr/bin/env Rscript
library(fitdistrplus)
library(MASS)


args=(commandArgs(TRUE))

data <- read.table(args[1]);


df<-fitdistr(data$V1, "lognormal")

print(df);
