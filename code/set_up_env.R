library(renv)

renv::init()
renv::install("tidyverse")
renv::install("Matrix")
renv::install("data.tree")
renv::snapshot()
