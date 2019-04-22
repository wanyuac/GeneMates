# Processing the output of function findPhysLink.
#
# Copyright 2019 Yu Wan (wanyuac@gmail.com)
# Licensed under the Apache License, Version 2.0
# First and the latest edition: 22 April 2019

out_pl <- readRDS("Example/physlink.rds")
pl <- out_pl$assoc
pl <- subset(pl, s_a == 1 & s_d > 0)  # strong evidence of physical linkage

out_lmm <- readRDS("Example/lmm.rds")
lmms <- out_lmm$assoc
lmms <- subset(lmms, s_a == 1)  # significant positive associations
