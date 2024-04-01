## FullStats.py

Used to recreate statistics at CBSA level and compile statistics at GBUA level

## ScalingStats.py
Used to calculate scaling laws for various population definitions (dasymetric refinement based on number of buildings, building footprint, or indoor area) and different geo and temporal completeness.

## ScalingStats_robust.py

We calculated the scaling laws for variations in how patches are defined (see paper for details).

## find_bufa_bui.py

Finds the total BUFA (footprint area) or BUI (indoor area) to calculate dasymetricly refined populations at GBUA or CBSA level.


## nl_fit.r

R code to calculate the non-linear fit of our data (y = A*P^B, where y is a statistic, such as indoor area, and P is population. A, and B are fitting parameters)
