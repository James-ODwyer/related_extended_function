This R script contains extentions to two functions within the related R package. 
`familysim` and `compareestimators`

The R script here extends out these two functions to allow for more than 100 loci to be analysed at once when simulating data using related. 
The script additions are very simple loops and have not been extensively tested. The extended functions work well on binary SNP data but have not been tested on other forms of
genetic marker e.g. microsatellites. 

This script also expands on the `compareestimators` function to allow for all relatedness estimators to be compared to each other and not just the four moment estimators the original
function allowed. 
