# Investigation_point_placement
Soriano-Elias, B.; Bach-Gómez, Anna & Molist Montañà, M. (In prep.) Investigating the influence of point placement strategies in the stability of modeled movement corridors: The case of the middle neolithic in the North-East of the Iberian Peninsula.

This repository contains all the data and scripts required to fully reproduce all analyses presented in the paper.

To get started:

Set the wd in a folder with all the data

&

Activate the following packages:

List of nescessary packages:

packages <- c("raster","leastcostpath","terra", "sf", "ggplot2","dplyr", "spatstat", "stars","corrplot","tidyr")

for (package in packages) { install.packages(package, character.only = TRUE) }

for (package in packages) { library(package, character.only = TRUE) }
