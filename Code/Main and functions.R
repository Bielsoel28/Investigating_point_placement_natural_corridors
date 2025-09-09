#################################################################################### 
###################### Main and functions for: ##################################### 
### Investigating the influence of point placement strategies in the stability of ## 
### modeled movement corridors: The case of the middle neolithic in the North-East #
### of the Iberian Peninsula. ######################################################
#### Biel Soriano-Elias, Anna Bach G?mez, Miquel Molist Montany? ###################
####################################################################################


# set_wd() in the desired folder
setwd("XXXXX/XXXXX/XXXXX")

# Setting up necessary packages
packages <- c("raster","leastcostpath","terra", "sf", "ggplot2","dplyr", "spatstat", "stars","corrplot","tidyr")

for (package in packages) {
  library(package, character.only = TRUE)
}

### Functions needed ###

Edge_model <- function(y,q,w,dir) { 
  
  #Creating the folder to store the shapes of the models
  dir.create(dir)
  
  #Loading of the input DEM or "y"
  elevacio <- y #raster for points
  res_rast <- ext(elevacio)
  
  #Creation of the limit edge points
  
  #Creation of the raster's polygon
  ext_red <- c(res_rast[1]+1,res_rast[2]-1,res_rast[3]+1,res_rast[4]-1)
  elevacio_pol <- crop(elevacio, ext_red)
  pol_ras <- st_union(st_as_sf(as.polygons(elevacio_pol, dissolve = TRUE)))
  
  #Converting raster's polygon to edge points
  mpol_ras <- st_cast(pol_ras, to = "POLYGON") #multipolygon to polygon 
  pol_line <- st_cast(mpol_ras, to = "LINESTRING") #polygon to line 
  densitat <- 1 / q #calculation of point density depending on q
  p_ras <- st_line_sample(pol_line, density = densitat) #line to points depending of density of q
  punts_loop <- p_ras[[1]] #loading and division of edge points for loop functionality
  punts_loop <- st_sf(geometry = st_sfc(punts_loop)) #conversion of array to sf
  punts_loop <- st_cast(punts_loop, "POINT")
  
  #Application of the cost function (creation of the cost surface layer from package "leastcostpath")
  mde_base <- create_slope_cs(elevacio, cost_function = w, neighbours = 16, exaggeration = FALSE) #cost calculation
  
  print("Step 1 completed") 
  
  
  #Computation of LCPs for natural corridors
  
  #Creation of empty object to store LCPs
  xarxa_camins_total <- st_sfc(st_linestring(), crs = crs(elevacio)) #creation of line string sfc object
  xarxa_camins_total <- st_sf(geometry = xarxa_camins_total) #creation of line string sf object
  xarxa_camins_total <- xarxa_camins_total[-1, ] #emptying of the line string sf object
  
  #Creation of LCPs 
  
  #Loop for LCP calculation
  for (i in 1:nrow(punts_loop)) {
    
    #Recovery of coordinates of origin point
    punt_objectiu <- punts_loop[i, ]
    
    #Recovery of the destination point
    punts_desti <- punts_loop[-i, ]
    
    #LCPs calculation 
    xarxa_camins <- create_lcp(mde_base, punt_objectiu, punts_desti, cost_distance = FALSE) #lcps calcultaion
    xarxa_camins <- st_set_crs(xarxa_camins,crs(elevacio)) #seting crs
    
    #Sumation of all LCPs without attributes
    xarxa_camins_total <- rbind(xarxa_camins_total, xarxa_camins[ , -c(1,2,3,4)])
    
    rm(xarxa_camins) #saving memory
    
  }
  
  print("Step 2 completed") 
  
  #Saving of the LCPs to Shapefile
  nom_arxiu <- paste(substitute(y),as.character(q),substitute(w),"xarxa_camins_total.gpkg") #File name creation
  xarxa_camins_total <- xarxa_camins_total[!st_is_empty(xarxa_camins_total), ] #filter of empty geometries
  st_write(xarxa_camins_total, dsn = file.path(dir,nom_arxiu)) #Saving of LCPs object from sf object to GPKG
  
  print("Process ended succesfully")
}
Edge_rotate_model <- function(y,q,q2,w,dir) { 
  
  #Creating the folder to store the shapes of the models
  dir.create(dir)
  
  #Loading of the input DEM or "y" 
  elevacio <- y
  res_rast <- ext(elevacio)
  
  #Creation of the limit edge points
  
  #Creation of the raster's polygon
  ext_red <- c(res_rast[1]+1,res_rast[2]-1,res_rast[3]+1,res_rast[4]-1)
  elevacio_pol <- crop(elevacio, ext_red)
  pol_ras <- st_union(st_as_sf(as.polygons(elevacio_pol, dissolve = TRUE)))
  
  #Converting raster's polygon to edge points
  pol_line <- st_cast(pol_ras, to = "LINESTRING") #polygon to line 
  densitat <- 1 / q2 #calculation of point density depending on q2
  p_ras <- st_line_sample(pol_line, density = densitat) #line to points depending of density of q
  punts_loop <- p_ras[[1]] #loading and division of edge points for loop functionality
  punts_loop <- st_sf(geometry = st_sfc(punts_loop)) #conversion of array to sf
  punts_loop <- st_cast(punts_loop, "POINT")
  
  #Application of the cost function (creation of the cost surface layer from package "leastcostpath")
  mde_base <- create_slope_cs(elevacio, cost_function = w, neighbours = 16, exaggeration = FALSE) #cost calculation
  
  print("Step 1 completed") 
  
  #Computation of LCPs for natural corridors
  
  #Setting up number of movements
  m <- q/q2
  
  #Selecting point for loop of movement
  
  for (p in 1:m) { 
    
    time_taken <- system.time({ 
      #Identifying rows of the points, subsisting and sf conversion
      maxm <- seq(from=1*p,to=nrow(punts_loop),by=m) #identifying rows
      punts_loop_mov <- punts_loop[maxm, ] #subletting points
      
      #Creation of empty object to store LCPs
      xarxa_camins_total <- st_sfc(st_linestring(), crs = crs(elevacio)) #creation of line string sfc object
      xarxa_camins_total <- st_sf(geometry = xarxa_camins_total) #creation of line string sf object
      xarxa_camins_total <- xarxa_camins_total[-1, ] #emptying of the line string sf object
      
      #Loop for LCP calculation
      for (i in 1:nrow(punts_loop_mov)) {
        
        #Recovery of coordinates of origin point
        punt_objectiu <- punts_loop_mov[i, ]
        
        #Recovery of the destination points
        punts_desti <- punts_loop_mov[-i, ]
        
        #LCPs calculation 
        xarxa_camins <- create_lcp(mde_base, punt_objectiu, punts_desti, cost_distance = FALSE) #lcps calcultaion
        xarxa_camins <- st_set_crs(xarxa_camins,crs(elevacio)) #seting crs
        
        #Sumation of all LCPs without attributes
        xarxa_camins_total <- rbind(xarxa_camins_total, xarxa_camins[ , -c(1,2,3,4)])
        
        rm(xarxa_camins) #saving memory
        gc()
      }
      
      #Saving of the LCPs to Shapefile
      nom_arxiu <- paste(substitute(y),as.character(q),substitute(w),"posicio",as.character(p),"xarxa_camins_total.gpkg", sep = "_") #File name creation
      xarxa_camins_total <- xarxa_camins_total[!st_is_empty(xarxa_camins_total), ] #dilter of empty geometries
      st_write(xarxa_camins_total, dsn = file.path(dir,nom_arxiu)) #Saving of LCPs object from sf object to Shapefile
      rm(xarxa_camins_total) #erasing the data frame where the paths are stored
      
    })
    
    print(paste("final model",as.character(p)))
    
    # Store the user time
    time_results <<- rbind(time_results, data.frame(time_in_seconds = time_taken["elapsed"]))
    
    # Set the row name as the current value of p
    rownames(time_results)[nrow(time_results)] <<- as.character(p)
    
  }
  
  print("Process ended succesfully")
  
}

Kernel_function <- function(x,y,dir) { 
  
  #Creating the folder to store the rasters
  dir.create(dir)
  
  # Obtaining raster extent from "y"
  ras_ex <- y |>
    sf::st_bbox()
  
  # Loading LCPs network
  lcps_reduced <- st_read(x)
  
  ## Compute density in spatstat (KDE) for each iteration 
  #Setting up the fuction
  
  lcp_kernel <- function(lcp) {
    
    # Define the window (class owin) for the density analysis    
    win <- spatstat.geom::owin(c(ras_ex[1],ras_ex[3]),
                               c(ras_ex[2],ras_ex[4]))
    
    # Define ppp object based on point locations
    lcps_pts_pp <- spatstat.geom::ppp(x = sf::st_coordinates(lcp)[,1], 
                                      y = sf::st_coordinates(lcp)[,2], 
                                      window = win)
    
    # Calculate density 
    lcps_density <- spatstat.explore::density.ppp(
      x = lcps_pts_pp,
      sigma = 20, 
      kernel = "gaussian",
      eps = 1 # cell size
    )
  }

  # Application of the function
  lcp_kernel_data <- lcp_kernel(lcps_reduced)
  
  ## Normalize values
  min_val <- min(lcp_kernel_data$v, na.rm = TRUE)
  max_val <- max(lcp_kernel_data$v, na.rm = TRUE)
  
  lcp_kernel_data$v <- (lcp_kernel_data$v - min_val) / (max_val - min_val)
  
  ## Save the density as GEOtiff for further analysis
  # convert it to a terra raster
  lcps_density_rast <- stars::st_as_stars(
    lcp_kernel_data) |> 
    sf::st_set_crs(25831) |>
    terra::rast()
  
  # Save
  terra::writeRaster(lcps_density_rast,filename = file.path(dir,paste(sub("\\.[^.]*$", "", basename(x)),"KDE_raster.tiff")))
  
}
Evaluation_function <- function(x,n,y,dir) {
  
  #Creating the folder to store the shapes of the models
  dir.create(dir)
  
  #Loading of the rasters
  raster_obj <- y
  
  rasters <- list.files(x, full.names = TRUE)
  
  #Creation of the raster's polygon
  raster_camins <- rast(rasters[1])
  pol_ras <- st_union(st_as_sf(as.polygons(raster_obj, dissolve = TRUE)))
  
  #Creation of n equidistant control points
  
  #Obtaining raster extension
  extent_raster <- ext(raster_camins) 
  
  #Identify and get Cells ID of n control points
  pc_y <- seq(from = extent_raster[3], to = extent_raster[4], length.out = n) #limit Y calculation
  pc_x <- seq(from = extent_raster[1], to = extent_raster[2], length.out = n) #limit X calculation
  xy_pc <- expand.grid(x = pc_x[2:(n-1)], y = pc_y[2:(n-1)]) #create a grid of control points
  pc_sf <- st_as_sf(xy_pc, coords = c("x","y"), crs = crs(raster_camins)) #control points to sf
  intersects_list <- st_intersects(pc_sf, pol_ras[[1]], sparse = TRUE) # logical for point within the raster
  pc_sf2 <- pc_sf[lengths(intersects_list) > 0, ] #extracting points
  
  #creation of an empty data frame to store data
  df_combinat <- data.frame(ID_CP = 1:nrow(pc_sf2)) #n rows
  
  #Loop to store all the raster values for the same control points
  for(i in rasters) { 
    
    #Loading of the input ZPM map or "x"
    raster_camins <- rast(i)
    
    #Extracting the values
    df <- terra::extract(raster_camins, pc_sf2) #Extracting the values
    df[is.na(df)] <- 0 #Converting NA to 0
    df_combinat <- cbind(df_combinat, df[ ,2]) #storing the values
    nom_col <- gsub("[a-zA-Z._/ ]|200_tobler|500_tobler|50_tobler|60_tobler||90_tobler||75_tobler|", "", i) #setting up column name
    colnames(df_combinat)[ncol(df_combinat)] <- nom_col #applying new name to column
    
    #Freeing memory
    rm(raster_camins)
    
    print(paste("final",substitute(i)))
    
  }
  
  
  #ordering data_frame to calculate difference of variation for point in pairs of models
  t_df <- t(df_combinat)
  t_df <- t_df[-1, ]
  t_df <- t_df[order(as.numeric(rownames(t_df))), ]
  df_2 <- data.frame(matrix(ncol = 0, nrow = nrow(t_df)-1)) 
  df_3 <- data.frame(matrix(ncol = 0, nrow = nrow(t_df)-1))   
  
  
  # Creation of boxplot of values #
  
  # Data frame to long format
  t_df_long <- as.data.frame(t_df) %>%
    pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value")
  
  # Setting up the name of the variable
  t_df_long$Variable <- gsub("V", "",  t_df_long$Variable, ignore.case = TRUE)
  
  # Convert Variable to numeric for proper ordering
  t_df_long$Variable <- as.numeric(t_df_long$Variable)
  
  # Plotting of the results
  boxplot <- ggplot(t_df_long, aes(x = factor(Variable, levels = sort(unique(Variable))), y = Value)) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = "Natural corridors values per control point", x = "Control Points", y = "Natural corridors values")
  
  ggsave(file.path(dir,"boxplot.tiff"),boxplot, width = 12, height = 8) # Saving the results
  
  #calculating the difference between values of points for each pair of models
  for (c in 1:ncol(t_df)) {
    
    t_df_in <- as.data.frame(t_df[ , c], colnames = colnames(t_df[ , 3]))
    df_2_in <- abs(diff(t_df_in[ ,1]))
    nom_col <- paste0("valors_diff_punt_", as.character(c))
    df_2[[nom_col]] <- df_2_in
    row.names(df_2) <- rownames(t_df)[-1]
    
  }
  
  #calculating the difference between values of points for every model and the first
  for (c3 in 1:ncol(t_df)) {
    
    t_df_in2 <- as.data.frame(t_df[ , c3], colnames = colnames(t_df[ , 3]))
    df_3_in <- abs(t_df_in2[-1,1] - t_df_in2[1,1])
    nom_col2 <- paste0("valors_dif_abs_punt_", as.character(c3))
    df_3[[nom_col2]] <- df_3_in
    row.names(df_3) <- rownames(t_df)[-1]
    
  }
  
  #creation of means of the substracted values
  df_mean2 <- colMeans(t(df_2))
  df_mean3 <- colMeans(t(df_3))
  df_sd2 <- apply(t(df_2), 2, sd)
  df_sd3 <- apply(t(df_3), 2, sd)
  df_sd2_sum <- df_mean2 + df_sd2
  df_sd2_res <- df_mean2 - df_sd2
  df_sd3_sum <- df_mean3 + df_sd3
  df_sd3_res <- df_mean3 - df_sd3
  
  # Storing the results in a data frame
  df_mean_final <- data.frame(M_diff = df_mean2, 
                              M_diff_abs = df_mean3,
                              SD_diff = df_sd2,
                              Sum_M_SD_diff = df_sd2_sum,
                              Res_M_SD_diff = df_sd2_res,
                              SD_diff_abs = df_sd3,
                              Sum_M_SD_diff_abs = df_sd3_sum,
                              Res_M_SD_diff_abs = df_sd3_res
  )
  
  write.csv(df_mean_final, file = file.path(dir,"means_pc.csv")) # Saving the results
  
  # Function to select unique values
  select_unique_rows <- function(data, column) {
    # Keep rows where the column value is not duplicated
    data[!duplicated(data[[column]]) & !duplicated(data[[column]], fromLast = TRUE), ]
  }
  
  # Sequential evaluation graphic creation #
  # Data frame to store results
  df_line_1_pre <- data.frame(Mean = df_mean2, 
                              M_SD_sum = df_sd2_sum,
                              M_SD_minus = df_sd2_res
  )
  
  df_line_1_pre$Time <- time_results$time_in_seconds[-1]  # Adding the time column to the data frame
  
  # Storing only unique values
  df_line_1_pre_2 <- select_unique_rows(df_line_1_pre, "Mean")
  df_time_pre_1 <- data.frame(time = df_line_1_pre_2[, 4], row.names = row.names(df_line_1_pre_2)) #picking up time
  df_line_1_pre_2 <- df_line_1_pre_2[ ,-4] #erasing time from df
  
  # Reordering values to plot
  df_line_1 <- data.frame(
    x= rep(row.names(df_line_1_pre_2)),
    y = c(df_line_1_pre_2[ ,1],df_line_1_pre_2[ ,2], df_line_1_pre_2[ ,3]),
    variable = rep(colnames(df_line_1_pre_2), each = nrow(df_line_1_pre_2))
  )
  
  # Creating df for time in ggplot  
  df_time_1 <-data.frame(
    x= rep(row.names(df_time_pre_1)),
    y = c(df_time_pre_1[ ,1]),
    variable = rep(colnames(df_time_pre_1), each = nrow(df_time_pre_1))
  )
  max_time_1 <- max(df_time_1$y, na.rm = TRUE) #establishing limits for plot
  df_time_1$x <- as.numeric(as.character(df_time_1$x))
  df_time_1$x <- -df_time_1$x 
  
  df_line_1$x <- as.numeric(as.character(df_line_1$x))
  df_line_1 <- df_line_1[order(df_line_1$x), ]
  df_line_1$x <- -df_line_1$x  # Multiply by -1 to invert values
  
  # Setting up X axis labels
  by <- as.numeric(row.names(df_line_1_pre))
  by_bo <- by[2] - by[1]
  
  # Ploting results
  df_line_plot_1 <- ggplot() +
    # Add secondary y-axis line first (so it stays in the background)
    geom_line(data = df_time_1, aes(x = x, y = y / max_time_1 * 0.25), 
              color = "red", linetype = "dashed", linewidth = 1) +
    
    # Add primary y-axis lines
    geom_line(data = df_line_1, aes(x = x, y = y, color = variable, group = variable), linewidth = 1) +
    
    # X-axis
    scale_x_continuous(breaks = seq(min(df_line_1$x), max(df_line_1$x), by = by_bo), labels = abs) + 
    
    # Primary Y-axis
    scale_y_continuous(
      limits = c(-0.01, 0.25),
      name = "Variation on Natural corridors values",
      sec.axis = sec_axis(~ . * (max_time_1 / 0.25), 
                          name = "Time (seconds)", 
                          breaks = scales::pretty_breaks(n = 5), 
                          labels = scales::comma)
    ) +
    
    # Colors
    scale_color_manual(values = c("grey", "grey", "black")) +
    
    # Labels
    labs(
      x = "Models (Numbers referring to in-between distance of points)", 
      y = "Variation on Natural corridors values", 
      color = "Legend"
    ) +
    
    # Theme
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # View the plot
  ggsave(file.path(dir,"Graphic_mean_diff.tiff"),df_line_plot_1, width = 12, height = 8) # Saving the files
  
  # Evaluation to fist point graphic creation #
  keep_first_unique <- function(df, column_name) {
    df[!duplicated(df[[column_name]]), ]
  }
  
  df_line_2_pre <- data.frame(Mean = df_mean3,
                              M_SD_sum = df_sd3_sum,
                              M_SD_minus = df_sd3_res
  )
  
  df_line_2_pre$Time <- time_results$time_in_seconds[-1]  # Adding the time column to the data frame
  
  # Storing only unique values
  df_line_2_pre_2 <- keep_first_unique(df_line_2_pre, "Mean")
  df_time_pre_2 <- data.frame(time = df_line_2_pre_2[, 4], row.names = row.names(df_line_2_pre_2)) #picking up time
  df_line_2_pre_2 <- df_line_2_pre_2[ ,-4] #erasing time from df
  
  # Reordering values to plot
  df_line_2 <- data.frame(
    x= rep(row.names(df_line_2_pre_2)),
    y = c(df_line_2_pre_2[ ,1],df_line_2_pre_2[ ,2], df_line_2_pre_2[ ,3]),
    variable = rep(colnames(df_line_2_pre_2), each = nrow(df_line_2_pre_2))
  )
  
  # Creating df for time in ggplot  
  df_time_2 <-data.frame(
    x= rep(row.names(df_time_pre_2)),
    y = c(df_time_pre_2[ ,1]),
    variable = rep(colnames(df_time_pre_2), each = nrow(df_time_pre_2))
  )
  max_time_2 <- max(df_time_2$y, na.rm = TRUE) #establishing limits for plot
  df_time_2$x <- as.numeric(as.character(df_time_2$x))
  df_time_2$x <- -df_time_2$x 
  
  df_line_2$x <- as.numeric(as.character(df_line_2$x))
  df_line_2 <- df_line_2[order(df_line_2$x), ]
  df_line_2$x <- -df_line_2$x  # Multiply by -1 to invert values
  
  # Ploting results
  df_line_plot_2 <- ggplot() +
    geom_line(data = df_time_2, aes(x = x, y = y / max_time_2 * 0.25), 
              color = "red", linetype = "dashed", linewidth = 1) +
    geom_line(data = df_line_2, aes(x = x, y = y, color = variable, group = variable), linewidth = 1) +
    scale_x_continuous(breaks = seq(min(df_line_2$x), max(df_line_2$x), by = by_bo), labels = abs) + 
    scale_y_continuous(
      limits = c(-0.01, 0.25),
      name = "Variation on Natural corridors values",
      sec.axis = sec_axis(~ . * (max_time_2 / 0.25), 
                          name = "Time (seconds)", 
                          breaks = scales::pretty_breaks(n = 5), 
                          labels = scales::comma)
    ) +
    scale_color_manual(values = c("grey", "grey", "black")) +
    labs(
      x = "Models (Numbers referring to in-between distance of points)", 
      y = "Variation on Natural corridors values", 
      color = "Legend"
    ) +
    
    # Theme
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # View the plot
  ggsave(file.path(dir,"Graphic_mean_diff_abs.tiff"),df_line_plot_2, width = 12, height = 8) 
  
  print("Proc?s finalitzat amb ?xit")
  
}

FETE_model <- function(y,q,w,dir) { 
  
  #Creating the folder to store the shapes of the models
  dir.create(dir)
  
  #Loading of the input DEM or "y" and rivers layers or "z" and convertion to terra
  elevacio <- y
  
  #Creation of the raster's polygon
  pol_ras <- st_union(st_as_sf(as.polygons(elevacio, dissolve = TRUE)))
  
    #Converting raster's polygon to middle points
    distancia <- q #loading distance between points
    punts_dins <- st_make_grid(pol_ras, cellsize = distancia, what = "centers") #creating points
    punts_dins2 <- punts_dins[st_intersects(punts_dins, pol_ras, sparse = FALSE)] #selecting point thata coincide with the raster
    punts_loop2 <- do.call(rbind, lapply(punts_dins2, st_coordinates)) #loading and division of edge points for loop functionality
  
  print("Step 1 completed")  
  
  
  #Application of the cost function (creation of the cost surface layer from package "leastcostpath")
  
  #Creation of base cs map 
  mde_base <- create_slope_cs(elevacio, cost_function = w, neighbours = 16, exaggeration = FALSE)
  
  #if condition for covering NA values for lcp calculation
  if (dis == "EXTERIOR"){  
    mde_base <- create_slope_cs(elevacio, cost_function = w, neighbours = 16, exaggeration = FALSE) #cost calculation
    buff_punts <- st_buffer(p_ras, dist = 25) #buffer for ensuring points values
    elevacio2 <- subst(elevacio, NA, 1) #new mde without NA for cost of points
    mde_base2 <- create_slope_cs(elevacio2, cost_function = w, neighbours = 16, exaggeration = FALSE) #creation of cost raster for points
    mde_base <- replace_values(mde_base, mde_base2, buff_punts) #replacing cost rasters values with point modified values
  }
  
  print("Step 2 completed") 
  
  
  #Computation of LCPs for ZPMs
  
  #Creation of empty object to store LCPs
  xarxa_camins_total <- st_sfc(st_linestring(), crs = crs(elevacio)) #creation of line string sfc object
  xarxa_camins_total <- st_sf(geometry = xarxa_camins_total) #creation of line string sf object
  xarxa_camins_total <- xarxa_camins_total[-1, ] #emptying of the line string sf object
  
  #Creation of LCPs 
    
    #Loop for LCP calculation
    for (i in 1:nrow(punts_loop2)) {
      
      #Recovery of coordinates of origin point
      punt_objectiu <- punts_loop2[i, ]
      
      #Recovery of the destination points
      punts_desti <- punts_loop2[-i, ]
      
      #LCPs calculation 
      xarxa_camins <- create_lcp(mde_base, punt_objectiu, punts_desti, cost_distance = FALSE) #lcps calcultaion
      xarxa_camins <- st_set_crs(xarxa_camins,crs(elevacio)) #seting crs
      
      #Sumation of all LCPs without attributes
      xarxa_camins_total <- rbind(xarxa_camins_total, xarxa_camins[ , -c(1,2,3,4)])
      
      rm(xarxa_camins) #saving memory
      
    }
  
  print("Step 3 completed") 
  
  #Saving of the LCPs to Shapefile
  nom_arxiu <- paste(substitute(y),substitute(q),substitute(w),"xarxa_camins_total.gpkg") #File name creation
  xarxa_camins_total <- xarxa_camins_total[!st_is_empty(xarxa_camins_total), ] #filter of empty geometries
  st_write(xarxa_camins_total, dsn = file.path(dir,nom_arxiu)) #Saving of LCPs object from sf object to GPKG
  
  print("Process ended succesfully")
  
}
FETE_vs_Edge_function <- function(x,x2,n,y,dir) {
  
  #Creating the folder to store the graphics
  dir.create(dir)
  
  # Loading the rasters
  raster_obj <- y
  
  raster_camins_1 <- x # Edge
  raster_camins_2 <- x2 # FETE
  
  #Creation of n equidistant control points
  
  #Poligonize raster objective
  pol_ras <- st_union(st_as_sf(as.polygons(raster_obj, dissolve = TRUE)))
  
  #Obtaining raster extension
  extent_raster <- ext(raster_obj) 
  
  #Identify and get Cells ID of n control points
  pc_y <- seq(from = extent_raster[3], to = extent_raster[4], length.out = n) #limit Y calculation
  pc_x <- seq(from = extent_raster[1], to = extent_raster[2], length.out = n) #limit X calculation
  xy_pc <- expand.grid(x = pc_x[2:(n-1)], y = pc_y[2:(n-1)]) #create a grid of control points
  pc_sf <- st_as_sf(xy_pc, coords = c("x","y"), crs = crs(raster_obj)) #control points to sf
  intersects_list <- st_intersects(pc_sf, pol_ras[[1]], sparse = TRUE) # logical for point within the raster
  pc_sf2 <- pc_sf[lengths(intersects_list) > 0, ] #extracting points
  
  #Extracting the values
  df_1 <- terra::extract(raster_camins_1, pc_sf2) #Extracting the values
  df_1[is.na(df_1)] <- 0 #Converting NA to 0
  df_2 <- terra::extract(raster_camins_2, pc_sf2) #Extracting the values
  df_2[is.na(df_2)] <- 0 #Converting NA to 0
  df_combinat <- cbind(df_1[ ,2],df_2[ ,2]) #storing the values
  colnames(df_combinat) <- c("Edge","FETE") #applying new name to column
  
  #Calculation of value
  Edge_vs_FETE <- df_combinat[ ,1] - df_combinat[ ,2]
  mat <- matrix(Edge_vs_FETE, nrow=7, byrow=TRUE)
  
  # Plooting the graphics #
  
  # Plot the matrix using the image function
  
  library(fields)
  colors <- colorRampPalette(c("blue", "white", "red"))(100) #Setting up color pallete
  
  tiff(file.path(dir,"heatmap.tiff"), width=7, height=7, units="in", res=300, compression="lzw")
  
  image.plot(1:7, 1:7, t(mat[nrow(mat):1, ]), col=colors, axes=FALSE, 
             zlim=c(-1,1), legend.mar=5)  
  box()
  
  dev.off()
  
  # Ploting Edge model with control points
  tiff(file.path(dir,"raster_edge_punts.tiff"), width=7, height=7, units="in", res=300, compression="lzw")
  
  plot(x) # Plotting raster
  plot(pc_sf2, add = TRUE, col = "White", bg = "white", pch = 21 ) #Plotting points
  
  dev.off()
  
  # Ploting FETE model with control points
  tiff(file.path(dir,"raster_FETE_punts.tiff"), width=7, height=7, units="in", res=300, compression="lzw")
  
  plot(x2)
  plot(pc_sf2, add = TRUE, col = "White", bg = "white", pch = 21 )
  
  dev.off()
  
}

Real_world_function <- function(y,y2,q,z,z2,w,dir) {
  
  #Creating the folders to store the results
  dir.create(dir)
  dir.create(file.path(dir, "SHAPES"))
  dir.create(file.path(dir, "RASTERS"))
  dir.create(file.path(dir, "RASTERS_RAW"))
  dir.create(file.path(dir, "RASTER_PARTS"))
  
  #Setting up functions to create KDE rasters for loop
  lcp_kernel <- function(lcp) {
    
    # Define the window (class owin) for the density analysis    
    win <- spatstat.geom::owin(c(ras_ex[1],ras_ex[3]),
                               c(ras_ex[2],ras_ex[4]))
    
    # Define ppp object based on point locations
    lcps_pts_pp <- spatstat.geom::ppp(x = sf::st_coordinates(lcp)[,1], 
                                      y = sf::st_coordinates(lcp)[,2], 
                                      window = win)
    
    # Calculate density 
    lcps_density <- spatstat.explore::density.ppp(
      x = lcps_pts_pp,
      sigma = 20 * res(elevacio)[1], 
      kernel = "gaussian",
      eps = 25 # cell size
    )
  }
  
  # Loading the raster to the function
  elevacio <- y 
  elevacio_back <- y2
  costa <- z2["geometry"]
  
  # Obtaining raster extent
  ext_r <- ext(elevacio)
  xmin <- ext_r[1]
  xmax <- ext_r[2]
  ymin <- ext_r[3]
  ymax <- ext_r[4]
  
  # Defining block size in cells
  block_cells <- 1250  # Number of cells per block in each dimension
  res_x <- res(elevacio)[1]  # Resolution in x direction
  res_y <- res(elevacio)[2]  # Resolution in y direction
  
  # Compute block size in spatial units
  block_size_x <- block_cells * res_x
  block_size_y <- block_cells * res_y
  
  # List to store the blocks
  blocks <- list()
  index <- 1
  
  # Divide the raster into blocks
  for (i in seq(xmin, xmax, by = block_size_x)) {
    for (j in seq(ymin, ymax, by = block_size_y)) {
      
      # Ensure last blocks fit within raster by adjusting size dynamically
      x_end <- min(i + block_size_x, xmax)
      y_end <- min(j + block_size_y, ymax)
      
      # Defining block extension
      ext_block <- ext(i, x_end, j, y_end)
      
      # Crop the raster block
      block_raster <- crop(elevacio, ext_block)
      
      # Write raster block to file
      is_empty <- all(is.na(values(block_raster))) # avoid saving empty rasters
      
      # Creation of the blocks
      if (is_empty) {
        
        print("Skipping raster: It is entirely NA")
        
      } else {
        
        blocks[[index]] <- block_raster
        writeRaster(blocks[[index]], filename = file.path(dir, "RASTER_PARTS", paste0("block_", index, ".tif")), overwrite = TRUE)
        
        print(paste("Block", index, "created with extension", ext_block))
        index <- index + 1
      }
    }
  }
  
  print("Step 1 completed")
  
  # Read all rasters
  for(r in seq_along(blocks)) {
    
    # Loading the raster's block
    elevacio_r <- blocks[[r]]
    
    #Identify valid part of the raster
    polygons_ras <- as.polygons(elevacio_r, dissolve = TRUE)
    names(polygons_ras)[1] <- "Layer_1" 
    polygons_ras_fil <- polygons_ras[polygons_ras$Layer_1 >= 0, ]
    polygons_ras_fil_sf <- st_as_sf(aggregate(polygons_ras_fil))
    
    #Creating buffer for paths creation and points locations
    buff_52 <- st_buffer(polygons_ras_fil_sf, 102 * res(elevacio)[1]) # for paths creation
    buff_50 <- st_buffer(polygons_ras_fil_sf, 100 * res(elevacio)[1]) #for points locations
    buff_ol <- st_buffer(polygons_ras_fil_sf, 50 * res(elevacio)[1]) #for points locations
    
    #Setting up the raster for path creation
    elevacio_r_52 <- crop(elevacio_back, buff_52)
    elevacio_r_52 <- ifel(is.na(elevacio_r_52), 0, elevacio_r_52)
    
    #Setting up points
    pol_line <- st_cast(buff_50, to = "MULTILINESTRING") #polygon to mutlilinestring
    pol_line <- st_cast(pol_line, to = "LINESTRING") #mutlilinestring to line 
    densitat <- 1 / q #calculation of point density depending on q
    p_ras <- st_line_sample(pol_line, density = densitat) #line to points depending of density of q
    punts_loop <- p_ras[[1]] #loading and division of edge points for loop functionality
    punts_loop <- st_sf(geometry = st_sfc(punts_loop)) #conversion of array to sf
    punts_loop <- st_cast(punts_loop, "POINT")
    
    #Application of the cost function (creation of the cost surface layer from package "leastcostpath")
    
    #Creation of base cs map 
    mde_base <- create_slope_cs(elevacio_r_52, cost_function = w, neighbours = 16, exaggeration = FALSE)
    
    #Creation of rivers' buffers based on "z"
    rius_r <- st_intersection(rius,polygons_ras_fil_sf)
    buff_rius <- st_buffer(rius_r, dist = 50) #buffer creation
    
    #Computing of the rivers' buffers cost multipliers 
    mde_base <- update_values(mde_base, buff_rius, FUN = function(mult) {mult / 1.5})
    
    #Computing of the rivers cost multipliers
    mde_base <- update_values(mde_base, rius_r, FUN = function(mult) {mult / 2})
    
    rm(rius_r)
    rm(buff_rius)
    
    #Extra slope cost multipliers calculation
    
    #Computing slope map
    pendent <- terrain(elevacio_r_52, v="slope", unit="radians", neighbors= 8) #Slope's raster calculation from "y"
    pendent <- tan(pendent) * 100 #Conversion of slope's raster values (radians for tan funcionality) to percentages
    
    #Computing of the > 60% slopes  cost multipliers  
    pendent_60 <- pendent #saving a copy of slope mde
    pendent_60[pendent > 60] <- 1 #identifying > 60% slopes
    pendent_60[pendent < 60] <- NA #erasing no > 60% slopes
    pendent60_pol <- st_as_sf(as.polygons(pendent_60, round =TRUE, aggregate= TRUE, values = FALSE)) #vectoring raster 
    mde_base <- update_values(mde_base, pendent60_pol, FUN = function(mult) { mult / 6}) #updating mde_base values
    
    #Computing of the > 100% slopes cost multipliers
    pendent_100 <- pendent
    pendent_100[pendent > 100] <- 1
    pendent_100[pendent < 100] <- NA
    pendent100_pol <- st_as_sf(as.polygons(pendent_100, round =TRUE, aggregate= TRUE, values = FALSE))
    mde_base <- update_values(mde_base, pendent100_pol, FUN = function(mult) { mult / 10})
    
    rm(pendent_100)
    rm(pendent_60)
    rm(pendent100_pol)
    rm(pendent60_pol)
    
    #Extra cost for wetlands
    # Logical condition: where both altitud < 5 and pendent < 5
    cond <- (elevacio_r_52 < 5) & (pendent < 5)
    # Create result raster: set to 1 where condition is TRUE, NA elsewhere
    altitud_10 <- ifel(cond, 1, NA)
    
    altitud_10_pol <- st_as_sf(as.polygons(altitud_10, round =TRUE, aggregate= TRUE, values = FALSE))
    mde_base <- update_values(mde_base, altitud_10_pol, FUN = function(mult) { mult / 1.5})
    
    rm(pendent)
    rm(altitud_10)
    rm(altitud_10_pol)
    
    gc()
    
    #Computation of LCPs for ZPMs
    
    #Creation of empty object to store LCPs
    xarxa_camins_total <- st_sfc(st_linestring(), crs = crs(elevacio)) #creation of line string sfc object
    xarxa_camins_total <- st_sf(geometry = xarxa_camins_total) #creation of line string sf object
    xarxa_camins_total <- xarxa_camins_total[-1, ] #emptying of the line string sf object
    
    #Creation of LCPs 
    
    #Loop for LCP calculation
    for (i in 1:nrow(punts_loop)) {
      
      #Recovery of coordinates of origin and end points
      punt_objectiu <- punts_loop[i, ]
      
      punts_desti <- punts_loop[-i, ]
      
      
      #LCPs calculation 
      xarxa_camins <- create_lcp(mde_base, punt_objectiu, punts_desti, cost_distance = FALSE) #lcps calcultaion
      xarxa_camins <- st_set_crs(xarxa_camins,crs(elevacio)) #seting crs
      
      #Sumation of all LCPs without attributes
      xarxa_camins_total <- rbind(xarxa_camins_total, xarxa_camins[ , -c(1,2,3,4)])
      
      rm(xarxa_camins) #saving memory
      
    }
    
    # Filtering paths with coast, buffers and erasing invalid geometries 
    xarxa_camins_total <- st_intersection(xarxa_camins_total,costa)
    xarxa_camins_total <- st_intersection(xarxa_camins_total,buff_ol)
    xarxa_camins_total %>% filter(!st_geometry_type(.) %in% c("POINT", "MULTIPOINT"))
    
    #Saving of the LCPs to Shapefile
    nom_arxiu <- paste("block",substitute(r),"xarxa_camins_total.gpkg")  #File name creation
    xarxa_camins_total <- xarxa_camins_total[!st_is_empty(xarxa_camins_total), ] #filter of empty geometries
    st_write(xarxa_camins_total, dsn = file.path(dir,"SHAPES",nom_arxiu)) #Saving of LCPs object from sf object to GPKG
    
    # Reducing raster to overlapping buffer extension
    elevacio_ol <- crop(elevacio_r_52, buff_ol)
    
    # Obtaining raster extent from r with overlaping areas
    ras_ex <- elevacio_ol |>
      sf::st_bbox()
    
    ## Compute density in spatstat (KDE) for each iteration 
    lcps_reduced <- sf::st_cast(xarxa_camins_total, "MULTIPOINT")
    
    lcp_kernel_data <- lcp_kernel(lcps_reduced)
    
    ## Save the density as GEOtiff for further analysis
    # convert it to a terra raster
    lcps_density_rast <- stars::st_as_stars(
      lcp_kernel_data) |> 
      sf::st_set_crs(25831) |>
      terra::rast()
    
    # Save
    terra::writeRaster(lcps_density_rast,filename = file.path(dir,"RASTERS_RAW",paste(as.character(r),"KDE_raster.tiff")))
    
    # Free memory
    rm(punt_objectiu)
    rm(punt_desti)
    rm(xarxa_camins_total)
    rm(lcps_density_rast)
    rm(lcps_reduced)
    rm(lcp_kernel_data)
    rm(buff_ol)
    rm(polygons_ras_fil_sf)
    rm(polygons_ras_fil)
    
    
    print(paste("Process ended succesfully in raster", substitute(r)))
    
  }
  
  # Normalize values of all the area #
  # List all TIFF files in the RASTERS_RAW folder
  llista_rast <- list.files(file.path(dir, "RASTERS_RAW"), pattern = "\\.tif[f]?$", full.names = TRUE)
  # Read all rasters
  rast_list <- lapply(llista_rast, terra::rast)
  
  # Identify MAX and MIN values
  #Setting up the function
  max_vals_set <- function(r) {
    
    max(values(r))  
    
  }
  
  min_vals_set <- function(r) {
    
    min(values(r))  
    
  }
  
  #Applying fuction to all the rasters
  max_vals <- sapply(rast_list,max_vals_set)
  min_vals <- sapply(rast_list,min_vals_set)
  
  #Normalizing values in regard MAX and MIN of the total area
  for(k in seq_along(rast_list)) { 
    
    k_rast <- rast_list[[k]]
    
    min_val <- min(min_vals)
    max_val <- max(max_vals)
    
    k_rast <- (k_rast - min_val) / (max_val - min_val) 
    
    terra::writeRaster(k_rast,filename = file.path(dir,"RASTERS",paste(as.character(k),"KDE_raster_norm.tiff"))) # Saving results
    
    rm(k_rast) 
    rm(vec_val) 
    
  }
  
  # Build the final mosaic #
  llista_rast_norm <- list.files(file.path(dir, "RASTERS"), pattern = "\\.tif[f]?$", full.names = TRUE)
  
  # Read all rasters
  rast_list_norm <- lapply(llista_rast_norm, terra::rast)
  rast_list_norm <- sprc(rast_list_norm)
  
  # Build the mosaic
  mosaic_rast <- mosaic(rast_list_norm, fun = mean)
  
  # Smooth borders between blocks
  mosaic_rast <- focal(mosaic_rast, w = matrix(1, nrow = 51, ncol = 51), fun = mean, na.policy = "omit")
  
  # Scaling results between 0 and 1
  min_val_v <- minmax(mosaic_rast)[1]
  max_val_v <- minmax(mosaic_rast)[2]
  
  mosaic_rast <- (mosaic_rast  - min_val_v) / (max_val_v - min_val_v) 
  
  # Save the output
  writeRaster(mosaic_rast, file.path(dir, "Real_case_mosaic.tiff"), overwrite = TRUE)
  
  print("Process endend succesfully") 
  
}

### Quantity series ###

# ---- Raster_A  ----

# Creation of the folder to store the models
dir.create("Quantity series")

# Loading of the raster and setting up the CRS
raster_a <- rast("Data/Rasters/Raster_A.tif")
crs(raster_a)<- "EPSG:25831"

# Creation of the folder to store the models
dir.create("Quantity series/Raster_a_Q")

# Loop for the creation of the models
distances <- seq(20, 1000, by = 20) # setting up the distances

time_results <- data.frame(time_in_seconds = numeric(0))

# Loop through distances
for (d in distances) {
  # Measure the time it takes to run the Edge_model
  time_taken <- system.time(Edge_model(raster_a, d, "tobler", "Quantity series/Raster_a_Q/SHAPES"))
  
  # Store the user time (you can use other time components like 'elapsed' or 'cpu' if needed)
  time_results <- rbind(time_results, data.frame(time_in_seconds = time_taken["elapsed"]))
  
  # Set the row name as the current value of d
  rownames(time_results)[nrow(time_results)] <- as.character(d)

}

# Kernel application in loop
llista_shp <- list.files("Quantity series/Raster_a_Q/SHAPES", pattern = "\\.gpkg$", full.names = TRUE) #Loading of the shapes models

for (i in llista_shp) {
  
  Kernel_function(i,raster_a,"Quantity series/Raster_a_Q/RASTERS")
  
}

# Comparison and graphics creation
Evaluation_function("Quantity series/Raster_a_Q/RASTERS", 9, raster_a,"Quantity series/Raster_a_Q/IMATGES")

rm(time_results)

# ---- Raster_B  ----

raster_b <- rast("Data/Rasters/Raster_B.tif")
crs(raster_b)<- "EPSG:25831"

dir.create("Quantity series/Raster_b_Q")

time_results <- data.frame(time_in_seconds = numeric(0))

# Loop through distances
for (d in distances) {
  # Measure the time it takes to run the Edge_model
  time_taken <- system.time(Edge_model(raster_b, d, "tobler", "Quantity series/Raster_b_Q/SHAPES"))
  
  # Store the user time (you can use other time components like 'elapsed' or 'cpu' if needed)
  time_results <- rbind(time_results, data.frame(time_in_seconds = time_taken["elapsed"]))
  
  # Set the row name as the current value of d
  rownames(time_results)[nrow(time_results)] <- as.character(d)
  
}

llista_shp <- list.files("Quantity series/Raster_b_Q/SHAPES", pattern = "\\.gpkg$", full.names = TRUE)

for (i in llista_shp) {
  
  Kernel_function(i,raster_b,"Quantity series/Raster_b_Q/RASTERS")
  
}

Evaluation_function("Quantity series/Raster_b_Q/RASTERS", 9, raster_b,"Quantity series/Raster_b_Q/IMATGES")

rm(time_results)

# ---- Raster_C  ----

raster_c <- rast("Data/Rasters/Raster_C.tif")
crs(raster_c)<- "EPSG:25831"

dir.create("Quantity series/Raster_c_Q")

time_results <- data.frame(time_in_seconds = numeric(0))

# Loop through distances
for (d in distances) {
  # Measure the time it takes to run the Edge_model
  time_taken <- system.time(Edge_model(raster_c, d, "tobler", "Quantity series/Raster_c_Q/SHAPES"))
  
  # Store the user time (you can use other time components like 'elapsed' or 'cpu' if needed)
  time_results <- rbind(time_results, data.frame(time_in_seconds = time_taken["elapsed"]))
  
  # Set the row name as the current value of d
  rownames(time_results)[nrow(time_results)] <- as.character(d)
  
}

llista_shp <- list.files("Quantity series/Raster_c_Q/SHAPES", pattern = "\\.gpkg$", full.names = TRUE)

for (i in llista_shp) {
  
  Kernel_function(i,raster_c,"Quantity series/Raster_c_Q/RASTERS")
  
}

Evaluation_function("Quantity series/Raster_c_Q/RASTERS", 9, raster_c,"Quantity series/Raster_c_Q/IMATGES")

rm(time_results)

### Spatial series ###

# Creation of the folder to store the models
dir.create("Spatial series")

# ---- Raster_A_A  ----

dir.create("Spatial series/Raster_a_S_a")

time_results <- data.frame(time_in_seconds = numeric(0))

# Creation of the rotate models
Edge_rotate_model(raster_a,50,5,"tobler","Spatial series/Raster_a_S_a/SHAPES")

llista_shp <- list.files("Spatial series/Raster_a_S_a/SHAPES", pattern = "\\.gpkg$", full.names = TRUE)

for (i in llista_shp) {
  
  Kernel_function(i,raster_a,"Spatial series/Raster_a_S_a/RASTERS")
  
}

Evaluation_function("Spatial series/Raster_a_S_a/RASTERS", 9, raster_a,"Spatial series/Raster_a_S_a/IMATGES")

rm(time_results)

# ---- Raster_A_B  ----

dir.create("Spatial series/Raster_a_S_b")

time_results <- data.frame(time_in_seconds = numeric(0))

Edge_rotate_model(raster_a,200,5,"tobler","Spatial series/Raster_a_S_b/SHAPES")

llista_shp <- list.files("Spatial series/Raster_a_S_b/SHAPES", pattern = "\\.gpkg$", full.names = TRUE)

for (i in llista_shp) {
  
  Kernel_function(i,raster_a,"Spatial series/Raster_a_S_b/RASTERS")
  
}

Evaluation_function("Spatial series/Raster_a_S_b/RASTERS", 9, raster_a,"Spatial series/Raster_a_S_b/IMATGES")

rm(time_results)

# ---- Raster_A_C  ----

dir.create("Spatial series/Raster_a_S_c")

time_results <- data.frame(time_in_seconds = numeric(0))

Edge_rotate_model(raster_a,500,5,"tobler","Spatial series/Raster_a_S_c/SHAPES")

llista_shp <- list.files("Spatial series/Raster_a_S_c/SHAPES", pattern = "\\.gpkg$", full.names = TRUE)

for (i in llista_shp) {
  
  Kernel_function(i,raster_a,"Spatial series/Raster_a_S_c/RASTERS")
  
}

Evaluation_function("Spatial series/Raster_a_S_c/RASTERS", 9, raster_a,"Spatial series/Raster_a_S_c/IMATGES")

rm(time_results)

### FETE vs Edge series ###

# Creation of the folder to store the models
dir.create("Edge vs FETE")

# ---- Edge_1,6_km  ----
raster_un_sis <- rast("Data/Rasters/Raster_un_sis_km.tif")
crs(raster_un_sis)<- "EPSG:25831"

dir.create("Edge vs FETE/Edge_un_sis_km")

Edge_model(raster_un_sis, 200, "tobler", "Edge vs FETE/Edge_un_sis_km/SHAPES")

llista_shp <- list.files("Edge vs FETE/Edge_un_sis_km/SHAPES", pattern = "\\.gpkg$", full.names = TRUE)

for (i in llista_shp) {
  
  Kernel_function(i,raster_un_sis,"Edge vs FETE/Edge_un_sis_km/RASTERS")
  
}

# ---- Edge_2_km  ----

raster_dos <- rast("Data/Rasters/Raster_dos_km.tif")
crs(raster_dos)<- "EPSG:25831"

dir.create("Edge vs FETE/Edge_dos_km")

Edge_model(raster_dos, 200, "tobler", "Edge vs FETE/Edge_dos_km/SHAPES")

llista_shp <- list.files("Edge vs FETE/Edge_dos_km/SHAPES", pattern = "\\.gpkg$", full.names = TRUE)

for (i in llista_shp) {
  
  Kernel_function(i,raster_dos,"Edge vs FETE/Edge_dos_km/RASTERS")
  
}

# ---- FETE_1_km  ----

dir.create("Edge vs FETE/FETE_un_km")

FETE_model(raster_a, 200, "tobler","Edge vs FETE/FETE_un_km/SHAPES")

llista_shp <- list.files("Edge vs FETE/FETE_un_km/SHAPES", pattern = "\\.gpkg$", full.names = TRUE)

for (i in llista_shp) {
  
  Kernel_function(i,raster_a,"Edge vs FETE/FETE_un_km/RASTERS")
  
}

# Loading of the desired raster for comparison
x <- rast("Edge vs FETE/Edge_un_km/RASTERS/raster_a 200 tobler xarxa_camins_total KDE_raster.tiff")
x2 <- rast("Edge vs FETE/FETE_un_km/RASTERS/raster_a 200 tobler xarxa_camins_total KDE_raster.tiff")

# Function to comparison and graphics creation
FETE_vs_Edge_function(x,x2,9, raster_a,"FETE_un_km")

# ---- FETE_1,6_km  ----

dir.create("Edge vs FETE/FETE_un_sis_km")

FETE_model(raster_un_sis, 200, "tobler","Edge vs FETE/FETE_un_sis_km/SHAPES")

llista_shp <- list.files("Edge vs FETE/FETE_un_sis_km/SHAPES", pattern = "\\.gpkg$", full.names = TRUE)

for (i in llista_shp) {
  
  Kernel_function(i,raster_un_sis,"Edge vs FETE/FETE_un_sis_km/RASTERS")
  
}

x <- rast("Edge vs FETE/Edge_un_sis_km/RASTERS/raster_un_sis_km 200 tobler xarxa_camins_total KDE_raster.tiff")
x2 <- rast("Edge vs FETE/FETE_un_sis_km/RASTERS/raster_un_sis_km 200 tobler xarxa_camins_total KDE_raster.tiff")

FETE_vs_Edge_function(x,x2,9, raster_un_sis,"FETE_un_sis_km")

# ---- FETE_2_km  ----

dir.create("Edge vs FETE/FETE_dos_km")

FETE_model(raster_dos, 200, "tobler","Edge vs FETE/FETE_dos_km/SHAPES")

llista_shp <- list.files("Edge vs FETE/FETE_dos_km/SHAPES", pattern = "\\.gpkg$", full.names = TRUE)

for (i in llista_shp) {
  
  Kernel_function(i,raster_dos,"Edge vs FETE/FETE_dos_km/RASTERS")
  
}

x <- rast("Edge vs FETE/Edge_dos_km/RASTERS/raster_dos 200 tobler xarxa_camins_total KDE_raster.tiff")
x2 <- rast("Edge vs FETE/FETE_dos_km/RASTERS/raster_dos 200 tobler xarxa_camins_total KDE_raster.tiff")

FETE_vs_Edge_function(x,x2,9, raster_dos,"FETE_dos_km")

# ---- Real world case ----

# Loading of the rasters and shapes necessary for the modelling
girona <- rast("Data/Rasters/girona.tif") # Raster of the case of study
girona_bg <- rast("Data/Rasters/girona_background.tif") # Background raster for overlapping areas
rius <- st_read("Data/Vectors/cat_ETRS89_31N_modificat.shp") # Shapes of the rivers for cost barriers
costa <- st_read("Data/Vectors/costa_girona_ETRS.shp") # Shapes of the cosatlines

# Application of the function to compute the real world case of study
Real_world_function(girona, girona_bg, 1250, rius, costa,"tobler","Real_case")

### Analysis of the real world case of study ### 

# Loading of the sites and natural corridors model
corridors <- rast("Real_case/Real_case_mosaic.tiff") # Natural corridors
jac <- st_read("Data/Vectors/punts_jaciments.shp") # Sites

# Obtaining point values
c_v_sites <- st_drop_geometry(jac[ ,2])
c_v_sites[ ,"Point value"] <- terra::extract(corridors, jac, ID = FALSE)

# Setting up focal window analysis
matrix_window <- matrix(1, nrow = 121, ncol = 121)

# Quantile 0.75
focal_q_75  <- focal(corridors, w = matrix_window, fun =function(x) quantile(x, 0.75, na.rm=TRUE), na.policy = "omit")

c_v_sites[ ,"25% higher values (3km buffer)"] <- terra::extract(focal_q_75, jac, ID = FALSE)

# Quantile 0.90
focal_q_90 <- focal(corridors, w = matrix_window, fun =function(x) quantile(x, 0.90, na.rm=TRUE), na.policy = "omit")

c_v_sites[ ,"10% higher values (3km buffer)"] <- terra::extract(focal_q_90, jac, ID = FALSE)

# Mean
focal_mean <- focal(corridors, w = matrix_window, fun = function(x) mean(x), na.policy = "omit")

c_v_sites[ ,"Mean value (3km buffer)"] <- terra::extract(focal_mean, jac, ID = FALSE)

# Elaborating the graphics #
# Reshape data to long format to combine to y values
df_long <- c_v_sites %>%
  pivot_longer(cols = c(2, 3, 4,5), names_to = "Variable", values_to = "Value")

# Plot the results
c_v_plot <- ggplot(df_long, aes(x = ID_jacimen, y = Value, color = Variable)) +
  geom_point(size = 1.5) +   # Points for each variable           # Line connecting points
  labs(title = "Sites and natural corridors: some metrics",
       x = NULL,
       y = "Natural corridors frequency (escalated 0 to 1)",
       color = "Metrics") +
  scale_color_brewer(palette = "Dark2") + 
  theme_classic() +
  theme(
    panel.grid.major.x = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.major.y = element_line(color = "gray80", linetype = "dashed"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Center title, make bold
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)# Add x-grid
  )

ggsave(file.path("Real_case","Corridors_values_girona.tiff"),c_v_plot, width = 12, height = 8) # save the file

# Normality test #
# Setting up function
normality_results <- lapply(c_v_sites, function(column) {
  if (is.numeric(column)) {  # Check if the column is numeric
    shapiro.test(column)$p.value  # Extract the p-value
  } else {
    NA  # If not numeric, return NA
  }
})

# Convert results to a data frame
normality_results_df <- data.frame(
  Column = names(normality_results),
  P_Value = unlist(normality_results)
)

# Add a column to interpret the results
normality_results_df$Normal <- ifelse(normality_results_df$P_Value > 0.05, "Yes", "No")

# Save results
write.csv(normality_results_df,file.path("Real_case","Normality_test.csv"),row.names = FALSE)

# Correlation test #
# Testing the data
cor_test <- cor(c_v_sites[ ,c(2,3,4,5)], method = "pearson") 

# Plotting and saving the data
tiff(file.path("Real_case","Correlation_plot_pearson.tiff"), width = 7*300, height = 4*300, res = 300) # Width and height in pixels
corrplot.mixed(cor_test,lower.col = "black", number.cex = 0.35, tl.pos="lt", tl.cex= 0.4, cl.cex = 0.4)
dev.off()


# ---- Supplementary: Raster_A_D  ----

dir.create("Spatial series_supp")

dir.create("Spatial series_supp/Raster_a_S_d")

time_results <- data.frame(time_in_seconds = numeric(0))

# Creation of the rotate models
Edge_rotate_model(raster_a,60,5,"tobler","Spatial series_supp/Raster_a_S_d/SHAPES")

llista_shp <- list.files("Spatial series_supp/Raster_a_S_d/SHAPES", pattern = "\\.gpkg$", full.names = TRUE)

for (i in llista_shp) {
  
  Kernel_function(i,raster_a,"Spatial series_supp/Raster_a_S_d/RASTERS")
  
}

Evaluation_function("Spatial series_supp/Raster_a_S_d/RASTERS", 9, raster_a,"Spatial series_supp/Raster_a_S_d/IMATGES")

rm(time_results)

# ---- Supplementary: Raster_A_E  ----

dir.create("Spatial series_supp/Raster_a_S_e")

time_results <- data.frame(time_in_seconds = numeric(0))

Edge_rotate_model(raster_a,75,5,"tobler","Spatial series_supp/Raster_a_S_e/SHAPES")

llista_shp <- list.files("Spatial series_supp/Raster_a_S_e/SHAPES", pattern = "\\.gpkg$", full.names = TRUE)

for (i in llista_shp) {
  
  Kernel_function(i,raster_a,"Spatial series_supp/Raster_a_S_e/RASTERS")
  
}

Evaluation_function("Spatial series_supp/Raster_a_S_e/RASTERS", 9, raster_a,"Spatial series_supp/Raster_a_S_e/IMATGES")

rm(time_results)

# ---- Supplementary: Raster_A_F  ----

dir.create("Spatial series_supp/Raster_a_S_f")

time_results <- data.frame(time_in_seconds = numeric(0))

Edge_rotate_model(raster_a,90,5,"tobler","Spatial series_supp/Raster_a_S_f/SHAPES")

llista_shp <- list.files("Spatial series_supp/Raster_a_S_f/SHAPES", pattern = "\\.gpkg$", full.names = TRUE)

for (i in llista_shp) {
  
  Kernel_function(i,raster_a,"Spatial series_supp/Raster_a_S_f/RASTERS")
  
}

Evaluation_function("Spatial series_supp/Raster_a_S_f/RASTERS", 9, raster_a,"Spatial series_supp/Raster_a_S_f/IMATGES")

rm(time_results)
