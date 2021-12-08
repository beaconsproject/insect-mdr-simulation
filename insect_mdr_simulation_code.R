# Spruce budworm simulation model demo using fictional data

library(raster)
library(foreign)
library(dplyr)
library(tidyr)
library(data.table)

# INPUT DATA #########################################################################

# Prepare reserve raster, bF raster, and bF age raster
# Here we use a fictional square ecoregion of 10,000 km2 and square reserve of 100km2. raster cell sizes are 1km2 and all rasters are snapped to the same grid.

# Make reserve map. Must have value 1 for the reserve, 0 for the study region, and NA outside the study region.
res <- raster(nrows = 102, ncols = 102, vals = 1)
res[c(2:10, 21:101),] <- 0
res[,c(2:10, 21:101)] <- 0
res[c(1, 102),] <- NA
res[,c(1, 102)] <- NA

# Create random bF map where 25% of ecoregion is bF cells
bF <- res
bF[!is.na(bF)] <- sample (c(NA,1), size=length(bF[!is.na(bF)]), replace=T, prob=c(.75,.25))

# randomly assign age classes 1:7 to bF cells. Age classes are 1: 0-20, 2: 20-40, 3: 40-60, 4: 60-80, 5: 80-100, 6: 100-120, 7: 120+.
bFage <- bF
bFage[!is.na(bFage)] <- sample (c(1:7), size=length(bFage[!is.na(bFage)]), replace=T)
cell_area <- 1


# PROCESS INPUT DATA #################################################################

# Convert bF map into stand map where every cell has a unique value
bFstand <- bF
values(bFstand)[values(bFstand)==1 & !is.na(values(bFstand))] <- 1:length(values(bFstand)[values(bFstand)==1 & !is.na(values(bFstand))])

# Make reserve stand map and calculate reserve stand areas. This is used later to calculate area of each age class in the reserve.
res_df <- as.data.frame(freq(res)) # used to calc reserve area
resStand <- res * bFstand
resStand[resStand < 1] <- NA # convert 0's to NA

# Convert reserve stand map to data frame
res_stand_df <- as.data.frame(freq(resStand)) # convert to data frame of cells
res_stand_df$Area_km2_reserve <- res_stand_df$count * cell_area # add area column
names(res_stand_df)[names(res_stand_df)=="value"] <- "stand_id" # rename
res_stand_df <- res_stand_df[!res_stand_df$stand_id %in% 0,] # remove zero row if it exists
res_stand_df <- res_stand_df[!is.na(res_stand_df$stand_id),] # remove NA row if it exists

# Convert full landscape stand map to data frame and assign ages
landscape_df <- as.data.frame(stack(bFstand, bFage)) # Combine age and bF stand rasters into a stack and convert to df
names(landscape_df) <- c("stand_id", "Age_orgnl") # rename
landscape_df <- landscape_df[!landscape_df$stand_id %in% 0 & !landscape_df$Age_orgnl %in% 0,] # remove zero row if it exists
landscape_df <- landscape_df[!is.na(landscape_df$stand_id) & !is.na(landscape_df$Age_orgnl),] # remove NA row if it exists
landscape_df$Area_km2 <- cell_area # add area column
landscape_df$age_at_year_0 <- landscape_df$Age_orgnl

# Calculate important areas
region_area <- sum(res_df$count[!is.na(res_df$value)]) * cell_area # ecoregion area
total_stand_area <- sum(landscape_df$Area_km2) # bF area in ecoregion
reserve_stand_area <- sum(res_stand_df$Area_km2_reserve) # bF area in reserve
reserve_area <- res_df$count[res_df$value==1 & !is.na(res_df$value)] * cell_area # reserve area

print(paste0("Ecoregion area: ", region_area))
print(paste0("bF area in ecoregion: ", total_stand_area))
print(paste0("Reserve area: ", reserve_area))
print(paste0("bF area in reserve: ", reserve_stand_area))


# MODEL PARAMETERS ###################################################################

reps <- 100 # how many times to repeat the simulation
run_length <- 260 # how long the simulation runs in years
ob_freq <- 40 # set frequency of outbreaks in years


# RUN MODEL ##########################################################################

# make stand tables directory
stand_tables_Dir <- paste0("./stand_tables/")
dir.create(stand_tables_Dir, recursive=TRUE)

# In order to alter outbreak frequencies in the model we need to update age every time step and trigger outbreaks based on requested frequency
landscape_dt <- as.data.table(landscape_df)
ob_years <- seq(20, run_length, ob_freq) # make vector of years that will have outbreaks

for(r in 1:reps){
  if(r==1 | r%%10==0){print(r)}
  for(yr in seq(20, run_length, 20)){ # at each 20 year time step...
    
    #print(yr)
    
    # add new column
    landscape_dt$old_age <- landscape_dt[[paste0("age_at_year_",yr-20)]]
    landscape_dt$new_age <- landscape_dt$old_age + 1
    
    if(yr %in% ob_years){ # if current year is an outbreak year, trigger outbreak. For outbreak stands, over ride the new age and return to age class 1.
      
      # set probability a cell experiencing outbreak dies. Randomly select between 0.78 and 0.89 based on mortality estimates from Ostaff and MacLean and MacLean 1980.
      landscape_dt$p_mortality <- sample(c(seq(0.78, 0.89, 0.01)), size = nrow(landscape_dt), replace = TRUE)
      landscape_dt$r_1 <- runif(nrow(landscape_dt), 0, 1) # generate a uniform random number
      
      # if the stand is over age 60, and the random number is less than the probability of mortality, the cell dies and age returns to 0. Using 60 to represent mature stands.
      # setting to zero instead of 1 means that growth will not begin until the next time step. This simulates the fact tht regrowth doesn't start until the outbreak is over. We are assuming the outbreak lasts 20 years (i.e. 1 time step).
      landscape_dt[new_age > 4 & r_1 < p_mortality, new_age := 0]
      
      # delete P columns
      landscape_dt[,p_mortality:=NULL]
      landscape_dt[,r_1:=NULL]
    }
    landscape_dt[[paste0("age_at_year_",yr)]] <- landscape_dt$new_age # name new age column
    landscape_dt[,old_age:=NULL] # delete temp columns
    landscape_dt[,new_age:=NULL]
  }
  write.csv(as.data.frame(landscape_dt), paste0(stand_tables_Dir, "stand_table_",r,".csv"), row.names = FALSE)
}
      

##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################



# Summarise data
# Summary table - proportion of ecoregion in each age class
# Calculate min area and min proportion, max area and max proportion, mean area and mean proportion
# Repeat for the reserve
counter <- 1
for(r in 1:reps){
  repDf <- read.csv(paste0(stand_tables_Dir, "stand_table_",r,".csv")) # open df
  
  #############################################
  # aggregate into age classes 1-40, 41-80, 80+
  for(col in colnames(repDf)[which(grepl("age_at_year", names(repDf)))]){
    repDf[[col]][repDf[[col]] %in% c(0,1,2)] <- 1 # 1-40
    repDf[[col]][repDf[[col]] %in% c(3,4)] <- 2 # 41-80
    repDf[[col]][repDf[[col]]>4] <- 3 # 80+
  }
  #############################################
  
  # FULL LANDSCAPE - make long format and select columns
  repDf_1 <- repDf %>%
    gather(timestep, age_class, which(grepl("age_at_year",names(repDf)))) %>%
    group_by(timestep, age_class) %>%
    summarise(Area2=sum(Area_km2)) %>%
    mutate(repetition=r, proportion_stand=Area2/total_stand_area, proportion_region=Area2/region_area) %>%
    as.data.frame()
  
  # merge
  if(counter == 1){
    landscapeDf <- repDf_1
  } else{
    landscapeDf <- rbind(landscapeDf, repDf_1)
  }
  
  # MDR RESERVE - make long format and select columns
  # proporions reported as proportion of total stand area and proportion of reserve area
  repDf_2a <- repDf %>%
    filter(stand_id %in% res_stand_df$stand_id) # subset by stands in reserve
  repDf_2b <- merge(repDf_2a, res_stand_df[c("stand_id","Area_km2_reserve")], by="stand_id")
  repDf_2 <- repDf_2b %>%
    gather(timestep, age_class, which(grepl("age_at_year",names(repDf)))) %>%
    group_by(timestep, age_class) %>%
    summarise(Area2=sum(Area_km2_reserve)) %>%
    mutate(repetition=r, proportion_stand=Area2/reserve_stand_area, proportion_reserve=Area2/reserve_area) %>%
    as.data.frame()
  
  # merge
  if(counter == 1){
    reserveDf <- repDf_2
  } else{
    reserveDf <- rbind(reserveDf, repDf_2)
  }
  counter <- counter + 1
}

# ISSUE - if age_class falls to zero, it doesn't get reported in output table. Search through all timesteps and reps and add in any missing rows with 0's
for(t in unique(landscapeDf$timestep)){
  for(r in unique(landscapeDf$repetition)){
    test_list <- landscapeDf$age_class[landscapeDf$timestep==t & landscapeDf$repetition==r]
    for(i in 1:3){
      if(!i %in% test_list){
        landscapeDf <- rbind(landscapeDf, data.frame(timestep=t, age_class=i, Area2=0, repetition=r, proportion_stand=0, proportion_reserve=0))
      }
    }
  }
}
for(t in unique(reserveDf$timestep)){
  for(r in unique(reserveDf$repetition)){
    test_list <- reserveDf$age_class[reserveDf$timestep==t & reserveDf$repetition==r]
    for(i in 1:3){
      if(!i %in% test_list){
        reserveDf <- rbind(reserveDf, data.frame(timestep=t, age_class=i, Area2=0, repetition=r, proportion_stand=0, proportion_reserve=0))
      }
    }
  }
}

# make landscape summary
summaryDf1 <- landscapeDf %>%
  group_by(age_class) %>%
  summarise(min_area_landscape=round(min(Area2),2),
            min_prop_landscape=round(min(proportion_region),3), 
            max_area_landscape=round(max(Area2),2), 
            max_prop_landscape=round(max(proportion_region),3), 
            mean_area_landscape=round(mean(Area2),2), 
            mean_prop_landscape=round(mean(proportion_region),3)) %>%
  as.data.frame()

# Make reserve summary
summaryDf2 <- reserveDf %>%
  group_by(age_class) %>%
  summarise(min_area_reserve=round(min(Area2),2), 
            min_prop_reserve=round(min(proportion_reserve),3), 
            max_area_reserve=round(max(Area2),2), 
            max_prop_reserve=round(max(proportion_reserve),3), 
            mean_area_reserve=round(mean(Area2),2), 
            stdev_area_reserve=round(sd(Area2),2), 
            mean_prop_reserve=round(mean(proportion_reserve),3)) %>%
  as.data.frame()

summaryDf <- merge(summaryDf1, summaryDf2, by="age_class")
print(summaryDf)