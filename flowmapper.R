## Created by Rohan Shad, MD ##
## Maps extracted points from LVAD flow waveforms to generate patient specific flow curves ###

library(ggplot2)
library(spectral)
library(tidyverse)
library(readxl)
library(dplyr)
library(egg)


# Import plot from LVAD screen
#HR = 110 / Pulse Pressure is 9
HVAD.flow <- read.csv("flow5.csv", header = F, col.names = c("T","Flow"))

# Import pressure readings / other data (Dummy variables here)
lvad.set <- read.csv("patient_data_file.csv")
lvad.set$avg.sys <- read.csv("systolic_pressure_reads")
lvad.set$avg.dia <- read.csv("diastolic_pressure_reads")

# Calculate Pulse Pressure
lvad.set$avg.pulse <- lvad.set$avg.sys - lvad.set$avg.dia


## This section is essential to make sure the HVAD flow is periodic and smooth ##
HVAD.flow <- HVAD.flow[31:46,]
HVAD.flow[1,2] <- 5.108579
HVAD.flow[16,2] <- 5.108579
HVAD.flow[15,2] <- 5.108579
HVAD.flow[14,2] <- 5.444906

HVAD.flow$T <- HVAD.flow$T - HVAD.flow[1,1]

HVAD.splines <- as.data.frame(spline(HVAD.flow$T,HVAD.flow$Flow, method = 'natural'))
HVAD.splines_2 <- as.data.frame(spline(HVAD.flow$T,HVAD.flow$Flow, method = 'natural', xmin = 0.0))

# Plot out splines for quick check
g <- ggplot(data = HVAD.splines, mapping = aes(x = x, y = y)) + theme_article()
g <- g + geom_line(data = HVAD.splines, aes(x=x,y=y)) + labs(x = "Time in seconds", y = "Flow in L/min")
g <- g + coord_cartesian(xlim = c(0,4), ylim = c(-2,12))
g

#### Map patient specific curves ####

  # Patient ID, Heart Rate, LVAD flow output, Pulse Pressure required for mapping

  flow_mapper <- function(ID,HR,CO,Pulse){
    new_flow <- data.frame("X" = 120/HR * HVAD.flow[,1], "Y" = Pulse/20 * HVAD.flow[,2] + CO - mean(HVAD.flow[,2])*Pulse/20)
  
    # Creates flow splines and converts everything by 16.66 to get cgs units (For simvascular) from L/min flow
    new_flow_splines <- as.data.frame(spline(new_flow$X,new_flow$Y*16.66, xmin = 0))
    write.table(new_flow_splines, sep = " ", file = paste0("flow_files/", as.character(ID),"_test",".flow"), col.names = FALSE, row.names = FALSE)
  }


#### Function that feeds everything into flow mapper ####

  for (i in seq(lvad.set$ID)){
    flow_mapper(lvad.set[i,1], lvad.set[i,9],lvad.set[i,16],lvad.set[i,21])
  }

# Produces Clinical Data set #

  write.csv(lvad.set, file = "clinical_data.csv")


#### Function to create table of time periods for each case ####

  time_period <- function(ID,HR,CO,Pulse){
    new_flow <- data.frame("X" = 120/HR * HVAD.flow[,1], "Y" = Pulse/20 * HVAD.flow[,2] + CO - mean(HVAD.flow[,2])*Pulse/20)
    
    #Creates flow splines and converts everything by 16.66 to get cgs units from L/min flow
    new_flow_splines <- as.data.frame(spline(new_flow$X,new_flow$Y*16.66))
    new_flow_splines[48,1]
  }

  time.dat <- data.frame()

  for (i in seq(lvad.set$ID)){
    a <- time_period(lvad.set[i,1], lvad.set[i,9],lvad.set[i,16],lvad.set[i,21])
    dat <- data.frame(lvad.set[i,1],a)
    time.dat <- rbind(time.dat, dat)
  }

  time.dat$ID <- time.dat$lvad.set.i..1.
  time.dat$T <- time.dat$a

  write.csv(time.dat, file = "time_periods.csv", row.names = F)
