#!/usr/bin/env Rscript

#>>>>>>>>>>>>>Arvis Sulovari<<<<<<<<<<<<<<<<<
#>>>>>>>>>>>>>>>06/27/2016<<<<<<<<<<<<<<<<<<<

args = commandArgs(trailingOnly=TRUE)


#Error handling
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
}


#Required packages
require(data.table)



#Call in required scripts
source("VIP_detection_simulation.R")
source("GC_profile_data.R")
source("VIP_position_sim_function.R")
source("Read_Position_Simulation.R")
source("LTR_profile.R")


#Reference files (GC and LTR)
gc_profile <- read.delim("empirical_data/pIRS_GC_coverage_VIPS_200bpWindow.txt",header=T)

#vip_gc_profile <- read.delim("HBV_GC_profile_array.txt",header=F)
#colnames(vip_gc_profile) <- c("GC.")

vip_gc_profile_array <- read.delim("empirical_data/HBV_GC_profile_array.txt",header=T)
#vip_gc_profile_array$P_vip <- vip_gc_profile_array$P_vip/sum(vip_gc_profile_array$P_vip)

ltr_profile <- fread("empirical_data//RepeatMask_repeatRegions",header=T)
ltr_vip_distance <- read.table("empirical_data//ltr_vip_distance.txt",header=T)


#Save each trailing argument into the appropriate variable
human_virus_length <- as.numeric(args[1])
vip_nr <- as.numeric(args[2])
vip_len_mean <- as.numeric(args[3])
vip_len_sd <- as.numeric(args[4])
vip_len_min <- as.numeric(args[5])
seq_depth <- as.numeric(args[6])
read_length <- as.numeric(args[7])
read_insert_mean <-as.numeric(args[8])
read_insert_sd <- as.numeric(args[9])
read_perVIP <- as.numeric(args[10])
min_read_length <-as.numeric(args[11])
read_repeat_freq <- as.numeric(args[12])
read_trim_prop <- as.numeric(args[13])
read_trim <- as.numeric(args[14])
cell_prop <- as.numeric(args[15])
seed_value <- as.numeric(args[16])

#Start creating  profile files required during simulation
gc_array <- gc_profile_builder(human_virus_length,gc_profile,vip_gc_profile,seed_value)

vip_coords <- vip_sim(human_virus_length,vip_nr,vip_len_mean,vip_len_min,vip_len_sd,gc_array,seed_value)

ltr_array <- ltr_simulator(human_virus_length,ltr_vip_distance,vip_nr,vip_coords,ltr_profile)

sim_out <- array(0,dim=c(seq_depth*2,vip_nr*4))
final_out_arr <- array(NA, dim=c(seq_depth,2))


sim_out <- splitReadCounts(human_virus_length=human_virus_length,vip_nr=vip_nr,vip_coords=vip_coords,seq_depth=seq_depth,read_length=read_length,read_insert_mean=read_insert_mean,read_insert_sd=read_insert_sd,min_read_length=min_read_length,read_repeat_freq=read_repeat_freq,read_trim_prop=read_trim_prop,read_trim=read_trim,cell_prop=cell_prop,gc_array=gc_array,ltr_array = ltr_array,seed_value = seed_value)
    

#Count number of successfully identified virus integrations
success_counter <- 0
for(j in 1:vip_nr){
s <- na.omit(c(sim_out[2,((4*j)-3):((4*j)-2)]))
e <- na.omit(c(sim_out[2,((4*j)-1):(4*j)]))
#DEFINITION OF SUCCESS: More than read_perVIP supporting the VIP integration
ifelse((sum(s)+sum(e))>=read_perVIP,success_counter<-(success_counter+1),success_counter<-(success_counter+0))
}
#Separate reads into type
split_start <- c()
chimera_start <- c()
split_end <- c()
chimera_end <- c()
for(j in 1:vip_nr){
s_s <- na.omit(c(sim_out[2,((4*j)-3)]))
c_s <- na.omit(c(sim_out[2,((4*j)-2)]))
s_e <- na.omit(c(sim_out[2,((4*j)-1)]))
c_e <- na.omit(c(sim_out[2,(4*j)]))
#Collect the number of each read type at both ends of the VIP
split_start <- append(split_start,s_s)
chimera_start <- append(chimera_start,c_s)
split_end <- append(split_end,s_e)
chimera_end <- append(chimera_end,c_e)
}

#START POPULATING FINAL ARRAY
#Create
final_array <- array(NA, dim=c(1,26+(4*vip_nr)))
#Populate: Simulation parameters
final_array[1,1] <-human_virus_length
final_array[1,2] <- vip_nr
final_array[1,3] <- vip_len_mean 
final_array[1,4] <- vip_len_sd
final_array[1,5] <- vip_len_min
final_array[1,6] <- seq_depth
final_array[1,7] <- read_length
final_array[1,8] <- read_insert_mean
final_array[1,9] <- read_insert_sd
final_array[1,10] <- read_perVIP 
final_array[1,11] <- min_read_length
final_array[1,12] <- read_repeat_freq 
final_array[1,13] <- read_trim_prop 
final_array[1,14] <- read_trim
final_array[1,15] <- cell_prop
final_array[1,16] <- seed_value 
#Populate: Simulation results
final_array[1,17] <- (success_counter/vip_nr)*100
final_array[1,18] <- sum(split_start)
final_array[1,19] <- mean(split_start)
final_array[1,20] <- sum(chimera_start)
final_array[1,21] <- mean(chimera_start)
final_array[1,22] <- sum(split_end)
final_array[1,23] <- mean(split_end)
final_array[1,24] <- sum(chimera_end)
final_array[1,25] <- mean(chimera_end)
final_array[1,26] <- sys_time[1]
final_array[1,27:(26+4*vip_nr)] <- sim_out[2,1:(4*vip_nr)]


#Save simulation results into a csv file
write.csv(final_array,"VIPC_simulation_output.csv",append=T)