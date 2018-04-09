#VACC jobs for parameter sweep

#Idea: save output from each simulation into an array (~20 columns)
##The array should have the following columns: every parameter being sweeped (~16, copied in two rows), ...
## ... number of split & chimeric reads in Start & End of the VIP (4 columns)


####################################STEP 0####################################
##LOAD previous WorkSpace file
#Fixed on May 16th, 2016
load("May16_WorkSpace.RData")

####################################STEP 1####################################



#Make each parameter an array of numbers

#Numeric arguments
human_virus_length <- c(100000,1000000,10000000,50000000)
vip_nr <- c(10,50,100)
vip_len_mean <- c(300)
vip_len_sd <- c(1000,10000)
vip_len_min <- c(0,50)
seq_depth <- c(1,2,3,4,5,6,7,8,9,10,20,30)
read_length <- c(50,75,100,120)
read_insert_mean <- 500
read_insert_sd <- read_insert_mean/20
read_perVIP <- 2
min_read_length <- 20
vip_repeat_freq <- c(0.2,0.80)
vip_clip_prop <- 0.2
#How many reads cover repeat regions (unaccounted for by RepeatMasker)
read_repeat_freq <- 0.01
read_trim_prop <- 0.2
read_level_qc_prop <- 0.8
read_trim <- 0.15
#seed_value <- sample(1:100000000,1)


#Create arrays for the array arguments
gc_array <- gc_profile_builder(human_virus_length,gc_profile,vip_gc_profile,seed_value)
vip_coords <- vip_sim(human_virus_length,vip_nr,vip_len_mean,vip_len_min,vip_len_sd,gc_array,seed_value)
ltr_array <- ltr_simulator(human_virus_length,ltr_vip_distance,vip_nr,vip_coords,ltr_profile)

#Array arguments
#vip_coords <- vip_coords
#gc_aray <- gc_array
#ltr_array <- ltr_array

#####################################SINGLE INSTANCE#####################################

sim_out <- (
  splitReadCounts(human_virus_length=human_virus_length,
                            
                  #vip_nr=100,
                  #vip_len_mean=200,
                  #vip_len_sd=20,
                  #vip_len_min=50,
                  vip_coords <- vip_coords,
                  #Assume that only 5% of cells (from tissue, blood, etc) have VIP integrations
                  #cell_prop=0.05,
                            
                  seq_depth=seq_depth,
                  read_length=read_length,
                  read_insert_mean=read_insert_mean,
                  read_insert_sd=read_insert_mean/20,
                  read_perVIP=read_perVIP,
                            
                  #Minimum length of a sequence that can be aligned (default BWA Q-value, PHRED scaled)
                  min_read_length=min_read_length,
                  #How frequently are VIPs clipped?
                  vip_repeat_freq=vip_repeat_freq,
                  #By how much are VIPs clipped (proportion of VIP average)?
                  vip_clip_prop=vip_clip_prop,
                  #How many reads are zero'ed?
                  read_repeat_freq=read_repeat_freq,
                  #How much are reads are trimmed?
                  read_trim_prop=read_trim_prop,
                  #How much of a read should have high quality nucleotides
                  read_level_qc_prop=read_level_qc_prop,
                  #By how much are the reads trimmed (proportion of read length)?
                  read_trim=read_trim,
                  #mismatch_rate=0.04,
                            
                  gc_aray <- gc_array,
                  #vip_gc_adj_factor=5,
                  
                  ltr_array <- ltr_array,
                  seed_value = sample(1:100000000,1)
)
)


#sim_out[((2*i)-1):(2*i),1:(vip_nr*4)] <- out_tmp
#colnames(sim_out) <- colnames(out_tmp)


#Count number of successfully identified
success_counter <- 0

for(j in 1:vip_nr){
  
  s <- na.omit(c(sim_out[(2*i),((4*j)-3):((4*j)-2)]))
  e <- na.omit(c(sim_out[(2*i),((4*j)-1):(4*j)])) 
  
  
  #DEFINITION OF SUCCESS: More than two reads supporting the VIP integration
  ifelse((sum(s)+sum(e))>=read_perVIP,success_counter<-(success_counter+1),success_counter<-(success_counter))
  
}

#Separate reads into type
split_start <- c()
chimera_start <- c()
split_end <- c()
chimera_end <- c()

for(j in 1:vip_nr){
  
  s_s <- na.omit(c(sim_out[(2*i),((4*j)-3)]))
  c_s <- na.omit(c(sim_out[(2*i),((4*j)-2)]))
  
  s_e <- na.omit(c(sim_out[(2*i),((4*j)-1)]))
  c_e <- na.omit(c(sim_out[(2*i),(4*j)]))
  
  #Collect the number of each read type at both ends of the VIP
  split_start <- append(split_start,s_s)
  chimera_start <- append(chimera_start,c_s)
  split_end <- append(split_end,s_e)
  chimera_end <- append(chimera_end,c_e)
  
}



#START POPULATING FINAL ARRAY
#Create
final_array <- array(NA, dim=c(1,27))

#Populate: Simulation parameters
final_array[1,1] <- human_virus_length
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
final_array[1,12] <- vip_repeat_freq 
final_array[1,13] <- vip_clip_prop 
final_array[1,14] <- read_repeat_freq 
final_array[1,15] <- read_trim_prop 
final_array[1,16] <- read_level_qc_prop 
final_array[1,17] <- read_trim 
final_array[1,18] <- seed_value 

#Populate: Simulation results
final_array[1,19] <- success_counter
final_array[1,20] <- sum(split_start)
final_array[1,21] <- mean(split_start)
final_array[1,22] <- sum(chimera_start)
final_array[1,23] <- mean(chimera_start)
final_array[1,24] <- sum(split_end)
final_array[1,25] <- mean(split_end)
final_array[1,26] <- sum(chimera_end)
final_array[1,27] <- mean(chimera_end)


#WRITE FINAL_ARRAY *AND* APPEND OUTPUT FOR PARALLEL WRITING
write.table(final_array,"FINAL_ARRAY_TABLE.txt",append = T)






#####################################CODE WRITING#####################################
#Write R code to be run - each permutation of parameters per line

#May 19th Update: 
#1) Replaced "seq_depth" with "read_nr" in VIP_detection_simulation_TEST



human_virus_length <- c(1000000)
vip_nr <- 50
vip_len_mean <- 500
vip_len_sd <- vip_len_mean*10
vip_len_min <- c(10,50,200)
seq_depth <- c(1,2,4,6,8,10,20,40)
#read_nr <- c(5000,10000,20000,30000,40000,50000,100000)
read_length <- c(75,100,120,300)
read_insert_mean <- c(300,1000,2000)
read_insert_sd <- read_insert_mean/20
read_perVIP <- c(2,4,6,8,10)
min_read_length <- c(20,40)
#vip_repeat_freq <- c(0.2,0.5)
#vip_clip_prop <- c(0.1,0.5)
read_repeat_freq <- 0.05
read_trim_prop <- 0.05
#read_level_qc_prop <- c(0.05,0.1)
read_trim <- 0.15
cell_prop <- c(0.01,0.1,0.2,1)
#seed_value <- sample(1:100000000,1)



for(human_virus_length_0 in human_virus_length){
  
  for(vip_nr_0 in vip_nr){
    
    for(vip_len_mean_0 in vip_len_mean){
      
  #    for(vip_len_sd_0 in vip_len_sd){
        
        #Out of 3
        for(vip_len_min_0 in vip_len_min[3]){
          
          #for(read_nr_0 in read_nr){
            
            for(seq_depth_0 in seq_depth){
            
              for(read_length_0 in read_length){
              
                for(read_insert_mean_0 in read_insert_mean){
                
                  for(min_read_length_0 in min_read_length){
                  #Out of 5
                   for(read_perVIP_0 in read_perVIP[5]){
                  
                     for(read_trim_prop_0 in read_trim_prop){
                    
                       for(read_trim_0 in read_trim){
                      
                         for(cell_prop_0 in cell_prop){
                
                                  code <- paste0("
                                                 sys_time <- system.time(
                                                 for(z in 1){
                                                 seed_value <- sample(1:100000000,1)
                                                 #gc_array <- gc_profile_builder(",human_virus_length_0,",gc_profile,vip_gc_profile,seed_value)
                                                 vip_coords <- vip_sim(",human_virus_length_0,",",vip_nr_0,",",vip_len_mean_0,",",vip_len_min_0,",",vip_len_mean_0*10,",gc_array,seed_value)
                                                 ltr_array <- ltr_simulator(",human_virus_length_0,",ltr_vip_distance,",vip_nr_0,",vip_coords,ltr_profile)
                                                 sim_out <- (splitReadCounts_test(human_virus_length=",human_virus_length_0,",vip_nr=",vip_nr_0,",vip_coords = vip_coords, seq_depth=",seq_depth_0,",read_length=",read_length_0,",read_insert_mean=",read_insert_mean_0,",read_insert_sd=",read_insert_mean_0/20,",read_perVIP=",read_perVIP_0,",min_read_length=",min_read_length_0,",read_repeat_freq=",0.05,",read_trim_prop=",read_trim_prop_0,",read_trim=",read_trim_0,",cell_prop=",cell_prop_0,",gc_array=gc_array,ltr_array=ltr_array,seed_value = seed_value))
                                                 #Count number of successfully identified
                                                 success_counter <- 0
                                                 for(j in 1:",vip_nr_0,"){
                                                 s <- na.omit(c(sim_out[2,((4*j)-3):((4*j)-2)]))
                                                 e <- na.omit(c(sim_out[2,((4*j)-1):(4*j)]))
                                                 #DEFINITION OF SUCCESS: More than read_perVIP supporting the VIP integration
                                                 ifelse((sum(s)+sum(e))>=",read_perVIP_0,",success_counter<-(success_counter+1),success_counter<-(success_counter+0))
                                                 }
                                                 #Separate reads into type
                                                 split_start <- c()
                                                 chimera_start <- c()
                                                 split_end <- c()
                                                 chimera_end <- c()
                                                 for(j in 1:",vip_nr_0,"){
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
                                                 }
                                                 )
                                                 #START POPULATING FINAL ARRAY
                                                 #Create
                                                 final_array <- array(NA, dim=c(1,26+4*",vip_nr_0,"))
                                                 #Populate: Simulation parameters
                                                 final_array[1,1] <-",human_virus_length_0,"
                                                 final_array[1,2] <- ",vip_nr_0,"
                                                 final_array[1,3] <- ",vip_len_mean_0," 
                                                 final_array[1,4] <- ",vip_len_mean_0*10,"
                                                 final_array[1,5] <- ",vip_len_min_0,"
                                                 final_array[1,6] <- ",seq_depth_0,"
                                                 final_array[1,7] <- ",read_length_0,"
                                                 final_array[1,8] <- ",read_insert_mean_0,"
                                                 final_array[1,9] <- ",read_insert_mean_0/20,"
                                                 final_array[1,10] <- ",read_perVIP_0," 
                                                 final_array[1,11] <- ",min_read_length_0,"
                                                 final_array[1,12] <- 0.05 
                                                 final_array[1,13] <- ",read_trim_prop_0," 
                                                 final_array[1,14] <- ",read_trim_0,"
                                                 final_array[1,15] <- ",cell_prop_0,"
                                                 final_array[1,16] <- seed_value 
                                                 #Populate: Simulation results
                                                 final_array[1,17] <- success_counter
                                                 final_array[1,18] <- sum(split_start)
                                                 final_array[1,19] <- mean(split_start)
                                                 final_array[1,20] <- sum(chimera_start)
                                                 final_array[1,21] <- mean(chimera_start)
                                                 final_array[1,22] <- sum(split_end)
                                                 final_array[1,23] <- mean(split_end)
                                                 final_array[1,24] <- sum(chimera_end)
                                                 final_array[1,25] <- mean(chimera_end)
                                                 final_array[1,26] <- sys_time[1]
                                                 final_array[1,27:(26+4*",vip_nr_0,")] <- sim_out[2,1:(4*",vip_nr_0,")]
                                                 #WRITE FINAL_ARRAY *AND* APPEND OUTPUT FOR PARALLEL WRITING
                                                 write.table(final_array,\"ReadNR_SIMULATION_OUTPUT_PART15_of_15.txt\",col.names = F, row.names = F,append = T)"
                                                 )
                                  
                                  cat(code)
                                  
                       }
                       
                     }
                     
                   }
                   
                 }
                  
                }
                
              }
              
            }
            
          }
          
        }
        
      }
  }
}





#s <- paste0(a," and ",b)
#print(s)
