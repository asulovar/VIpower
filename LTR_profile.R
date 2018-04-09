#LOAD the necessary LTR-related files
#require(data.table)
#ltr_profile <- fread("LTR_folder/RepeatMask_repeatRegions",header=T)

#Distribute simulated LTR regions across the simulated sequence

ltr_simulator <- function(human_virus_length,ltr_vip_distance,vip_nr,vip_coords,ltr_profile,seed_value){
  
  #Need VIP coords so that we can place LTR regions around them using *Empirical data*
  
  #vip_coords <- vip_sim(human_virus_length,vip_nr,vip_len_mean,vip_len_sd,gc_array,seed_value=seed_value)
  
  
  #Estimate number of LTRs using genome-wide LTR data
  #Using a fixed (upper limit) for the genome length
  gw_length <- 3234000000
  gw_ltr_nr <- nrow(ltr_profile)
  gw_ltr_lengths <- ltr_profile$genoEnd-ltr_profile$genoStart
  
  genome_prop <- human_virus_length/gw_length
  ltr_nr <- round(genome_prop*gw_ltr_nr)
  
  #Simulate LTR start positions (without replacement)
  #ltr_start <- sample(1:human_virus_length,ltr_nr,replace = F)
  
  ltr_start <- c()
  ltr_end <- c()
  
  #*#*#*#*#*#*#*#LTR-adjustment#*#*#*#*#*#*#*#
  #Use LTR-VIP distance information (distibutions) to postion LTRs accordingly
  
  #Step 1: create a position for LTRs based on VIP's Start & End coordinates
  ##Step 1.1: Distance between VIP_start and LTR_end should be drown from the ltr_vip_distance[,2]
  ltr_end_vip_start_dist_index <- sample(1:nrow(ltr_vip_distance),vip_nr,replace = F)
  ltr_end_vip_start_dist <- ltr_vip_distance[ltr_end_vip_start_dist_index,2]
  
  ###Assign END & START coordinates to LTRs upstream of the VIPs
  ltr_end_up <- vip_coords[,1]-ltr_end_vip_start_dist
  ltr_lengths_up <- sample(gw_ltr_lengths,vip_nr,replace = F)
  ltr_start_up <- ltr_end_up-ltr_lengths_up
  
  #Add to the final count array
  ltr_start <- append(ltr_start,ltr_start_up)
  ltr_end <- append(ltr_end,ltr_end_up)
    
  ##Step 1.2: Distance between VIP_end and LTR_start should be drown from ltr_vip_distance[,1]
  ltr_start_vip_end_dist_index <- sample(1:nrow(ltr_vip_distance),vip_nr,replace = F)
  ltr_start_vip_end_dist <- ltr_vip_distance[ltr_start_vip_end_dist_index,1]
  
  ###Assign START & END coordinates to LTRs downstream of the VIPs
  ltr_start_down <- vip_coords[,2]+ltr_start_vip_end_dist
  ltr_lengths_down <- sample(gw_ltr_lengths,vip_nr,replace = F)
  ltr_end_down <- ltr_start_down+ltr_lengths_down
  
  #Add to the final count array
  ltr_start <- append(ltr_start,ltr_start_down)
  ltr_end <- append(ltr_end,ltr_end_down)
  
  #Simulate remainder of the LTR positions using empirical LTR length distribution;
  #No. of rest of LTRs is: (ltr_nr - 2*vip_nr)
  rest_ltr <- ltr_nr - (length(ltr_start))
  
  if(rest_ltr>0) {
         #The rest of LTRs will be sampled from the rest of genome that has no LTRs yet, to *avoid overlap*
         ltr_start_rest <- sample(1:human_virus_length,rest_ltr,replace = F)
         ltr_lengths_rest <- sample(gw_ltr_lengths,rest_ltr,replace = F)
         ltr_end_rest <- ltr_start_rest+ltr_lengths_rest
         
         ltr_start <- append(ltr_start,ltr_start_rest)
         ltr_end <- append(ltr_end,ltr_end_rest)
  
  } else {
    ltr_out <- cbind(ltr_start,ltr_end)
    #print(ltr_out)
    #break
  }

    #Done
  ltr_out <- cbind(ltr_start,ltr_end)
  
  return(ltr_out)
  
}

