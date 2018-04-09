###FUNCTION START

vip_sim <- function(human_virus_length,
                    vip_nr,
                    vip_len_mean,
                    vip_len_min,
                    vip_len_sd,
                    gc_array,
                    seed_value=123456789)
  
  {
  
  
  #Probability of choosing a VIP start point is determined by GC content of the region
  set.seed(seed_value)
  vip_start <- sample(human_virus_length,vip_nr,prob = gc_array[,3],replace = F)
  
  #Use LTR information (Start & End) to postion VIPs accordingly
  #vip_start_tmp <- sample(human_virus_length,vip_nr,prob = gc_array[,3],replace = F)
  
  #Use LTR information (Start & End) to postion VIPs accordingly
  #vip_start
  
  #NEGATIVE BINOMIAL FORMULA
  #tmp_vip_start <- rnbinom(n = vip_nr,size = 1,prob = vip_nr/human_virus_length)
  #vip_start <- round(tmp_vip_start*human_virus_length/max(tmp_vip_start))
  
  vip_len_tmp <- c()
  
  l <- 0
  
  while(l<vip_nr) {
    
    #set.seed(seed_value)
    len <- round(abs(rnorm(1,vip_len_mean,vip_len_sd)))
    #len <- round(abs(rbeta(1,vip_len_mean,vip_len_sd)))
    
    ifelse(len >= vip_len_min,vip_len_tmp <- append(vip_len_tmp,len),l <- length(vip_len_tmp))
    
    l <- length(vip_len_tmp)
    #l <- 0
    
    #print(l)
    
  }
  
  
  vip_end <- vip_start+vip_len_tmp
  
  vip_coords <- cbind(vip_start,vip_end)
  vip_coords <- vip_coords[order(vip_coords[,1]),]
  
  return(vip_coords)
}



