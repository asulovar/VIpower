
##############FINAL FUNCTION##################

splitReadCounts <- function(human_virus_length,
                            vip_nr,
                            #vip_len_mean,
                            #vip_len_sd,
                            #vip_len_min,
                            vip_coords,
                            seq_depth,
                            #read_nr,
                            read_length,
                            read_insert_mean,
                            read_insert_sd,
                            read_perVIP,
                            min_read_length,
                            #vip_repeat_freq,
                            #vip_clip_prop,
                            read_repeat_freq,
                            #repeat_freq,
                            read_trim,
                            read_trim_prop,
                            #read_level_qc_prop,
                            gc_array,
                            ltr_array,
                            cell_prop,
                            seed_value=123456789) {
  
  
  
  #seq_depth <- ((read_nr*2)*read_length)/human_virus_length
  
  #Calculate read_nr, so that read_coords can use it in the next step
  read_nr <- round((seq_depth*human_virus_length)/(read_length*2))
  
  #gc_array <- gc_profile_builder(human_virus_length,gc_profile,vip_gc_profile,seed_value)
  
  #Simulate VIP coordinates
  #vip_coords <- vip_sim(human_virus_length,vip_nr,vip_len_mean,vip_len_sd,gc_array,seed_value)
  
  #Simulate read coordinates
  read_coords <- read_simulator(human_virus_length,read_nr,read_length,read_insert_mean,read_insert_sd,gc_array,seed_value)
  
  
  
  ######################
  ######QC Block########
  ######################
  
  #How many reads to trim?
  #set.seed(seed_value*4)
  reads_to_clip <- sample(read_nr,read_trim_prop*read_nr)
  
  #How many reads to zero?
  #set.seed(seed_value*5)
  reads_to_zero <- sample(read_nr,read_repeat_freq*read_nr)
  
  #How many VIPs to clip?
  #set.seed(seed_value*6)
  #vips_to_clip <- sample(vip_nr,vip_repeat_freq*vip_nr)
  
  
  
  #*****************Nucleotide-level QC (Reads)*****************
  
  #[1] Trim end of the reads (both pairs) based on Q30 (BWA)
  for (i in 1:length(reads_to_clip)) {
    
    read_coords[reads_to_clip[i],] <- read_coords[reads_to_clip[i],]-(c(0,read_length*read_trim,0,read_length*read_trim))
    
  }
  
  
  
  #*******************Read-level QC*********************
  #MOST (>50%) of the failed reads are found in the regions where VIP are integrated!
  # 
  
  for (i in 1:length(reads_to_zero)) {
    
    col <- sample(c(1,3),1)
    
    read_coords[reads_to_zero[i],col:(col+1)] <- NA
    
  }
  
  
  
  #READ-level QC - final step: Remove "PCR duplications"
  read_coords <- unique(read_coords)
  read_nr <- nrow(read_coords)
  
  
  
  #**********************TEST FOR OVERLAP OF Virus insertions with PE reads 
  #Test for overlap between VIP start/end and read 1 start/end or read2 start/end
  vip_out_arr <- array(0,dim=c(2,(vip_nr*4)))
  
  #split_reads_vec_VIPstart <- c()
  #split_reads_vec_VIPend <- c()
  
  #chimeric_reads_vec_VIPstart <- c()
  #chimeric_reads_vec_VIPend <- c()
  
  for (i in 1:vip_nr) {
    
    split_reads_vec_VIPstart <- c()
    split_reads_vec_VIPend <- c()
    
    chimeric_reads_vec_VIPstart <- c()
    chimeric_reads_vec_VIPend <- c()
    
   # for (j in 1:read_nr) {
      
      #     split_reads <- (vip_coords[i,1] >= read_coords[j,1] && vip_coords[i,1] <= read_coords[j,2]) || 
      #     (vip_coords[i,1] >= read_coords[j,3] && vip_coords[i,1] <= read_coords[j,4]) ||
      #     (vip_coords[i,2] >= read_coords[j,1] && vip_coords[i,2] <= read_coords[j,2]) ||
      #     (vip_coords[i,2] >= read_coords[j,3] && vip_coords[i,2] <= read_coords[j,4])
      #     
      
      #######################
      ###Split-Reads Block###
      #######################
    
      #Count number of split-reads that fall on the VIP baundaries (STart and End)
      
      
        ##VIP Start
        ####Read 1 (end): 
        virus_start_read1_end <- which((((read_coords[,2] - vip_coords[i,1]) >= min_read_length) & 
                                          ((read_coords[,2] - vip_coords[i,1]) <= (read_length-min_read_length)))=="TRUE")
           
        ####Read 2 (end)
        virus_start_read2_end <- which((((read_coords[,4]-vip_coords[i,1]) >= min_read_length) & 
                                          ((read_coords[,4]-vip_coords[i,1]) <= (read_length-min_read_length)))=="TRUE")
        
        
        #####Check if the START of the supporting split READ 1 on VIP start overlaps with repeat region
        #####Count how many of these supporting READ 1 need to be removed
        bad_reads_virus_start_read1_end <- c()
        bad_reads_virus_start_read2_end <- c()
          
        #FIND BAD SPLIT READS: inside LTRs and < min_mappable region on human
        
        if(length(virus_start_read1_end)>0) {
          
            for(j in 1:length(virus_start_read1_end)) {
      
              ######Make sure that Read 1 (start) is not inside an LTR = NOT(INSIDE_LTR)
              #m <- length(which(!(((read_coords[virus_start_read1_end[j],1] >= ltr_array[,1]) & (read_coords[virus_start_read1_end[j],1] <= ltr_array[,2]))=="TRUE")))
  
              m <- length(which(
                
                (
                  
                  (
                    (ltr_array[,1] <= (read_coords[virus_start_read1_end[j],1]-min_read_length)) & 
                      
                      (ltr_array[,2] >= read_coords[virus_start_read1_end[j],1]) &
                        #Gap between LTR_end and vip_start is less than min_read_length
                        (ltr_array[,2] > (vip_coords[i,1]-min_read_length))
                  )
                ) =="TRUE"
              )
              )
      
              ifelse(m>0,bad_reads_virus_start_read1_end <- append(bad_reads_virus_start_read1_end,1),
                     bad_reads_virus_start_read1_end <- append(bad_reads_virus_start_read1_end,0)
                     )

            }
        }
        
  
        
        if(length(virus_start_read2_end)>0) {
          
          for(j in 1:length(virus_start_read2_end)) {
            
            #FIND BAD READS: inside LTRs and < min_mappable region on human
            
            n <- length(which(
              
            (
              
              (
                (ltr_array[,1] <= (read_coords[virus_start_read2_end[j],3]-min_read_length)) & 
                  
                  (ltr_array[,2] >= read_coords[virus_start_read2_end[j],3]) &
                  
                    #Gap between LTR_end and vip_start is less than min_read_length
                    (ltr_array[,2] > (vip_coords[i,1]-min_read_length))
                
                )
              ) =="TRUE"
                             )
                        )
            
            ifelse(n>0,bad_reads_virus_start_read2_end <- append(bad_reads_virus_start_read2_end,1),
                   bad_reads_virus_start_read2_end <- append(bad_reads_virus_start_read2_end,0))
            
           }
        }
        
        
        
        #Add up number of CORRECTED supporting split reads
        split_reads_nr_VIPstart <-   length(c(virus_start_read1_end,virus_start_read2_end)) - 
                                       sum(c(bad_reads_virus_start_read1_end,bad_reads_virus_start_read2_end))
        
        ###CELLULARITY ADJUSTMENT
        split_reads_nr_VIPstart <- floor(split_reads_nr_VIPstart*cell_prop)
        ########
        


        ##VIP End
        virus_end_read1_end <- c()
        virus_end_read2_end <- c()
        
        
        ####Read 1 (end)
        virus_end_read1_end <- which(
                                      (
                                        ((read_coords[,2]-vip_coords[i,2]) >= min_read_length) & 
                                       ((read_coords[,2]-vip_coords[i,2]) <= (read_length-min_read_length))
                                      ) =="TRUE"
                                   )


        ####Read 2 (end)
        virus_end_read2_end <- which(
                                      (
                                        ((read_coords[,4]-vip_coords[i,2]) >= min_read_length) & 
                                      ((read_coords[,4]-vip_coords[i,2]) <= (read_length-min_read_length))
                                     ) == "TRUE"
                                   )
      
        
        #####Check if the START of the supporting split READ 1 on VIP end overlaps with repeat region
        #####Count how many of these supporting READ 1 need to be removed (i.e., bad_reads_*)
        bad_reads_virus_end_read1_end <- c()
        
        ###Make sure that Read 1 (end) is not inside an LTR = NOT(INSIDE_LTR)
        #***
        if(length(virus_end_read1_end)>0) {
          
          for(j in 1:length(virus_end_read1_end)) {
          
            m <- length(which(
              
              (
                
                (
                  (ltr_array[,1] <= read_coords[virus_end_read1_end[j],2]) & 
                    
                    (ltr_array[,2] >= (read_coords[virus_end_read1_end[j],2]-min_read_length)) &
                  
                      #Gap between vip_end and LTR_start is less than min_read_length
                      (ltr_array[,1] <= (vip_coords[i,2]+min_read_length))
                )
              ) =="TRUE"
            )
            )
            
            ifelse(m>0,bad_reads_virus_end_read1_end <- append(bad_reads_virus_end_read1_end,1),
                   bad_reads_virus_end_read1_end <- append(bad_reads_virus_end_read1_end,0))
            
          }
        }
        
        
        
        bad_reads_virus_end_read2_end <- c()
        ###Make sure that Read 2 (end) is not inside an LTR = NOT(INSIDE_LTR)
        #***
        if(length(virus_end_read2_end)>0) {
          
          for(j in 1:length(virus_end_read2_end)) {
            
            ######Make sure that Read 2 (END) is not inside an LTR = NOT(INSIDE_LTR)
            #m <- length(which(!(((read_coords[virus_start_read1_end[j],1] >= ltr_array[,1]) & (read_coords[virus_start_read1_end[j],1] <= ltr_array[,2]))=="TRUE")))
            
            n <- length(which(
              
              (
                
                (
                  (ltr_array[,1] <= read_coords[virus_end_read2_end[j],4]) & 
                    
                    (ltr_array[,2] >= read_coords[virus_end_read2_end[j],4]-min_read_length) &
                    
                    #Gap between vip_end and LTR_start is less than min_read_length
                    (ltr_array[,1] <= (vip_coords[i,2]+min_read_length))
                )
              ) =="TRUE"
            )
            )
            
            ifelse(n>0,bad_reads_virus_end_read2_end <- append(bad_reads_virus_end_read2_end,1),
                   bad_reads_virus_end_read2_end <- append(bad_reads_virus_end_read2_end,0))
            
          }
        }
        
        
        ###Make sure that Read 2 (end) is not inside an LTR = NOT(INSIDE_LTR)
        #!(read_coords[j,4] >= ltr_array[,1] & read_coords[j,4] <= ltr_array[,2])
      
        #Add up number of CORRECTED supporting split reads
        split_reads_nr_VIPend <- length(c(virus_end_read1_end,virus_end_read2_end)) - 
          sum(c(bad_reads_virus_end_read1_end,bad_reads_virus_end_read2_end))
        
        
        ###CELLULARITY Adjustment
        split_reads_nr_VIPend <- floor(split_reads_nr_VIPend*cell_prop)
        ########
        
      ###TODO: Remove?
      #Append all non-zero split-read counts
      split_reads_vec_VIPstart <- append(split_reads_vec_VIPstart,c(split_reads_nr_VIPstart))
      split_reads_vec_VIPend <- append(split_reads_vec_VIPend,c(split_reads_nr_VIPend))
      
      
      
      ##########################
      ###Chimeric Reads Block### 
      ##########################
      
      #IMPORTANT: To make sure we do not include split reads into chimeric reads
      #we will only count those reads found entirely inside a virus insertion region.
      
      #####Chimeric at VIP Start
      
      #READ 2 is mapped entirely inside the virus insertion region 
      #AND read 1 is not inside the virus insertion region
      #Only unique reads are kept in case read 1 is one insertion 
      #and read 2 is in another
      
      chimeric_reads_location_VIPstart <-which(
        (
        #Read 2 START vs VIP start
        (read_coords[,3] - vip_coords[i,1]) > (-min_read_length) &
          # Read 2 END vs VIP END
          (read_coords[,4] - vip_coords[i,2]) < min_read_length &
            #Read 1 END is not inside the virus insertion
            (read_coords[,2] - vip_coords[i,1]) < min_read_length

        )=="TRUE"
      )

      chimeric_reads_location_VIPstart <- unique(chimeric_reads_location_VIPstart)
      
      
      #Placeholder for COUNT of bad chimeric reads on Viruses' start position
      bad_chimeric_reads_VIPstart <- c()
      
      #IDEA: bad chimeric read = The *OTHER* read (i.e. READ 1) of the pair is 'inside' an LTR
      if(length(chimeric_reads_location_VIPstart)>0) {
        
        for(j in 1:length(chimeric_reads_location_VIPstart)) {
          
          ######Make sure that Read 1 (END) is not inside an LTR = NOT(INSIDE_LTR)
          #m <- length(which(!(((read_coords[virus_start_read1_end[j],1] >= ltr_array[,1]) & (read_coords[virus_start_read1_end[j],1] <= ltr_array[,2]))=="TRUE")))
          
          m <- length(which(
            
            (
              (
                #END of Read 1 is inside an LTR
                #(read_coords[chimeric_reads_location_VIPstart[j],2] >= ltr_array[,1]) & 
                  
                  #(read_coords[chimeric_reads_location_VIPstart[j],2] <= ltr_array[,2]) &
                    
                
                    #START of Read 1 is "inside" an LTR *AND* END of Read 1 is is inside the LTR
                    #LTR_start-Read_1_start < (min_map) 'AND' Read_1_end-LTR_end < (min_map)
                    #(
                      (ltr_array[,1] - read_coords[chimeric_reads_location_VIPstart[j],1]) < min_read_length &
                       (read_coords[chimeric_reads_location_VIPstart[j],2]-ltr_array[,2]) < min_read_length
                     #)
                  
              )
            ) =="TRUE"
          )
          )
          
          ifelse(m>0,bad_chimeric_reads_VIPstart <- append(bad_chimeric_reads_VIPstart,1),
                 bad_chimeric_reads_VIPstart <- append(bad_chimeric_reads_VIPstart,0))
          
        }
      }
      
      #***
      ##ADJUST FOR THE BAD READS
      chimeric_reads_nr_VIPstart <- length(chimeric_reads_location_VIPstart)-sum(bad_chimeric_reads_VIPstart)
      
      
      ###CELLULARITY
      chimeric_reads_nr_VIPstart <- floor(chimeric_reads_nr_VIPstart*cell_prop)
      ########
      
      
      
      
      ##Chimeric at VIP End
      
      chimeric_reads_location_VIPend <- which(
        (
          ####Read 1 END vs VIP END
          (read_coords[,2] - vip_coords[i,2]) < min_read_length &
            ####Read 1 START vs VIP START
            (read_coords[,1] - vip_coords[i,1]) > (-min_read_length) &
              ####Read 2 START is NOT inside the insertion
              (read_coords[,3] - vip_coords[i,2]) > min_read_length
        )=="TRUE"
      )
      
      
      #Placeholder for COUNT of bad chimeric reads on Viruses' end position
      bad_chimeric_reads_VIPend <- c()
      
      #***
      if(length(chimeric_reads_location_VIPend)>0) {
        
        for(j in 1:length(chimeric_reads_location_VIPend)) {
          
          ######Make sure that ###Read 2 (end) is not inside an LTR
          
          n <- length(which(
            
            (
              
              (
                #(read_coords[chimeric_reads_location_VIPend[j],4] >= ltr_array[,1]) & 
                  
                  #(read_coords[chimeric_reads_location_VIPend[j],4] <= ltr_array[,2]) &
                  
                  #START of Read 2 is 'inside' LTR Start *AND* END of Read 2 is 'inside' LTR End
                  #LTR_start-Read_2_start < (min_map) 'AND' Read_2_end-LTR_end < (min_map)
                  #(
                    (ltr_array[,1] - read_coords[chimeric_reads_location_VIPstart[j],3]) < min_read_length &
                      (read_coords[chimeric_reads_location_VIPstart[j],4]-ltr_array[,2]) < min_read_length
                  #)
              )
            ) =="TRUE"
          )
          )
          
          ifelse(n>0,bad_chimeric_reads_VIPend <- append(bad_chimeric_reads_VIPend,1),
                 bad_chimeric_reads_VIPend <- append(bad_chimeric_reads_VIPend,0))
          
        }
      }
      
      
      
      ##ADJUST FOR THE BAD READS
      chimeric_reads_nr_VIPend <- length(chimeric_reads_location_VIPend)-sum(bad_chimeric_reads_VIPend)
      
      ###CELLULARITY
      chimeric_reads_nr_VIPend <- floor(chimeric_reads_nr_VIPend*cell_prop)
      ########
      
      #Append all non-zero chimeric read counts
      chimeric_reads_vec_VIPstart <- append(chimeric_reads_vec_VIPstart,c(chimeric_reads_nr_VIPstart))
      chimeric_reads_vec_VIPend <- append(chimeric_reads_vec_VIPend,c(chimeric_reads_nr_VIPend))
      
      
      #*******###END OF CHIMERIC READS*******#
      
      
    #}
      ###END OF for(j in read_nr) loop
    
    #Add number of reads that hit the VIP boundries (start and end) 
    #0, 1 and 2 times (found in table[1], table[2] and table[3])
    
    #0 READS
    #+vip_out_arr[1,(2*i)-1]; +vip_out_arr[1,2*i]
    vip_out_arr[1,(4*i)-3] <- table(split_reads_vec_VIPstart)[1]
    vip_out_arr[1,(4*i)-2] <- table(chimeric_reads_vec_VIPstart)[1]
    
    vip_out_arr[1,(4*i)-1] <- table(split_reads_vec_VIPend)[1]
    vip_out_arr[1,4*i] <- table(chimeric_reads_vec_VIPend)[1]
    
    #1 READ
    #+vip_out_arr[2,(2*i)-1]); +vip_out_arr[2,2*i])
    vip_out_arr[2,(4*i)-3] <- split_reads_nr_VIPstart
    vip_out_arr[2,(4*i)-2] <- chimeric_reads_nr_VIPstart
    
    vip_out_arr[2,(4*i)-1] <- split_reads_nr_VIPend
    vip_out_arr[2,4*i] <- chimeric_reads_nr_VIPend
    
    #2 READS: SANITY CEHCK - IMPOSSIBLE FOR TWO READS OF THE SAME PAIR TO BE IN THE SAME PLACE!
    #+vip_out_arr[3,(2*i)-1]; +vip_out_arr[3,2*i]
    #vip_out_arr[3,(2*i)-1] <- table(split_reads_vec_VIPstart)[3]
    #vip_out_arr[3,2*i] <- table(split_reads_vec_VIPend)[3]
    
  
  #colnames(vip_out_arr) <- rep(c("split-reads_VIP_start","chimeric-reads_VIP_start","split-reads_VIP_end",
  #                               "chimeric-reads_VIP_end"),vip_nr)
    
  }
  
  return(vip_out_arr)
  
}


