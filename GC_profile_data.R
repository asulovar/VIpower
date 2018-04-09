#LOAD NECESSARY GC-related FILES
# 
# #Coverage - GC file
# gc_profile <- read.delim("pIRS_GC_coverage_VIPS_200bpWindow.txt",header=T)
# 
# #VIP position - GC file (raw - many rows)
# vip_gc_profile <- read.delim("HBV_GC_profile.txt",header=F)
# colnames(vip_gc_profile) <- c("GC.")
# 
# #VIP position - GC file (table)
# vip_gc_profile_array <- read.delim("HBV_GC_profile_array.txt",header=T)
# vip_gc_profile_array$P_vip <- vip_gc_profile_array$P_vip/sum(vip_gc_profile_array$P_vip)

#Builds a GC-profile for a given human-virus sequence
#Output has three columns:
#1) GC percentage (for 200bp tiling windows)
#2) Probability of a read being present - GC matched (#reads in a given 200bp region/#total reads)
#3) Probability of a VIP integration - GC adjusted (#VIP integrations in a given 200bp region/#total reads)

gc_profile_builder <- function(human_virus_length,gc_profile,vip_gc_profile,seed_value=1234567){


#############
######1######
#############
#Create Array where each nucleotide has a %GC value (first col)


gc_array <- array(NA, dim=c(human_virus_length,3))
colnames(gc_array) <- c("GC.","P_coverage_tmp","P_VIP")

#FIXED: Created array where each window (200bp length) has a %GC value (first col)
#Populate with %GC values, by sampling from the pool of %GC

#set.seed(seed_value)

#Number of tiling windows across the simulated sequence (200bp tiling windows -- same as input data)
gw_tiling_windows <- round(human_virus_length/200)

#Pool of GC-levels, sampled using the genome-wide GC distribution
s <- (sample(gc_profile[-c(1),1],gw_tiling_windows,replace=T,prob = gc_profile[-c(1),]$RefCnt/sum(gc_profile[-c(1),]$RefCnt)))

#GC levels - each pool sampled above is repeated 200 times     
#gc_array[,1] <- ((rep(sort(s),200)))

pool <- c()
window_order <- sample(length(s),length(s),replace=F)
  
#for (i in 1:length(s)) {
for (i in window_order) {
  
#  tmp <- ((rep((s[i]),200)))
  tmp <- ((rep((s[i]),200)))
  pool <- append(pool,c(tmp))

}

gc_array[,1] <- pool



#gc_array[,1] <- (sample(gc_profile[,1],human_virus_length,replace=T,prob = (gc_profile[,2])/sum(gc_profile[,2])))


#############
######2######
#############
#Add column of Coverage Probability (second col)

#SUM OF ALL COVERAGES...
gc_profile$P_coverage <- gc_profile$DepthCnt/sum(gc_profile$DepthCnt)

gc_profile_v2 <- merge(gc_array,gc_profile,by="GC.",sort = F)

#gc_profile_v2$P_coverage <- (gc_profile_v2$DepthCnt)/sum((gc_profile_v2$DepthCnt))

#MAX COVERAGE...
#gc_profile_v2$P_coverage <- gc_profile_v2$SmoothedMean/max(gc_profile_v2$SmoothedMean)


#############
######3######
#############
#Add column of VIP integrations' probability (third col)


gc_profile_v3 <- merge(gc_profile_v2,vip_gc_profile_array,by="GC.",sort = F)


gc_array_v4 <- cbind(gc_profile_v3$GC.,gc_profile_v3$P_coverage,gc_profile_v3$P_vip.x)

#colnames(gc_array_v4) <- c("GC","P_coverage","P_vip")

return(gc_array_v4)

}




#Histograms

#plot(gc_profile$GC.,gc_profile$SmoothedMean,xlab="%GC (200bp window)",ylab="Coverage",type="b")

#plot(gc_profile$GC.,gc_profile$RefCnt/sum(gc_profile$RefCnt),xlab="%GC (200bp window)",ylab="Counts(Genome_wide)",type="b")
