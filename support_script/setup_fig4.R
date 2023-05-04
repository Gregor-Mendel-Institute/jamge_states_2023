source("support_script/functions.R")
source("support_script/settings.R")

########################
## read in files
########################

# states across genomes chr1-5
col=read_and_fix_states("data/leaf_16/Col_16_ChromHMM_ChromHMM_concat_AT_segments.bed") #small WT
ddm=read_and_fix_states("data/leaf_16/ddm1_16_ChromHMM_ChromHMM_concat_AT_segments.bed") #small ddm

########################
## fix names
########################
# use colors and state from map (defined in settings.R)
col=left_join(col,map_mod2,by="state")%>%mutate(state=ns)%>%select(-ns)
ddm=left_join(ddm,map_mod2,by="state")%>%mutate(state=ns)%>%select(-ns)

########################
## some checks:
########################
stopifnot(length(setdiff(col$chrom,c("1","2","3","4","5")))==0)
stopifnot(length(setdiff(ddm$chrom,c("1","2","3","4","5")))==0)

########################
## figure specific 
## settings, colors etc
######################## 
qc=brewer.pal(9,"YlGnBu")[1:5]

########################
## clean up 
########################
rm(list=c("map_mod2","map_to_paper","m_order","rord","sampleo","cs"))

