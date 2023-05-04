source("support_script/settings.R")
source("support_script/functions.R")

# in addition to files from fig3 we need the data from the seedling vs leaf model 
# and the state assignment of the extensive seedlingd model

##############################
## read in files and reformat
##############################

# states across genomes, function only keeps chr1-5 
# leaf
lcol = read_and_fix_states("data/seedling_leaf_15/Col_15_ChromHMM_ChromHMM_reduced_AT_segments.bed")
# seedling
scol = read_and_fix_states("data/seedling_leaf_15/Col_S_15_ChromHMM_ChromHMM_reduced_AT_segments.bed")

# emission matrix
raw = read_delim("data/seedling_leaf_15/emissions_15_ChromHMM_ChromHMM_reduced_AT.txt")
mat = raw_to_matrix(raw,sampleo) 

## reorder states
mat=mat[rord,sampleo]
rownames(mat)=NULL

# state assignment file from extensive model 
large_col=read_and_fix_states("data/seedling_26/AT_26_ChromHMM_ChromHMM_July_AT_segments.bed") # large WT 

# Change state names to fit paper (full model) using map_to_paper defined in sett_n.R
large_col=mutate(large_col,state=str_remove(state,"^E"))%>%left_join(map_to_paper,by=c("state"="from"))
large_col=select(large_col,-state)%>%mutate(state=to)

########################
## clean up 
########################
rm(list=c("map_mod2","map_to_paper","raw","sampleo","cs"))

