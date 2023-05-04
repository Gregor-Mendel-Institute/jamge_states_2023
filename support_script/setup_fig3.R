source("support_script/functions.R") # contains functions used across figures 
source("support_script/settings.R") # sets themes, defines the colors to be used etc 

##############################
## read in files and reformat
##############################

# states across genomes, function only keeps chr1-5 
col=read_and_fix_states("data/leaf_16/Col_16_ChromHMM_ChromHMM_concat_AT_segments.bed") #small WT
ddm=read_and_fix_states("data/leaf_16/ddm1_16_ChromHMM_ChromHMM_concat_AT_segments.bed") #small ddm

# use colors and state names from map (defined in settings.R)
col=left_join(col,map_mod2,by="state")%>%mutate(state=ns)%>%select(-ns)
ddm=left_join(ddm,map_mod2,by="state")%>%mutate(state=ns)%>%select(-ns)

# annotations
gff = get_our_gff()

# gene positions
gene_pos=get_gene_pos(gff)

# te gene positions
te_gene_pos=get_te_gene_pos(gff)

# extra tests: only chr 1-5
stopifnot(length(setdiff(gene_pos$chrom,c("1","2","3","4","5")))==0)
stopifnot(length(setdiff(te_gene_pos$chrom,c("1","2","3","4","5")))==0)

# add introns
gff_w_introns=add_intron(gff)

# select features of interest and rename them as we like to plot
feature=gff_to_feature(gff_w_introns)

# emission matrix small model
mat_raw=read_delim("data/leaf_16/emissions_16_ChromHMM_ChromHMM_concat_AT.txt")

# Make matrix with order of marks from emission tibble
mat = raw_to_matrix(mat_raw,sampleo) #sampleo defined (settings.R)
rownames(mat)=map_mod2$ns # the state names to use for small matrix (settings.R)

########################
## clean up 
########################
rm(list=c("map_mod2","map_to_paper","mat_raw","rord","sampleo","cs","gff","gff_w_introns"))
