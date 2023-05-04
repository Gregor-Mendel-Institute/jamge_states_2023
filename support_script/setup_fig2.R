source("support_script/functions.R")
source("support_script/settings.R")

##############################
## read in files and reformat
#############################
# states across genome, function only keeps chr1-5 
large_col=read_and_fix_states("data/seedling_26/AT_26_ChromHMM_ChromHMM_July_AT_segments.bed") # large WT 
 
# Change state names to fit paper (full model) using map_to_paper defined in settings.R
large_col=mutate(large_col,state=str_remove(state,"^E"))%>%
    left_join(map_to_paper,by=c("state"="from"))
large_col=select(large_col,-state)%>%mutate(state=to)
large_col=large_col%>%mutate(state=factor(state,levels = c("H1","H2","H3","H4","H5","H6",
                                                           "F1","F2","F3","F4","F5","F6" ,
                                                           "I1","I2","I3",
                                                           "E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11" )))

# annotations
gff = get_our_gff()
# add introns
gff_w_introns=add_intron(gff)
# genes
gene_pos=get_gene_pos(gff)

# expression data
load("data/expression/tpm_seedlings.Rdata")
tpm_seedlings=left_join(gene_pos,tpm_seedlings)

## methylation, mnase and dnase data
cg = read_and_fix_states("data/metyl_dnase_mnase/S_26_CG_weightedmean.bed")
chg=read_and_fix_states("data/metyl_dnase_mnase/S_26_CHG_weightedmean.bed")
chh=read_and_fix_states("data/metyl_dnase_mnase/S_26_CHH_weightedmean.bed")
mnase=read_and_fix_states("data/metyl_dnase_mnase/S_26_mnaseSig.bed")
dnase=read_and_fix_states("data/metyl_dnase_mnase/S_26_dnaseSig.bed")

# Change state names to fit paper (full model) using map_to_paper defined in sett_n.R, change . to NA
cg=mutate(cg,state=str_remove(state,"^E"))%>%left_join(map_to_paper,by=c("state"="from"))%>%
  mutate(mean=ifelse(X7==".",NA,X7))%>%mutate(mean=as.numeric(mean))
chg=mutate(chg,state=str_remove(state,"^E"))%>%left_join(map_to_paper,by=c("state"="from"))%>%
  mutate(mean=ifelse(X7==".",NA,X7))%>%mutate(mean=as.numeric(mean))
chh=mutate(chh,state=str_remove(state,"^E"))%>%left_join(map_to_paper,by=c("state"="from"))%>%
  mutate(mean=ifelse(X7==".",NA,X7))%>%mutate(mean=as.numeric(mean))
mnase=mutate(mnase,state=str_remove(state,"^E"))%>%left_join(map_to_paper,by=c("state"="from")) # no . in this file
dnase=mutate(dnase,state=str_remove(state,"^E"))%>%left_join(map_to_paper,by=c("state"="from"))%>%
  mutate(value=ifelse(X7==".",NA,X7))%>%mutate(value=as.numeric(value))

## Jaccard Index files (JI is the diagonal of the confusion matrix)
files =c("allMarks.txt","noH2B.txt","noH3.txt","noH2A.txt","noH3mod.txt","noVariants.txt","onlySeq_Mendes.txt")
mat=list()
jindex=matrix(NA,ncol=length(files),nrow=nrow(map_to_paper))
for (i in 1:length(files)){
  file=files[i]
  con_mat1_raw=read_delim(paste0("data/EvalSubset/",file),skip=1)
  mat[[i]] = as.matrix(con_mat1_raw[,-1])
  rownames(mat[[i]])=str_remove(con_mat1_raw$...1,"^[A-Z]")
  colnames(mat[[i]])=str_remove(colnames(mat[[i]]),"^[A-Z]")
  mat[[i]] = mat[[i]][map_to_paper$from,map_to_paper$from]
  colnames(mat[[i]])=map_to_paper$to
  rownames(mat[[i]])=map_to_paper$to
  jindex[,i]=diag(mat[[i]])
}
colnames(jindex)=c("full dataset","no H2B","no H3", "no H2A", "no H3 PTMs", "no H2A/H3/H2B","2014")
rownames(jindex)=rownames(mat[[i]])  

########################
## clean up 
########################
rm(list=c("map_mod2","map_to_paper","m_order","rord","sampleo","cs2","mat","con_mat1_raw","files","i","file","gene_pos"))

