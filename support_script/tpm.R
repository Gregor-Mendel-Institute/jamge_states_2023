source("support_script/functions.R")

## script that loads the expression data for wt and ddm1 leaves. And defines the bins.

### data and annotations:
load("data/expression/tpm_leaves.Rdata")
gff=get_our_gff() 
te_gene_pos=get_te_gene_pos(gff)
# extra test
stopifnot(length(setdiff(te_gene_pos$chrom,c("1","2","3","4","5")))==0)

### TE genes only
te_tpm = left_join(te_gene_pos,tpm_leaves)
stopifnot(all_equal(te_tpm%>%select(-col_tpm,-ddm_tpm),te_gene_pos))

### expressed in ddm1 not in col
bins=te_tpm%>%filter(col_tpm==0 & ddm_tpm>0.1)%>%mutate(bin = cut_number(ddm_tpm, n = 4))%>%
  select(geneId,bin)
### not expressed (zbin = zero bin)
zbin=te_tpm%>%filter(col_tpm==0 & ddm_tpm<=0.1)%>%mutate(bin = "[0]")%>%
  select(geneId,bin)
### combine bins
bins=rbind(bins,zbin)%>%mutate(bin=fct_relevel(bin,"[0]",after=0))
# new names
levels(bins$bin)=c("no exp.","1st","2nd","3rd","4th")
table(bins$bin)

tpm_te_table = full_join(te_tpm,bins)
# test, that fct ordered correctly. NA means other, remove. 
test = tpm_te_table%>%filter(!is.na(bin))%>%group_by(bin)%>%summarise(min=min(ddm_tpm),max=max(ddm_tpm),n=n())

# testing that the order is right : no exp < 1st < 2nd etc
stopifnot( all_equal(test,test%>%arrange(min),ignore_row_order = F))
stopifnot( all_equal(test,test%>%arrange(max),ignore_row_order = F))

# save as bed file for e.g. deeptool plots)
tobed =tpm_te_table%>%select(chrom,start,end ,geneId ,bin,strand)
write_delim(tobed,file="tegenes_tpm_binned.bed",delim = "\t")

# clean up
rm(list=c("bins","zbin","te_tpm","tpm_leaves","test","gff","te_gene_pos","tobed"))


