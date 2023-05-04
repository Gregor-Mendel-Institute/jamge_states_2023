library("tidyverse")
library("ComplexHeatmap")
library("valr")
library("ggpubr")
library("circlize")
library("ggalluvial")

## source settings, functions and read in files
source("support_script/setup_fig3.R")

# Important data:
# col = state assignment col 
# ddm = state assignment ddm
# mat = emission matrix ddm and col model 
# gene_pos = genomic position of the genes
# te_gene_pos = genomic positopn of the te_genes
# feature = subset of gff file with the features we want to plot and more human friendly names
# other:
# cs2 = color scheme for states
# m_order = the order of states in emission matrix


######################################################
######################################################
# STEP 1 Generate data frames
######################################################
######################################################

## generate two data frames needed throughout fig3 B-D, sup D, E
# 1, overlap_models: this is the intersect between the col and ddm state assignments, summarized per state combination. 
# That is, overlaps between all existing state combinations:
# state in col, state in ddm, overlap, size of col state, size of ddm size

overlap_models=bed_intersect(col, ddm,suffix = c(".col",".ddm"))%>%
  group_by(state.col,state.ddm)%>%summarise(tot=sum(.overlap),.groups='drop_last')%>%
  mutate(size_col=sum(tot))%>%group_by(state.ddm)%>%mutate(size_ddm=sum(tot))

# the total genome size should be the same 
stopifnot(sum(filter(overlap_models,state.col==state.ddm )$size_col)==sum(overlap_models$tot) 
          & sum(overlap_models$tot)==sum(filter(overlap_models,state.col==state.ddm )$size_ddm))


# 2, from overlap_models, the statistics plotted in fig 3 are extracted.
# Note for this only the cases where the state is the same in col and ddm are of interest.
state_stats=overlap_models%>%ungroup%>%filter(state.col==state.ddm)%>%
  mutate(state=factor(state.col,levels=m_order))%>%
  mutate(logFC=log2(size_ddm/size_col))%>%
  mutate(conserv=tot/size_col)%>%
  mutate(pro_gm_col=100*(size_col/sum(size_col)))%>%
  mutate(index=tot/(size_col+size_ddm-tot))

###############################################
###############################################
## FIGURE 3A emission matrix from ddm + wt concat. model
###############################################
###############################################

## stats. part of genome covered by each state, for col and ddm
size_col=col%>%mutate(length=end-start)%>%group_by(state)%>%
  summarise(cov=sum(length))%>%arrange(state)
size_ddm=ddm%>%mutate(length=end-start)%>%group_by(state)%>%
  summarise(cov=sum(length))%>%arrange(state)
state_size=full_join(size_ddm,size_col,by="state",suffix=c(".ddm",".col"))%>%
  mutate(state=factor(state,levels=m_order))%>%arrange(state) # !! state_size is now order same way as matrix !!

## Heatmap annotations
row_ha = rowAnnotation(bp =anno_barplot(matrix(nc=2,c(state_size$cov.ddm,state_size$cov.col)), # works since state_size is in same order
                                        beside=T,gp = gpar(fill = c("red","darkgreen"))),
                       width=unit(1.5, "cm"))
lgd = Legend(at = 1:2, title = "genotype", legend_gp = gpar(fill = c("darkgreen","red")),labels=c("wt","ddm1"))
row_ha2 = rowAnnotation(state=m_order, col=list(state=cs2),show_legend = c(FALSE),show_annotation_name = FALSE)

# extra check!
stopifnot(sum(rownames(mat[m_order,])!=state_size$state)==0)

fig3a=Heatmap(mat[m_order,],
           cluster_rows=FALSE,
           cluster_columns = FALSE,
           col=cf(10),
           right_annotation = row_ha,
           row_names_side="left",
           left_annotation = row_ha2,
           heatmap_legend_param = list(title="Emission probability",legend_height = unit(10, "cm"), at = seq(0,1,by=0.2),
                                       title_position="leftcenter-rot"))

## plot:
pdf(paste(out_folder,"Fig3A.pdf",sep="/"))
plot(fig3a)
draw(lgd, x = unit(17.3, "cm"), y = unit(0.5, "cm"), just = c("right", "bottom"))
dev.off()

###############################################
###############################################
## FIGURE 3B-D 
###############################################
###############################################

########### Jaccard Index

fig3b=state_stats%>%mutate(state=fct_rev(state))%>%
  ggplot(aes(x=state,y=index,fill=state))+
  theme(legend.position = "none")+
  geom_col()+
  xlab("")+ylab("Jaccard Index")+
  scale_fill_manual(values=cs2)+coord_flip()


########### states with more or less overlap compared to genomewide

# genomewide overlap
gw_overlap=overlap_models%>%filter(state.col==state.ddm)%>%ungroup()%>%
  summarise(overlap=sum(tot),col=sum(size_col),ddm=sum(size_ddm))%>%mutate(same=overlap/col)

fig3c= state_stats%>%mutate(state=fct_rev(state))%>%mutate(conserv=conserv*100)%>%
  ggplot(aes(x=state,y=conserv,fill=state))+
  theme(legend.position = "none",
        axis.text.y=element_blank())+
  geom_col()+
  geom_hline(yintercept =gw_overlap$same*100,col="grey" )+ 
  xlab("")+ylab("Overlap (%)")+
  scale_fill_manual(values=cs2)+coord_flip()

########### differences in size

fig3d=state_stats%>%mutate(state=fct_rev(state))%>%
  ggplot(aes(x=state,y=logFC,fill=state))+
  geom_col()+
  theme(legend.position = "none",
        axis.text.y=element_blank())+
  scale_fill_manual(values=cs2)+
  xlab("")+ylab("logFC(ddm1/wt)")+
  coord_flip()

## plot:
pdf(paste(out_folder,"Fig3BCD.pdf",sep="/"))
plot(ggarrange(fig3b,fig3c,fig3d, 
          ncol = 3))
dev.off()
###############################################
###############################################
## FIGURE 3E-F feature overlap plots
###############################################
###############################################

# overlaps
f_ddm = bed_intersect(ddm,feature)%>%mutate(genotype="ddm")
f_col=bed_intersect(col,feature)%>%mutate(genotype="col")

fig3e = f_col%>%group_by(feature.y,state.x)%>%summarize(ov=sum(.overlap))%>%
  mutate(n=sum(ov),proc=(ov/n)*100)%>%
  ggplot(aes(x=feature.y,y=proc,fill=state.x))+geom_col()+scale_fill_manual(values=cs2)+
  theme(legend.position = "none",
        axis.text.x = element_blank())+ggtitle("WT")+ylab("% covered by each state")+xlab("")

fig3f = f_ddm%>%group_by(feature.y,state.x)%>%summarize(ov=sum(.overlap))%>%
  mutate(n=sum(ov),proc=(ov/n)*100)%>%
  ggplot(aes(x=feature.y,y=proc,fill=state.x))+geom_col()+scale_fill_manual(values=cs2)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle("ddm1")+ylab("% covered by each state")+xlab("")


## plot:
pdf(paste(out_folder,"Fig3EF.pdf",sep="/"))
plot(ggarrange(fig3e,fig3f, heights=c(1.1,1.45),
          nrow = 2))
dev.off()
#########################################################
#########################################################
## FIGURE 3 sup D Jaccard Index vs log2FC
#########################################################
#########################################################

figS3d=state_stats%>%ggplot(aes(x=index,y=logFC,fill=state.col))+geom_point(size=3,shape=21)+
  scale_fill_manual(values=cs2)+theme(legend.position="top")+
  #legend.text=element_text(family="Times New Roman",face = "bold"))+
  xlab("Jaccard index")+ylab("logFC(ddm1/wt)")+
  guides(fill=guide_legend(title="",nrow=2))

## plot:

pdf(paste(out_folder,"FigS3D.pdf",sep="/"))
plot(figS3d)
dev.off()

#########################################################
#########################################################
## FIGURE 3 sup E Jaccard Index vs log2FC per PC and TE genes
#########################################################
#########################################################
## This is a bit ad hoc analysis as the TE/PC genes positions need to be included in the analysis 
# This is with all TE genes... in the first version of the paper only the selected TE genes (used for quantiles) were used

bpc=c("darkgreen","red")
both_gene_pos=rbind(te_gene_pos%>%mutate(type="TE"),gene_pos%>%mutate(type="PC"))

## generate empty bins to fill (this are the ones chrom hmm used)
ch_info=col%>%group_by(chrom)%>%summarise(max=max(end))
hmmbins = tibble(chrom=rep(ch_info$chrom,ch_info$max/200),
                 start=unlist(sapply(ch_info$max,function(x) seq(0,x-200,200))),
                 end=unlist(sapply(ch_info$max,function(x) seq(200,x,200))))

## fill bins with genotype specific states from above
hbin_col=bed_intersect(hmmbins,col)%>%
  filter(.overlap>0)%>%
  dplyr::rename(start=start.x,end=end.x,state=state.y)%>%
  select(chrom,start,end,state)
hbin_ddm=bed_intersect(hmmbins,ddm)%>%
  filter(.overlap>0)%>%
  dplyr::rename(start=start.x,end=end.x,state=state.y)%>%
  select(chrom,start,end,state)


both_state=bed_intersect(full_join(hbin_ddm,hbin_col,suffix=c(".ddm",".col"),by=c("chrom","start","end")),both_gene_pos,suffix=c(".state",".gene"))

gene_state_stats=both_state%>%group_by(type.gene,state.col.state,state.ddm.state)%>%summarise(tot=sum(.overlap))%>%
  mutate(size_col=sum(tot))%>%ungroup()%>%group_by(type.gene,state.ddm.state)%>%mutate(size_ddm=sum(tot))%>%filter(state.col.state==state.ddm.state)%>%
  mutate(logFC=log2(size_ddm/size_col))%>% mutate(index=tot/(size_col+size_ddm-tot))
  

figS3e1 =gene_state_stats%>%filter(type.gene=="PC")%>%ggplot(aes(x=index,y=logFC,fill=state.col.state))+geom_point(size=3,shape=21)+
  scale_fill_manual(values=cs2)+theme(legend.position="bottom")+
  geom_hline(yintercept = 0,color="grey")+
  xlab("Jaccard index")+ylab("logFC(ddm1/wt)")+
  guides(fill=guide_legend(title="",nrow=2))+
  xlim(c(0,1))+
  ylim(c(-6,6))+
  ggtitle("Protein coding genes")

figS3e2 = gene_state_stats%>%filter(type.gene=="PC")%>%ungroup()%>%select(state.col.state,size_col,size_ddm)%>%
  pivot_longer(-state.col.state)%>%group_by(name)%>%mutate(proc=value/sum(value)*100)%>%
  ggplot(aes(x=state.col.state,y=proc,fill=name))+geom_col(position="dodge")+
  scale_fill_manual(values=bpc,labels=c("WT","ddm1"))+ggtitle("Protein coding genes")+
  theme(legend.position = "bottom",legend.title = element_blank())+ylab("% overlap")+xlab("")
  

figS3e3 =gene_state_stats%>%filter(type.gene=="TE")%>%ggplot(aes(x=index,y=logFC,fill=state.col.state))+geom_point(size=3,shape=21)+
  scale_fill_manual(values=cs2)+theme(legend.position="bottom")+
  geom_hline(yintercept = 0,color="grey")+
  xlab("Jaccard index")+ylab("logFC(ddm1/wt)")+
  guides(fill=guide_legend(title="",nrow=2))+
  xlim(c(0,1))+
  ylim(c(-6,6))+
  ggtitle("TE genes")+theme(legend.position = "none")

figS3e4 = gene_state_stats%>%filter(type.gene=="TE")%>%ungroup()%>%select(state.col.state,size_col,size_ddm)%>%
  pivot_longer(-state.col.state)%>%group_by(name)%>%mutate(proc=value/sum(value)*100)%>%
  ggplot(aes(x=state.col.state,y=proc,fill=name))+geom_col(position="dodge")+
  scale_fill_manual(values=bpc,labels=c("WT","ddm1"))+ggtitle("TE genes")+
  theme(legend.position = "none",legend.title = element_blank())+ylab("% overlap")+xlab("")

## plot:

pdf(paste(out_folder,"FigS3E.pdf",sep="/"))
plot(ggarrange(figS3e4,figS3e3,figS3e2,figS3e1,heights = c(1,1.3),
          ncol = 2,nrow=2))
dev.off()



#######################################################
#########################################################
## FIGURE 3 sup F Alluvium 
#########################################################
#########################################################
# code inspired by:
## https://stackoverflow.com/questions/68487536/how-to-align-and-label-the-stratum-in-ggalluvial-using-ggrepel-or-otherwise

data_long = select(overlap_models,state.col,state.ddm,tot )
## only changing?
#data_long = filter(data_long,state.col!=state.ddm)
names(data_long)=c("state.col","state.ddm","tot")
levs1 <- levels(data_long$state.col) 
levs2 <- levels(data_long$state.ddm)
res1 <- unique(data_long$state.col)
res2 <- unique(data_long$state.ddm)
cond1_cols <- cs2[levs1[levs1 %in% res1]]
cond2_cols <- cs2[levs2[levs2 %in% res2]]
columnCols <- c(cond1_cols, cond2_cols)
stratCols <- c(rev(cond1_cols), rev(cond2_cols))
figS3f=ggplot(data = data_long,
       aes(axis1 =  state.col, axis2 =  state.ddm, y = tot)) +
  geom_alluvium(width = 0,aes(fill=state.col)) +
  geom_stratum(fill = paste0(stratCols),col="white",width = 1/3,)+
  scale_x_discrete(limits = c("state.col", "state.ddm"),expand = c(1,1),labels=c("WT", "ddm"))+
  scale_fill_manual(values=cs2)+
  scale_y_continuous(breaks = NULL,name = "")+
  theme(legend.position = "top")+
  guides(fill=guide_legend(title="",nrow=2))

pdf(paste(out_folder,"FigS3F.pdf",sep="/"))
plot(figS3f)
dev.off()
#
