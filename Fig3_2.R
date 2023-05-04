library("tidyverse")
library("ComplexHeatmap")
library("valr")
library("ggpubr")

## source settings, functions and read in files
source("support_script/setup_fig3.R")
# plus extra files from leaf-seedling model and large seedling model
source("support_script/setup_fig3S.R")

# Important data:
# col = state assignment col 
# ddm = state assignment ddm
# large_col = state assignment extensive model wt
# lcol = state assignment col leaves
# scol = state assignment col seedling
# mat = emission matrix ddm and col model 
# gene_pos = genomic position of the genes
# te_gene_pos = genomic positopn of the te_genes
# feature = subset of gff file with the features we want to plot and more human friendly names
# other:
# cs2 = color scheme for states
# m_order = the order of states in emission matrix



# some extra colors only needed in this script
gcl=c("darkgreen","steelblue","red","grey")
names(gcl)=c("E","F","H","I")

#########################################################
## FIGURE 3 sup A Emission matrix leaf+seedling
#########################################################
lsize=lcol%>%mutate(length=end-start)%>%mutate(state=as.numeric(gsub("[A-Z]","",state)))%>%group_by(state)%>%
  summarise(cov=sum(length))%>%arrange(state)
ssize=scol%>%mutate(length=end-start)%>%mutate(state=as.numeric(gsub("[A-Z]","",state)))%>%group_by(state)%>%
  summarise(cov=sum(length))%>%arrange(state)

row_ha = rowAnnotation(bp =anno_barplot(matrix(nc=2,c(lsize$cov[rord],ssize$cov[rord])),
                                        beside=T,gp = gpar(fill = c("darkolivegreen4","darkolivegreen1"))),width=unit(1.5, "cm"))
row_ha2 = rowAnnotation("chromatin type"=c(rep("H",5),rep("F",5),rep("E","4"),"I"),col=list("chromatin type"=c("E"="darkgreen","F"="steelblue","H"="red","I"="grey")))

lgd = Legend(at = 1:2, title = "tissue", legend_gp = gpar(fill = c("darkolivegreen4","darkolivegreen1")),labels=c("leaf","seedling"))
figS3a=Heatmap(mat, 
          cluster_rows=FALSE,
          cluster_columns = FALSE,
          col=cf(10),
          heatmap_legend_param = list(title="Emission\nprobability",legend_height = unit(10, "cm"), at = seq(0,1,by=0.2),
                                      title_position="leftcenter-rot"),right_annotation = row_ha,
          left_annotation = row_ha2)

#plot:
pdf(paste(out_folder,"FigS3A.pdf",sep="/"))
plot(figS3a)
draw(lgd, x = unit(13, "cm"),y = unit(0.5, "cm"), just = c("right", "bottom"))
dev.off()

#########################################################
## FIGURE 3 sup B Map to extensive model chromatin types from leaf-seedling model
#########################################################
# the order we want the states
scol = mutate(scol,state=factor(gsub("[A-Z]","",state),levels = rord))%>%mutate(state=as.numeric(state))
model_comp=
  bed_intersect(scol, large_col,suffix = c(".sub",".large"))%>%
  mutate(state_group=gsub("[0-9]","",state.large))%>%
  group_by(state.sub,state_group)%>%
  summarise(tot=sum(.overlap))%>%
  mutate(state_length=sum(tot))%>%
  group_by(state_group)%>%
  mutate(type_length=sum(tot))%>%
  mutate(proc = (tot/state_length)*100)

stopifnot(sum(model_comp$tot)==sum(unique(model_comp$state_length)) & 
            sum(model_comp$tot)==sum(unique(model_comp$type_length)))

figS3b=model_comp%>%select(state.sub,state_group,proc)%>%distinct()%>%
  ggplot(aes(x=state_group,y=proc,fill=state_group))+
  geom_col()+facet_wrap(~state.sub)+
  scale_fill_manual(values=gcl,name = "", labels = c("Euchromatin", "Facultative", "Constitutive","Intergenic"))+
  theme(legend.position = "top",axis.text.x = element_blank())+ylab("Overlap (%)")+xlab("")

#plot:
pdf(paste(out_folder,"FigS3B.pdf",sep="/"))
plot(figS3b)
dev.off()  
#########################################################
## FIGURE 3 sup C Map to extensive model chromatin types from ddm1-col model
#########################################################

model_comp=
  bed_intersect(col, large_col,suffix = c(".sub",".large"))%>%
  mutate(state_group=gsub("[0-9]","",state.large))%>%
  group_by(state.sub,state_group)%>%
  summarise(tot=sum(.overlap))%>%
  mutate(state_length=sum(tot))%>%
  group_by(state_group)%>%
  mutate(type_length=sum(tot))%>%
  mutate(proc = (tot/state_length)*100)

stopifnot(sum(model_comp$tot)==sum(unique(model_comp$state_length)) & 
            sum(model_comp$tot)==sum(unique(model_comp$type_length)))

figS3c=model_comp%>%select(state.sub,state_group,proc)%>%distinct()%>%
  mutate(state=fct_relevel(state.sub,"fhI",after=13))%>%
  ggplot(aes(x=state_group,y=proc,fill=state_group))+
  geom_col()+facet_wrap(~state)+
  scale_fill_manual(values=gcl,name = "", labels = c("Euchromatin", "Facultative", "Constitutive","Intergenic"))+
  theme(legend.position = "top",axis.text.x = element_blank())+ylab("Overlap (%)")+xlab("")

#plot:
pdf(paste(out_folder,"FigS3C.pdf",sep="/"))
plot(figS3c)
dev.off()


