library("tidyverse")
library("valr")
library("ComplexHeatmap")
library("ggpubr")

## source settings, functions and read in files
source("support_script/setup_fig2.R")
# Important data:
# large_col = state assignment extensive model
# cg/chg/chh = processed methylation data
# dnase = processed dnase data
# mnase = processed mnase data
# gff = gff file from araport11 downloaded June 2016
# gff_w_introns = gff file with introns include
# jindex = those are the jaccard indices generated from the confusion matrices 
# tpm_seedlings = proceessed experssion data

# other:
# cs = color scheme for states

###############################################
###############################################
## FIGURE 2B
###############################################
###############################################

## extract features of interest
feature = gff_to_feature(gff_w_introns)

## overlap with states
f_col=bed_intersect(large_col,feature)

## create plot
fig2b = f_col%>%group_by(feature.y,state.x)%>%
    summarize(ov=sum(.overlap))%>%
    mutate(n=sum(ov),proc=(ov/n)*100)%>%
    ggplot(aes(x=feature.y,y=proc,fill=state.x))+
    geom_col()+scale_fill_manual(values=cs)+
    theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ylab("% covered by each state")+xlab("")

pdf(paste(out_folder,"Fig2B.pdf",sep="/"))
plot(fig2b)
dev.off()

###############################################
###############################################
## FIGURE 2C-E
###############################################
###############################################

state_tpm=bed_intersect(large_col,tpm_seedlings)%>%rename(to=to.x,tpm=tpm.y)
rrange=state_tpm%>%group_by(to)%>%mutate(qq=find_outlier_range(tpm))%>%ungroup()%>%summarise(max=max(qq))

fig2c=state_tpm%>%
  ggplot(aes(x=to,y=tpm,fill=to))+
  geom_boxplot(outlier.shape = NA,)+
  scale_fill_manual(values=cs)+
  scale_y_continuous(limits=c(0,0.8*rrange$max),breaks = seq(0,rrange$max,10))+
  coord_flip()+
  ylab("Gene expression (TPM)")+
  xlab("")+
  theme(legend.position = "none", 
        axis.text.y=element_blank(),
        panel.border=element_rect(colour = "black", fill=NA, size=1),
        axis.ticks.y=element_blank())


fig2d=cg%>%mutate(mean=mean*100)%>% #to %
  ggplot(aes(x=to,y=mean,fill=to))+
  geom_boxplot(outlier.shape = NA,)+
  scale_fill_manual(values=cs)+
  coord_flip()+ylab("% CG methylation")+xlab("")+
  theme(legend.position = "none", 
        panel.border=element_rect(colour = "black", fill=NA, size=1),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

rrange=dnase%>%group_by(to)%>%mutate(qq=find_outlier_range(value))%>%ungroup()%>%summarise(max=max(qq))

fig2e=dnase%>% 
  ggplot(aes(x=to,y=value,fill=to))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=cs)+
  ylim(c(0,rrange$max))+
  coord_flip()+
  ylab("DNase 1 signal")+xlab("")+
  theme(legend.position = "none", 
        axis.text.y=element_blank(),
        panel.border=element_rect(colour = "black", fill=NA, size=1),
        axis.ticks.y=element_blank())


pdf(paste(out_folder,"Fig2CDE.pdf",sep="/"))
plot(ggarrange(fig2c,fig2d,fig2e,nrow =1))
dev.off()

###############################################
###############################################
## FIGURE 2F
###############################################
###############################################

row_ha = rowAnnotation(state=names(cs), col=list(state=cs),show_legend = c(FALSE),show_annotation_name = FALSE)

fig2f = Heatmap(jindex,cluster_rows = FALSE,cluster_columns =  FALSE,row_names_gp = gpar(fontsize = 12),
        row_names_side="left",
        left_annotation = row_ha,
        col=  brewer.pal(9, "YlOrBr"),
        row_split = factor(str_extract(rownames(jindex),"^[A-Z]"),levels= unique(str_extract(rownames(jindex),"^[A-Z]"))),
        heatmap_legend_param = list(title="Jaccard\nIndex",legend_height = unit(3, "cm"), at = seq(0,1,by=0.2)),
        row_gap = unit(1.5, "mm"),border = F)

pdf(paste(out_folder,"Fig2F.pdf",sep="/"))  
plot(fig2f)
dev.off()

###############################################
###############################################
## FIGURE Sup 2C (A and B - Bhagyshree)
###############################################
###############################################

figS2c=chg%>%mutate(to=fct_rev(to),mean=mean*100)%>% # to %
  ggplot(aes(x=to,y=mean,fill=to))+
  geom_boxplot(outlier.shape = NA,varwidth = T)+
  scale_fill_manual(values=cs)+
  ylab("% CHG methylation")+xlab("")+
  scale_y_continuous(breaks=seq(0,100,20))+
  theme(legend.position = "none",
        panel.border=element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(paste(out_folder,"FigS2C.pdf",sep="/"))
plot(figS2c)
dev.off()

###############################################
###############################################
## FIGURE Sup 2D
###############################################
###############################################

rrange=chh%>%group_by(to)%>%mutate(qq=find_outlier_range(mean))%>%ungroup()%>%summarise(max=max(qq))

figS2d=chh%>%mutate(to=fct_rev(to),mean=mean*100)%>% # to in %
  ggplot(aes(x=to,y=mean,fill=to))+
  geom_boxplot(outlier.shape = NA,varwidth = T)+
  scale_fill_manual(values=cs)+
  scale_y_continuous(limits=c(0,rrange$max*100),breaks=seq(0,100,20))+
  ylab("% CHH methylation")+xlab("")+
  theme(legend.position = "none",
        panel.border=element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


pdf(paste(out_folder,"FigS2D.pdf",sep="/"))
plot(figS2d)
dev.off()

###############################################
###############################################
## FIGURE Sup 2E
###############################################
###############################################

rrange=mnase%>%group_by(to)%>%mutate(qq=find_outlier_range(X7))%>%ungroup()%>%summarise(max=max(qq))

figS2e=mnase%>%mutate(to=fct_rev(to),X7=X7)%>%
  ggplot(aes(x=to,y=X7,fill=to))+
  geom_boxplot(outlier.shape = NA,varwidth = T)+
  scale_fill_manual(values=cs)+
  ylab("MNase-seq Read Density")+xlab("")+
  scale_y_continuous(limits=c(0,0.8*rrange$max),breaks = seq(0,rrange$max,200))+
  theme(legend.position = "none",
        panel.border=element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(paste(out_folder,"FigS2E.pdf",sep="/"))
plot(figS2e)
dev.off()

###############################################
###############################################
## FIGURE Sup 2F
###############################################
###############################################

genome_dist = large_col%>%mutate(length=end-start)%>%group_by(state)%>%summarise(tot=sum(length))%>%
  ungroup()%>%mutate(all=sum(tot),proc=tot/all*100)

pdf(paste(out_folder,"FigS2F.pdf",sep="/"))
pie(genome_dist$proc,labels = paste(genome_dist$state,paste0(round(genome_dist$proc,1),"%"),sep=":"), col = cs)
dev.off()

###############################################
###############################################
## FIGURE Sup 2G
###############################################
###############################################

rna_feature=gff_to_rnatype(gff)
r_col=bed_intersect(large_col,rna_feature)

figS2g = r_col%>%group_by(feature.y,state.x)%>%
  summarize(ov=sum(.overlap))%>%
  mutate(n=sum(ov),proc=(ov/n)*100)%>%
  ggplot(aes(x=feature.y,y=proc,fill=state.x))+
  geom_col()+scale_fill_manual(values=cs)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("% covered by each state")+xlab("")


pdf(paste(out_folder,"FigS2G.pdf",sep="/"))
plot(figS2g)
dev.off()

###############################################
###############################################
## FIGURE Sup 2H
###############################################
###############################################

rrange=large_col%>%mutate(length=end-start)%>%group_by(state)%>%
  mutate(qq=find_outlier_range(length))%>%ungroup()%>%summarise(max=max(qq),min=min(length))

figS2h=large_col%>%mutate(length=end-start)%>%
  ggplot(aes(x=state,y=length,fill=state))+geom_boxplot(varwidth = T,outlier.shape = NA)+
  scale_fill_manual(values=cs)+
  scale_y_continuous(limits=c(rrange$min,0.7*rrange$max),breaks=seq(0,rrange$max,500))+
  ylab("State length (bp)")+xlab("")+
  theme(legend.position = "none",
        panel.border=element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(paste(out_folder,"FigS2H.pdf",sep="/"))
plot(figS2h)
dev.off()

