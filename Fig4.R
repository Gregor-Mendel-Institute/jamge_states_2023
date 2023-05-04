library("tidyverse")
library("ComplexHeatmap")
library("valr")
library("ggpubr")

## source settings, functions and read in files
source("support_script/setup_fig4.R")
source("support_script/tpm.R")

# Important data:
# col = state assignment col 
# ddm = state assignment ddm
# tpm_te_table = table with each TE genes,positions, id, col/ddm tpm and bin type (no.expression, 1:4 quantile,NA=expressed in col)
# TE genes expressed in Col (tpm>0) are removed before plotting!

######################################################
######################################################
# STEP 1 Generate data frames
######################################################
######################################################

te_col = bed_intersect(col,tpm_te_table,suffix=c(".state",".te"))%>%
  filter(!is.na(bin.te)) # use only TE genes not expressed in col
  
te_col_states=te_col%>%group_by(bin.te,state.state)%>%
  summarise(tot=sum(.overlap))%>%
  mutate(tot_l=sum(tot))%>%
  mutate(proc=(tot/tot_l)*100)
   

te_ddm = bed_intersect(ddm,tpm_te_table,suffix=c(".state",".te"))%>%
  filter(!is.na(bin.te)) # use only TE genes not expressed in col
te_ddm_states=te_ddm%>%group_by(bin.te,state.state)%>%
  summarise(tot=sum(.overlap))%>%
  mutate(tot_l=sum(tot))%>%
  mutate(proc=(tot/tot_l)*100)

## 
te_comp_stats=full_join(te_col%>%group_by(bin.te,geneId.te)%>%summarise(types=length(unique(state.state)),n=n()),
                        te_ddm%>%group_by(bin.te,geneId.te)%>%
                          summarise(types=length(unique(state.state)),n=n()),by=c("bin.te","geneId.te"),suffix=c(".col",".ddm"))%>%
              pivot_longer(-c(bin.te,geneId.te))



######################################################
######################################################
# Figure 4B States covering TE genes per bin DDM and Col
######################################################
######################################################

fig4b1=ggplot(te_col_states,aes(x=bin.te,y=proc,fill=state.state))+geom_col()+scale_fill_manual(values=cs2)+
  theme(legend.position = "none",legend.title = element_blank(),axis.text.x=element_blank())+
  ggtitle("WT")+
  ylab("% covered by each state")+
  xlab("")+
  guides(fill =guide_legend(nrow=2))+
  scale_y_continuous(breaks = c(0,25,50,75,100),limits = c(0,110))+
  annotate("text",
           x = c(1:5),
           y = 95,
           label =paste("n",sep="=",dplyr::count(tpm_te_table%>%filter(!is.na(bin)),bin)$n),
           col = "black",
           vjust = - 1)

fig4b2=ggplot(te_ddm_states,aes(x=bin.te,y=proc,fill=state.state))+geom_col()+scale_fill_manual(values=cs2)+
  theme(legend.position = "bottom",legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("ddm1")+ylab("% covered by each state")+xlab("")+guides(fill =guide_legend(nrow=2))

#plot:
pdf(paste(out_folder,"Fig4B.pdf",sep="/"))
plot(ggarrange(fig4b1,fig4b2, heights=c(1,1.35),
          nrow = 2))
dev.off()

######################################################
######################################################
# Figure 4C
######################################################
######################################################

find_outlier_range <- function(x) {
  return( quantile(x, .75,na.rm=T) + 1.5*IQR(x,na.rm=T))
}
rrange=te_ddm%>%group_by(state.state)%>%mutate(qq=find_outlier_range(ddm_tpm.te))%>%ungroup()%>%summarise(max=max(qq))


fig4c=te_ddm%>%ggplot(aes(x=state.state,y=ddm_tpm.te,fill=state.state))+
  geom_boxplot(outlier.shape = NA)+scale_fill_manual(values=cs2)+
  scale_y_continuous(limits=c(0,0.8*rrange$max),breaks=seq(0,rrange$max,5))+  ## made the y range slightly shorter to avoid empty space on top
  ylab("Expression of TE genes in TPM")+xlab("")+
  theme(legend.position = "none",
        panel.border =element_rect(colour = "black", fill=NA, size=1))+ggtitle("ddm1")


pdf(paste(out_folder,"Fig4C.pdf",sep="/"))
plot(fig4c)
dev.off()
######################################################
######################################################
# Figure 4 sup C 
######################################################
######################################################
rrange=te_col%>%mutate(length=end.te-start.te)%>%
  select(geneId.te,bin.te,length)%>%distinct()%>%
  group_by(bin.te)%>%mutate(qq=find_outlier_range(length))%>%ungroup()%>%summarise(max=max(qq))


figS4c=te_col%>%mutate(length=end.te-start.te)%>%
  select(geneId.te,bin.te,length)%>%distinct()%>%
  ggplot(aes(x=bin.te,y=length,fill=bin.te))+geom_boxplot(outlier.shape = NA)+
    scale_fill_manual(values=qc)+
    scale_y_continuous(limits=c(0,rrange$max))+
    ylab("TE gene length (bp)")+xlab("")+theme(legend.position = "none")


pdf(paste(out_folder,"FigS4C.pdf",sep="/"))
plot(figS4c)
dev.off()
######################################################
######################################################
# Figure 4 sup D
######################################################
######################################################

rrange=te_comp_stats%>%filter(name%in%c("types.col","types.ddm"))%>%
  group_by(bin.te)%>%mutate(qq=find_outlier_range(value))%>%ungroup()%>%summarise(max=max(qq))


figS4d1=te_comp_stats%>%
  filter(name%in%c("types.col","types.ddm"))%>%
  mutate(genotype=ifelse(name=="types.ddm","ddm1","WT"))%>%
  ggplot(aes(x=bin.te,y=value,fill=bin.te))+
    geom_boxplot(outlier.shape = NA)+
    facet_wrap(~genotype)+
    scale_fill_manual(values=qc)+
    scale_y_continuous(limits=c(0,rrange$max))+
    theme(legend.position = "none",
          legend.title = element_blank(),
          axis.text.x=element_blank(),
          strip.text=element_text(face = "bold",size = 12),
          axis.title.y = element_text(margin = margin(t = 0, r = 22, b = 0, l = 0)))+
  ylab("# state types")+
  xlab("")

rrange=te_comp_stats%>%filter(name%in%c("n.col","n.ddm"))%>%
  group_by(bin.te)%>%mutate(qq=find_outlier_range(value))%>%ungroup()%>%summarise(max=max(qq))


figS4d2=te_comp_stats%>%filter(name%in%c("n.col","n.ddm"))%>%
  mutate(genotype=ifelse(name=="n.ddm","ddm1","WT"))%>%
  ggplot(aes(x=bin.te,y=value,fill=bin.te))+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~genotype)+
  scale_fill_manual(values=qc)+
  scale_y_continuous(limits=c(0,rrange$max),breaks=seq(0,rrange$max,2))+
  theme(legend.position = "none",legend.title = element_blank(),
        axis.text.x=element_blank(),strip.text = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 17, b = 0, l = 0)))+
  ylab("# state changes")+xlab("")

rrange=rbind(te_col%>%select(bin.te,.overlap)%>%mutate(genotype="WT"),
             te_ddm%>%select(bin.te,.overlap)%>%mutate(genotype="ddm"))%>%
  group_by(bin.te)%>%mutate(qq=find_outlier_range(.overlap))%>%ungroup()%>%summarise(max=max(qq),min=min(.overlap))

figS4d3 = rbind(te_col%>%select(bin.te,.overlap)%>%mutate(genotype="WT"),
             te_ddm%>%select(bin.te,.overlap)%>%mutate(genotype="ddm"))%>%
  ggplot(aes(x=bin.te,y=.overlap,fill=bin.te))+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~genotype)+
  scale_fill_manual(values=qc)+
  scale_y_continuous(limits=c(0,rrange$max),breaks=seq(0,rrange$max,500))+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.x=element_blank(),
        strip.text = element_blank())+
  ylab("State length (bp)")+xlab("")
  

#plot:
pdf(paste(out_folder,"FigS4D.pdf",sep="/"))
plot(ggarrange(figS4d1,figS4d2,figS4d3, heights = c(0.9,1,1.3),
          nrow = 3))
dev.off()
######################################################
######################################################
# Figure 4 sup E
######################################################
######################################################
rrange=te_ddm%>%group_by(geneId.te)%>%
  mutate(part=ifelse(strand.te=="-",n():1,1:n()))%>%
  filter(part==1)%>%group_by(state.state)%>%mutate(qq=find_outlier_range(ddm_tpm.te))%>%ungroup()%>%summarise(max=max(qq))

figS4e=te_ddm%>%group_by(geneId.te)%>%
  mutate(part=ifelse(strand.te=="-",n():1,1:n()))%>%
  filter(part==1)%>%
  ggplot(aes(x=state.state,y=ddm_tpm.te,fill=state.state))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=cs2)+
  scale_y_continuous(limits=c(0,0.8*rrange$max),breaks=seq(0,rrange$max,5))+  ## made the y range slightly shorter to avoid empty space on top
  ylab("ddm1 expression (tpm)")+
  xlab("")+
  theme(legend.position ="none",
        panel.border = element_rect(colour = "black", fill=NA, size=1),)+
  ggtitle("State at TSS in ddm1")

pdf(paste(out_folder,"FigS4E.pdf",sep="/"))
plot(figS4e)
dev.off()
