library("RColorBrewer")

## each pdf to a4 format, to be outputed to folder
pdf.options(paper = "a4")
out_folder="pdf_figures"
if (!dir.exists(out_folder)){
  dir.create(out_folder)
}

## minimal theme as default
theme_set(theme_minimal())

## emission heatmaps in grey scale 
cf = colorRampPalette(c("grey100","grey10")) 

####################################################
## Extensive 26 state model
####################################################
# color scheme as in paper first version
cs=c(brewer.pal(7,"Reds")[7:2],brewer.pal(7,"BuPu")[2:7],brewer.pal(4,"Greys")[2:4],brewer.pal(3,"YlOrBr")[2:1],brewer.pal(9,"Greens"))

## the order that is used in the paper
map_to_paper = data.frame(
  from=c("23","25","24","22","21","26",
         "1","3","2","4","5","6",
         "18","19","20",
         "7","8","9","10","11","12","13","14","17","15","16"),
  to=factor(c("H1","H2","H3","H4","H5","H6",
              "F1","F2","F3","F4","F5","F6",
              "I1","I2","I3",
              "E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11"),
            levels = c("E11","E10","E9","E8","E7","E6","E5","E4","E3","E2","E1",
                       "I3","I2","I1",
                       "F6","F5","F4","F3","F2","F1",
                       "H6","H5","H4","H3","H2","H1")),
  cols=cs)

names(cs)=map_to_paper$to # to get the right colors in the plots

sam_ord=c("HTB2","HTB1","HTB4","HTB3","H1",
          "H3K27me1","H3K9me2","H2AW6","H3K9me1",
          "H2AW7", "H4K20me1","H3","HTR13",
          "H3K4me1","HTR5","H3K27me3","H2AK121ub",
          "H2AZ11","H2AZ9","H3K4me3","H3K14Ac",
          "H3K9Ac","H2A13","H2A2","H2AX","H2Bub","H3K36me3")

####################################################
## ddm vs col concatenated 
####################################################
# colors scheme for ddm/col data
cs2 =c( "#be4535", "#d6776b", "#e7b1aa","#db24d7","#ade2fd","#69cbfc","#05a8fa","#5c3dfd",
        "#7b09fd","#5802bc","#969696","#525252","#CCCCCC","#E8DA8D","#D3E3CF","#87B092")
names(cs2)=c(
  "hI","hII","hIII","fhI","fI","fII","fIII","fIV",
  "fV","fVI","ihI","ihII","iI","eI","eII","eIII"
)
## map chromHMM names to ours
map_mod2=data.frame(state=paste0("E",1:16),
                    ns=factor(c("fII","fIII","fI","fhI",
                      "ihI","hIII","hI","hII",
                      "ihII","fIV","fVI","fV",
                      "iI","eI","eII","eIII"),
                      levels=c(
                          "eI","eII","eIII","fI",
                          "fII","fIII","fIV","fV",
                          "fVI","fhI","hI","hII","hIII",
                          "iI","ihI","ihII")))

# matrix order for plotting
m_order =c("hI","hII","hIII","ihI","ihII","fI","fII","fIII","fIV","fV","fVI","fhI","iI","eI","eII","eIII")
# order of marks for heatmap
sampleo=c("H3K9me1","H3K9me2","H2AW6","H2AW7","H2AZ11","H2AZ",
          "H2AX","H3K27me3","H3K36me3","H3", "H1")

####################################################
## leaf vs seedling concatenated 
####################################################
# matrix row order for plotting
rord=c(1:8,10,9,12:15,11) # order of states 

