library(tidyverse)
library(valr)

#####################################
## functions used in several figures
#####################################

## function takes emission matrix as read in (raw) and rearanges the columns according to sample order (samo)
## outputs a matrix that can be plotted 
raw_to_matrix=function(raw,samo){
  mat = raw%>%select(-`state (Emission order)`)
  mat = as.matrix(mat)
  rownames(mat)=raw$`state (Emission order)`
  mat=mat[,samo]
}

## reads in ChromHMM output bed files, renames column X4 to state and excludes mitochondria and chloroplast
## outputs tibble with genomic coordinates and state assignment
read_and_fix_states=function(file,crm=c("Pt","Mt")){
  states = read_bed(file)%>%dplyr::rename(state=X4)
  states = filter(states,!chrom%in%crm)
  states
  }

## reads in the gff file used in this paper (included in the data subfolder), set useful colnames, excludes mitochondria and chloroplast
## extracts gene id
get_our_gff=function(crm=c("ChrC","ChrM")){
  gff = read_delim("data/gff/Araport11_GFF3_genes_transposons.201606.gff",delim="\t",comment = "#",col_names = FALSE)
  colnames(gff)=c("chrom","Source","feature","start","end","X6","strand","X8","info")
  gff = filter(gff,!chrom%in%c(crm))
  # files from chromHMM uses 1:5, Mt, Pt rather than Chr1:5, ChrC, ChrM,
  gff = gff%>%mutate(chrom=str_match(chrom,"^Chr(.+)")[,2])
  gff= gff%>%mutate(geneId=str_match(info,"^ID=([A-Za-z0-9]+)")[,2])
  gff
}

## adds introns to gff. introns= regions overlapping mRNA but not exons, assigns a geneId to each intro
## output tibble with extra introns 
add_intron=function(gff){
  intron_strand=list()
  g_size = gff%>%group_by(chrom)%>%summarize(size=max(end)) # chrom size
  for (strand in c("+","-")){
    intergenic=bed_complement(gff%>%filter(feature=="mRNA" & strand==!!strand),g_size)
    exon = gff%>%filter(feature=="exon" & strand==!!strand)%>%dplyr::select(chrom,start,end) # exons on + strand
    both = rbind(intergenic,exon) 
    intron_strand[[strand]]= bed_complement(both,g_size)%>%mutate(feature="intron",strand=strand)
  }
  introns = do.call("rbind",intron_strand)
  genes = filter(gff,feature=="mRNA")%>%select(chrom,start,end,feature,strand,geneId)%>%distinct()
  an_introns = bed_intersect(introns,genes)%>%filter(.overlap>0)%>%select(-start.y,-end.y)%>%distinct()%>%
    filter(strand.x==strand.y)%>%rename(start=start.x,end=end.x,feature=feature.x,geneId=geneId.y,strand=strand.x)%>%
    select(-feature.y,-strand.y,-.overlap)%>%mutate(Source="-",X6=".",X8=".",info="")
  gfft = full_join(gff,an_introns)
  gfft
}

## extract protein coding genes from gff
get_gene_pos=function(gff){
  gene_pos=gff%>%filter(feature=="gene",str_detect(info,"protein_coding"))%>%
    select(chrom,start,end,strand,geneId)
  gene_pos
}
## extract TE genes from gff
get_te_gene_pos=function(gff){
  te_gene_pos = gff%>%filter(feature=="transposable_element_gene")%>%
    select(chrom,start,end,strand,geneId)
  te_gene_pos
}

## keep only features of interest, remove te gene exons and intro, rename the features and relevel them to give preferred order
gff_to_feature=function(gff){
    feature=gff%>%filter(feature%in%c("gene","exon","five_prime_UTR",
                                      "three_prime_UTR","intron","CDS",
                                      "transposable_element","pseudogene",
                                      "transposon_fragment",
                                      "transposable_element_gene"))%>%
      mutate(feature=gsub("five_prime_UTR","5'UTR",feature))%>%
      mutate(feature=gsub("three_prime_UTR","3'UTR",feature))%>%
      mutate(feature=gsub("transposable_element$","TE",feature))%>%
      mutate(feature=gsub("transposon_fragment","TE fragment",feature))%>%
      mutate(feature=gsub("transposable_element_gene","TE gene",feature))%>%
      mutate(feature=factor(feature,levels=c("gene","5'UTR","exon","intron",
                                             "CDS","3'UTR","TE","TE fragment",
                                             "TE gene","pseudogene")))
    te_gene_id=filter(feature,feature=="TE gene")$geneId
    feature = feature%>%filter((feature%in%c("exon","intro"))+geneId%in%te_gene_id<2)
    feature
}

## keep only rna features of interest, relevel them to give preferred order
gff_to_rnatype=function(gff){
  feature=gff%>%filter(feature%in%c("mRNA","snoRNA","ncRNA","snRNA","rRNA",
                                    "lnc_RNA","miRNA","antisense_lncRNA",
                                    "antisense_RNA","tRNA",
                                    "pseudogenic_tRNA"))%>%
    mutate(feature=factor(feature,levels=c("mRNA","snoRNA","ncRNA","snRNA",
                                           "rRNA","lnc_RNA","miRNA",
                                           "antisense_lncRNA","antisense_RNA",
                                           "tRNA","pseudogenic_tRNA")))
  feature
}  

## how to set range when plotting without outliers  
find_outlier_range <- function(x) {
  return( quantile(x, .75,na.rm=T) + 1.5*IQR(x,na.rm=T))
}

