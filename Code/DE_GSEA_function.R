library(org.Hs.eg.db)
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(annotate)
library("Biobase")
library("matrixStats")
library("DESeq2")
library("BiocParallel")
register(MulticoreParam(30))
library(fgsea)
library(ggplot2)
library(data.table)
library(ggrepel)
library(gridExtra)
library(grid)

####################################################################
# Gene symbol search function
ENTREZtoSYMBOL <- function(x){
  vec=mapIds(org.Hs.eg.db,keys=x,column="SYMBOL",keytype="ENTREZID",multiVals="first")
  repl=which(is.na(vec))
  for(y in repl){vec[y]=x[y]}
  return(vec)
}

####################################################################
# Get DE analysis and rnk
# Only compare Response vs. No Response
getDEGs_Resp <- function(data, sample.info, filename, rnk=F, up=T, norm.out=F){
  # compare Response vs. No Response
  # Find the overlapping samples
  inter <- intersect(colnames(data),sample.info$Sample)
  info <- sample.info[sample.info$Sample %in% inter,]
  data <- data[,match(inter,colnames(data))]

  # Create the DESeq2 object
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  ebg <- exonsBy(txdb, by="gene")
  intersection <- intersect(rownames(data),ebg@partitioning@NAMES)
  ebg2 <- ebg[ebg@partitioning@NAMES %in% intersection]

  # sort by ID
  ebg2 <- ebg2[order(names(ebg2))]

  # Sort by gene model order
  data <- data[match(names(ebg2), rownames(data)),]

  # Create DESeq2 object
  ddsMat <- DESeqDataSetFromMatrix(countData = data,
                                   colData = DataFrame(info),
                                   design = ~ 1)
  # Get the gene symbol
  mcols(ddsMat)$symbols <-  mapIds(org.Hs.eg.db,
                                       keys=rownames(ddsMat),
                                       column="SYMBOL",
                                       keytype="ENTREZID",
                                       multiVals="first")
  # Use ALL samples in analysis
  dds.resp <- ddsMat

  # Set No_Response as the reference level
  dds.resp$Resp_NoResp <- factor(dds.resp$Resp_NoResp,levels = c("No_Response", "Response"))

  # Drop patient ID level which is not included
  dds.resp$Sample <- factor(dds.resp$Sample)

  # Set the design with controling patient
  design(dds.resp)<- ~ Resp_NoResp

  # Calculate the sample size factor
  dds.resp <- estimateSizeFactors(dds.resp)
  
  # Run the DEG analysis
  dds.resp <- DESeq(dds.resp)
  norm.ds <- counts(dds.resp,normalized = TRUE)
  colnames(norm.ds) <- info$Sample

  print(resultsNames(dds.resp))
  
  # Get the result
  res <-  results(dds.resp, contrast=c("Resp_NoResp","Response","No_Response"))
  res <- na.omit(res)
  res$symbol <- ENTREZtoSYMBOL(row.names(res))
  res <- res[,c(7,1:6)]
  norm.ds <- norm.ds[row.names(res),] #get same order
  rownames(norm.ds) <- res$symbol
  
  #Regulate Up/Down
  if(up){
    idx = sort(res$stat, decreasing = T, index.return=T)$ix
  }else{
    idx = sort(res$stat, decreasing = F, index.return=T)$ix
  }
  write.csv(as.data.frame(res[idx,]),paste(filename,".csv",sep=""),row.names = FALSE)
  
  #output rnk
  if(rnk){
    res_rnk <- res[idx,]
    rnk <- data.frame(Gene=res_rnk$symbol,stat = res_rnk$stat)
    write("#    rnk",paste(filename,".rnk",sep=""))
    write.table(rnk,paste(filename,".rnk",sep=""),sep="\t",
                row.names=F,col.names = F,quote = F,append=T)
  }
  if(norm.out){
    write.csv(as.matrix(norm.ds),paste(filename,"_norm.csv",sep=""))
  }
  print('finish!')
  return(res[idx,])
}

getDEGs_Time <- function(data, sample.info, filename, rnk=F, up=T, norm.out=F){
  # compare Response vs. No Response
  # Find the overlapping samples
  inter <- intersect(colnames(data),sample.info$Sample)
  info <- sample.info[sample.info$Sample %in% inter,]
  data <- data[,match(info$Sample,colnames(data))]

  # Create the DESeq2 object
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  ebg <- exonsBy(txdb, by="gene")
  intersection <- intersect(rownames(data),ebg@partitioning@NAMES)
  ebg2 <- ebg[ebg@partitioning@NAMES %in% intersection]

  # sort by ID
  ebg2 <- ebg2[order(names(ebg2))]

  # Sort by gene model order
  data <- data[match(names(ebg2), rownames(data)),]

  # Create DESeq2 object
  ddsMat <- DESeqDataSetFromMatrix(countData = data,
                                   colData = DataFrame(info),
                                   design = ~ 1)
  # Get the gene symbol
  mcols(ddsMat)$symbols <-  mapIds(org.Hs.eg.db,
                                   keys=rownames(ddsMat),
                                   column="SYMBOL",
                                   keytype="ENTREZID",
                                   multiVals="first")
  # Use ALL samples in analysis
  dds.resp <- ddsMat

  # Set PRE as the reference level (PreOn)
  dds.resp$PreOn <- factor(dds.resp$PreOn,levels = c("Pre", "On"))

  # Drop patient ID level which is not included
  dds.resp$Sample <- factor(dds.resp$Sample)

  # factor Patient ID
  dds.resp$PatientID <- factor(dds.resp$PatientID)

  # Set the design with controling patient
  design(dds.resp)<- ~ PatientID + PreOn

  # Calculate the sample size factor
  dds.resp <- estimateSizeFactors(dds.resp)
  
  # Run the DEG analysis
  dds.resp <- DESeq(dds.resp)
  norm.ds <- counts(dds.resp,normalized = TRUE)
  colnames(norm.ds) <- info$Sample

  print(resultsNames(dds.resp))
  
  # Get the result
  res <-  results(dds.resp, contrast=c("PreOn","On","Pre"))
  res <- na.omit(res)
  res$symbol <- ENTREZtoSYMBOL(row.names(res))
  res <- res[,c(7,1:6)]
  norm.ds <- norm.ds[row.names(res),] #get same order
  rownames(norm.ds) <- res$symbol
  
  #Regulate Up/Down
  if(up){
    idx = sort(res$stat, decreasing = T, index.return=T)$ix
  }else{
    idx = sort(res$stat, decreasing = F, index.return=T)$ix
  }
  write.csv(as.data.frame(res[idx,]),paste(filename,".csv",sep=""),row.names = FALSE)
  
  #output rnk
  if(rnk){
    res_rnk <- res[idx,]
    rnk <- data.frame(Gene=res_rnk$symbol,stat = res_rnk$stat)
    write("#    rnk",paste(filename,".rnk",sep=""))
    write.table(rnk,paste(filename,".rnk",sep=""),sep="\t",
                row.names=F,col.names = F,quote = F,append=T)
  }
  if(norm.out){
    write.csv(as.matrix(norm.ds),paste(filename,"_norm.csv",sep=""))
  }
  print('finish!')
  return(res[idx,])
}

####################################################################
# GSEA (use fgsea library)
fast_gsea<- function(rnk,pathways,min=15, max=500, np=10000,filename=''){
  fgseaRes<- fgsea(pathways, rnk, minSize=min, maxSize=max, nperm=np)
  fgseaRes <- fgseaRes[order(NES,decreasing=TRUE), ]
  fwrite(data.frame(fgseaRes), file=filename)
  return(fgseaRes)
}

get_fgsea<-function(rnk_file,pathways=reactome){
  df <- read.table(rnk_file,colClasses = c("character", "numeric"))
  colnames(df) <- c('symbol','stat')
  df.2 <- setNames(df$stat, df$symbol) #reshape rank
  str(df.2)
  filename <- unlist(strsplit(rnk_file,'\\.'))[1]
  print(filename)
  df.gsea <- fast_gsea(rnk = df.2, 
                       pathways = pathways,
                       filename = paste(filename,'_fgsea.csv',sep =''))
  return(df.gsea)
}

####################################################################
filter_DE <- function(DEGs,logFC,p.adj){
  t <- DEGs[abs(DEGs$log2FoldChange) > logFC,]
  t <- t[t$padj < p.adj,]
  return(t)
}

output_filter_DE <- function(filter_list,file_path){
  for(i in 1:length(filter_list)){
    output_DE <- filter_list[[i]]
    csv_path <- file.path(file_path,paste(names(filter_list)[i],'.csv',sep=''))
    write.csv(output_DE,file=csv_path, row.names = F)
  }
}

####################################################################
filter_GSEA <- function(GSEA,ES.Value,p.adj){
  t <- GSEA[GSEA$ES > ES.Value,]
  t <- t[t$padj < p.adj,]
  t <- t[order(t$NES,decreasing = TRUE),]
  t$pathway <- gsub('REACTOME_', '', t$pathway)
  return(t)
}

output_filter_GSEA <- function(filter_list,file_path){
  for(i in 1:length(filter_list)){
    output_GSEA <- data.frame(filter_list[[i]])
    csv_path <- file.path(file_path,paste(names(filter_list)[i],'.csv',sep=''))
    fwrite(output_GSEA,file=csv_path)
  }
}

####################################################################


####################################################################
#Probably No Need to Publish
#Plot DE & GSEA
####################################################################
##################
#DE Volcano plot
DE_Volcano <- function(DEGs_res, p.value = 0.05, FC.value = 1, xlimit = c(-10,10), ylimit = c(0,20), plot_title,pointsize, labsize, output_file, w, h, no_posLab, no_neglab){
  require(EnhancedVolcano)
  rownames(DEGs_res) <- DEGs_res$symbol
  
  FC <- FC.value
  p <- p.value

  keyvals <- rep('grey75', nrow(DEGs_res))
  names(keyvals) <- rep('NS', nrow(DEGs_res))

  keyvals[which(DEGs_res$log2FoldChange < -FC & DEGs_res$pvalue < p)] <- "#1F78B4" 
  names(keyvals)[which(DEGs_res$log2FoldChange  < -FC & DEGs_res$pvalue < p)] <- 'Signif. down-regulated'
  keyvals[which(DEGs_res$log2FoldChange > FC & DEGs_res$pvalue < p)] <- "#E31A1C"
  names(keyvals)[which(DEGs_res$log2FoldChange > FC & DEGs_res$pvalue < p)] <- 'Signif. up-regulated'

  selectLab.vec <- DEGs_res[which(keyvals != 'grey75'),]
  selectLab.vec.pos <- selectLab.vec[order(selectLab.vec$log2FoldChange,decreasing = T),]$symbol[1:no_posLab]
  selectLab.vec.neg <- selectLab.vec[order(selectLab.vec$log2FoldChange,decreasing = F),]$symbol[1:no_neglab]
    
  EV <- EnhancedVolcano(DEGs_res,
                        lab = rownames(DEGs_res),
                        selectLab = c(selectLab.vec.pos,selectLab.vec.neg),
                        pCutoff = p,
                        FCcutoff = FC,
                        cutoffLineType = 'blank',
                        title = plot_title,
                        subtitle = NULL,
                        xlim = xlimit,
                        ylim = ylimit,
                        x = 'log2FoldChange',
                        y = 'pvalue',
                        pointSize = pointsize, #2.5
                        labSize = labsize, #4.5
                        caption = paste("FC cutoff ", FC," ; p-value cutoff ", p,sep=''),
                        captionLabSize = 15,
                        colCustom = keyvals,
                        colAlpha = 0.75,
                        legendPosition = 'top',
                        legendLabSize = 15,
                        legendIconSize = 5.0,
                        drawConnectors = FALSE,  
                        widthConnectors = 0.5,    
                        colConnectors = 'grey50',         
                        gridlines.major = TRUE,
                        gridlines.minor = FALSE,     
                        border = 'partial',
                        borderWidth = 1.5,
                        borderColour = 'black'
                        )
  tiff(output_file, w, h, units = 'in', compression = 'lzw', res=150)
  print(EV)
  dev.off()
}

##################
#GSEA table plot
read_rnk <- function(rank_path){
  rnk <- read.table(rank_path,colClasses = c("character", "numeric"))
  colnames(rnk) <- c('symbol','stat')
  rnk2 <- setNames(rnk$stat, rnk$symbol) #reshape rank
  str(rnk2)
  return(rnk2)
}

plotGseaTable <- function(pathways, stats, fgseaRes, gseaParam=1,
                          #colwidths=c(5, 3, 0.8, 1.2, 1.2),
                          colwidths=c(5, 3, 0.8, 1.2, 1.2, 1.2),
                          render=TRUE) {
  
  fgseaRes$pathway <- gsub('REACTOME_', '', fgseaRes$pathway)
  rnk <- rank(-stats)
  ord <- order(rnk)
  
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
  statsAdj <- statsAdj / max(abs(statsAdj))
  
  pathways <- lapply(pathways, function(p) { #View(pathways)
    unname(as.vector(na.omit(match(p, names(statsAdj))))) #unname(na.omit(match(pathways[[1]], names(statsAdj))))
  })
  
  # fixes #40
  pathways <- pathways[sapply(pathways, length) > 0]
  
  ps <- lapply(names(pathways), function(pn) { #names(pathways)[1]
    p <- pathways[[pn]] 
    annotation <- fgseaRes[match(pn, fgseaRes$pathway), ]
    if(pn == "REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS"){
      pn = "REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR\n_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE\n_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS"
    }
    if(pn == "IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL"){
      pn = "IMMUNOREGULATORY_INTERACTIONS\n_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL"
    }
    if(pn == "SYNTHESIS_OF_BILE_ACIDS_AND_BILE_SALTS_VIA_7ALPHA_HYDROXYCHOLESTEROL"){
      pn = "SYNTHESIS_OF_BILE_ACIDS_AND_BILE_SALTS\n_VIA_7ALPHA_HYDROXYCHOLESTEROL"
    }
    if(pn == "ANTIGEN_ACTIVATES_B_CELL_RECEPTOR_BCR_LEADING_TO_GENERATION_OF_SECOND_MESSENGERS"){
      pn = "ANTIGEN_ACTIVATES_B_CELL_RECEPTOR_BCR\n_LEADING_TO_GENERATION_OF_SECOND_MESSENGERS"
    }
    if(pn == "TNF_RECEPTOR_SUPERFAMILY_TNFSF_MEMBERS_MEDIATING_NON_CANONICAL_NF_KB_PATHWAY"){
      pn = "TNF_RECEPTOR_SUPERFAMILY_TNFSF_MEMBERS\n_MEDIATING_NON_CANONICAL_NF_KB_PATHWAY"
    }
    if(pn == "CASPASE_ACTIVATION_VIA_EXTRINSIC_APOPTOTIC_SIGNALLING_PATHWAY"){
      pn = "CASPASE_ACTIVATION_VIA_EXTRINSIC\n_APOPTOTIC_SIGNALLING_PATHWAY"
    }
    if(pn == "CASPASE_ACTIVATION_VIA_DEATH_RECEPTORS_IN_THE_PRESENCE_OF_LIGAND"){
      pn = "CASPASE_ACTIVATION_VIA_DEATH_RECEPTORS\n_IN_THE_PRESENCE_OF_LIGAND"
    }
    
    list(
      textGrob(pn, just="right", x=unit(0.95, "npc"),check.overlap=TRUE,gp = gpar(fontsize = 4.5)), #pathway names
      ggplot() +
        geom_segment(aes(x=p, xend=p,
                         y=0, yend=statsAdj[p]),
                     size=0.35, lineend = "round", color="blue") +  #Gene ranks; width size=1; lineend = "round"
        scale_x_continuous(limits=c(0, length(statsAdj)),
                           expand=c(0, 0)) +
        scale_y_continuous(limits=c(-1, 1),
                           expand=c(0, 0)) +
        xlab(NULL) + ylab(NULL) +
        theme(panel.background = element_blank(),
              axis.line=element_blank(),
              axis.text=element_blank(),
              axis.ticks=element_blank(),
              panel.grid = element_blank(),
              axis.title=element_blank(),
              plot.margin = rep(unit(0,"null"),4),
              panel.spacing = rep(unit(0,"null"),4) #panel.border = element_rect(linetype = "dashed", fill = NA)
        ),
      textGrob(sprintf("%.2f", annotation$NES),check.overlap=TRUE,gp = gpar(fontsize = 4.5)),
      textGrob(sprintf("%.1e", annotation$pval),check.overlap=TRUE,gp = gpar(fontsize = 4.5)), #Change
      textGrob(sprintf("%.1e", annotation$padj),check.overlap=TRUE,gp = gpar(fontsize = 4.5)),
      textGrob(annotation$size,check.overlap=TRUE, gp = gpar(fontsize = 4.5)) #size – size of the pathway after removing genes not present in 'names(stats)'
    )
  })
  
  
  rankPlot <-
    ggplot() +
    geom_blank() +
    scale_x_continuous(limits=c(0, 20000), #length(statsAdj) #number of genes
                       expand=c(0, 0)) +
    scale_y_continuous(limits=c(-1, 1),
                       expand=c(0, 0)) +
    xlab(NULL) + ylab(NULL) +
    theme(panel.background = element_blank(),
          axis.text.x = element_text(size=6),
          axis.line=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid = element_blank(),
          axis.title=element_blank(),
          plot.margin = unit(c(0,0,0,0), "npc"),#c(0,0,0.5,0)
          panel.spacing = unit(c(0,0,0,0), "npc")
    )
  
  grobs <- c(
    list(textGrob("Pathway", x=unit(0.80, "npc"),check.overlap=TRUE,just="right",gp = gpar(fontsize = 4.5))),
    lapply(c("Gene ranks", "NES", "P-value","FDR","Size"), function(x)textGrob(x,x=unit(0.5, "npc"),check.overlap=TRUE,gp = gpar(fontsize = 4.5))),
    #lapply(c("Gene ranks", "NES", "pval", "padj"), function(x)textGrob(x,x=unit(0.5, "npc"),check.overlap=TRUE,gp = gpar(fontsize = 9))),
    unlist(ps, recursive = FALSE),
    list(nullGrob(),
         rankPlot,
         nullGrob(),
         nullGrob(),
         nullGrob()))
  
  # not drawing column if corresponding colwidth is set to zero
  grobsToDraw <- rep(colwidths != 0, length(grobs)/length(colwidths))
  
  p <- arrangeGrob(grobs=grobs[grobsToDraw],
                   ncol=sum(colwidths != 0),
                   widths=colwidths[colwidths != 0])
  
  if (render) {
    grid.draw(p)
  } else {
    p
  }
}
