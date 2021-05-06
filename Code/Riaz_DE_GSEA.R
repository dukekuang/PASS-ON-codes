library(here)
library(openxlsx)                  
source('code/DE_GSEA_function.R') ##DE analysis

here::here()
################################################################################
#load Riaz Data (ipi-PROG-NAIVE Sample info & RNA-seq Count)
load('data/Riaz_data.RData')

Sample.info <- Riaz_data$Samples
count <- Riaz_data$Count

#Response vs. Non Response Analysis condition on Pre/On therapy
Sample.PreOn.split <-  split(Sample.info,Sample.info$PreOn)
Sample.pre <- Sample.PreOn.split$Pre
Sample.on <- Sample.PreOn.split$On

#Time Biopsy PRE vs.ON analysis Condition on Response/Non Response Sample
Sample.RNR.split <-  split(Sample.info,Sample.info$Resp_NoResp)
Sample.R <- Sample.RNR.split$Response
Sample.NR <- Sample.RNR.split$No_Response

################################################################################
#-------------------------  DE analysis  ---------------------------------------
################################################################################
set.seed(12345)
#Pre condition, DE RespANLS(R vs. NR) 
Riaz.Pre.DE <- getDEGs_Resp(data = count,
                            sample.info = Sample.pre,
                            filename='data/Riaz_PRE_Sample_RespANLS_DE',
                            rnk=T) 

#On condition, DE RespANLS(R vs. NR) 
Riaz.On.DE <- getDEGs_Resp(data = count,
                           sample.info = Sample.on,
                           filename='data/Riaz_ON_Sample_RespANLS_DE',
                           rnk=T) 

#Response condition, Time Biopsy Pre vs. On DE
Riaz.R.DE <- getDEGs_Time(data = count,
                          sample.info = Sample.R,
                          filename = 'data/Riaz_Resp_Sample_TimeANLS_DE',
                          rnk=T) 
#Non Response condition, Time Biopsy Pre vs. On DE
Riaz.NR.DE <- getDEGs_Time(data = count,
                           sample.info = Sample.NR,
                           filename = 'data/Riaz_NoResp_Sample_TimeANLS_DE',
                           rnk=T)

################################################################################
#-------------------------     GSEA     ----------------------------------------
################################################################################
set.seed(779948250)
#get Pathways
reactome <- gmtPathways("data/c2.cp.reactome.v7.1.symbols.gmt")
str(head(reactome),2)

#GSEA @ Pre RespANLS(R vs NR) DE
Riaz.Pre.GSEA <-  get_fgsea(rnk_file = 'data/Riaz_PRE_Sample_RespANLS_DE.rnk', pathways=reactome)

#GSEA @ On RespANLS(R vs NR) DE
Riaz.On.GSEA <-  get_fgsea(rnk_file = 'data/Riaz_ON_Sample_RespANLS_DE.rnk', pathways=reactome)

#GSEA @ Response TimeANLS(Pre vs On) DE
Riaz.R.GSEA <-  get_fgsea(rnk_file = 'data/Riaz_Resp_Sample_TimeANLS_DE.rnk', pathways=reactome)

#GSEA @ Non-Response TimeANLS(Pre vs On) DE
Riaz.NR.GSEA <-  get_fgsea(rnk_file = 'data/Riaz_NoResp_Sample_TimeANLS_DE.rnk', pathways=reactome)

#################################################################################################################################################################################################
#From here: This part code no need to publish !!!!!!!
#-----------------------    Filter DE   ----------------------------------------
#Cutoff: |logFoldChange| > 1 && P.adj < 0.05
################################################################################
#Filter Pre/On condition, DE RespANLS(R vs. NR) 
Riaz.Pre.DE.FL <- filter_DE(Riaz.Pre.DE, logFC=1, p.adj = 0.05)
dim(Riaz.Pre.DE.FL) #[1] 190   7

Riaz.On.DE.FL <- filter_DE(Riaz.On.DE, logFC=1, p.adj = 0.05)
dim(Riaz.On.DE.FL) #[1] 1078    7

#Filter Response/Non-Response condition, DE TimeANLS(Pre vs. ON) 
Riaz.R.DE.FL <- filter_DE(Riaz.R.DE, logFC=1, p.adj = 0.05)
dim(Riaz.R.DE.FL) #[1] 60  7

Riaz.NR.DE.FL <- filter_DE(Riaz.NR.DE, logFC=1, p.adj = 0.05)
dim(Riaz.NR.DE.FL) #[1] 0 7

filtered_DEs <- list('Riaz_PRE_Sample_RespANLS_Filter_DE' = Riaz.Pre.DE.FL, 
                     'Riaz_ON_Sample_RespANLS_Filter_DE' = Riaz.On.DE.FL,
                     'Riaz_Resp_Sample_TimeANLS_Filter_DE' = Riaz.R.DE.FL,
                     'Riaz_NoResp_Sample_TimeANLS_Filter_DE' = Riaz.NR.DE.FL
                     )
output_filter_DE(filtered_DEs,'output/Filter_DE/')


################################################################################
#-----------------------    Filter GSEA   --------------------------------------
#select pathway from fGSEA result
#ES>0, padj < 0.05, order(NES)
################################################################################
#Filter Pre/On condition, GSEA RespANLS(R vs. NR) 
Riaz.Pre.GSEA.FL <- filter_GSEA(Riaz.Pre.GSEA, ES.Value=0, p.adj = 0.05)
dim(Riaz.Pre.GSEA.FL) #[1] 98  8

Riaz.On.GSEA.FL <- filter_GSEA(Riaz.On.GSEA, ES.Value=0, p.adj = 0.05)
dim(Riaz.On.GSEA.FL) #[1] 111   8

#Filter Response/Non-Response condition, GSEA TimeANLS(Pre vs. ON) 
Riaz.R.GSEA.FL <- filter_GSEA(Riaz.R.GSEA, ES.Value=0, p.adj = 0.05)
dim(Riaz.R.GSEA.FL) #[1] 110   8

Riaz.NR.GSEA.FL <- filter_GSEA(Riaz.NR.GSEA, ES.Value=0, p.adj = 0.05)
dim(Riaz.NR.GSEA.FL) #[1] 87  8

filtered_GSEAs <- list('Riaz_PRE_Sample_RespANLS_Filter_GSEA' = Riaz.Pre.GSEA.FL, 
                       'Riaz_ON_Sample_RespANLS_Filter_GSEA' = Riaz.On.GSEA.FL,
                       'Riaz_Resp_Sample_TimeANLS_Filter_GSEA' = Riaz.R.GSEA.FL,
                       'Riaz_NoResp_Sample_TimeANLS_Filter_GSEA' = Riaz.NR.GSEA.FL
  )
output_filter_GSEA(filtered_GSEAs,'output/Filter_GSEA/')

################################################################################
#--------------------   DE Volcano Plot   --------------------------------------
################################################################################
#Pre condition, DE RespANLS(R vs. NR) 
DE_Volcano(DEGs_res = Riaz.Pre.DE,
           p.value = 0.05, FC.value = 1,
           plot_title = 'Riaz et al. Pre-treatment Samples (R vs. NR)',
           pointsize = 2.5, labsize = 4.5,
           w = 10, h = 10,
           no_selectLab = 30,
           output_file = 'output/Riaz_RespANLS_PRE_DEplot2.tiff')

#On condition, DE RespANLS(R vs. NR) 


DE_Volcano(DEGs_res = Riaz.On.DE,
           p.value = 0.05, FC.value = 1,
           plot_title = 'Riaz et al. On-treatment Samples (R vs. NR)',
           pointsize = 2.5, labsize = 4.5, 
           xlimit = c(-12,15), ylimit = c(0,30),
           w = 10, h = 10,
           no_posLab = 35, no_neglab = 15,
           output_file = 'output/Riaz_RespANLS_ON_DEplot2.tiff')

#Response condition, Time Biopsy Pre vs. On DE
DE_Volcano(DEGs_res = Riaz.R.DE,
           p.value = 0.05, FC.value = 1,
           #selectlab_genes = unlist(ipiPROG_NAIVE_TimeANLS_sig),
           plot_title = 'Riaz et al. Responders\n(Pre-treatment Time vs. On-treatment Time)',
           pointsize = 2.5, labsize = 4.5, 
           xlimit = c(-6,6),ylimit = c(0,7.5),
           w = 10, h = 10,
           no_posLab = 35, no_neglab = 20,
           output_file = 'output/Riaz_TimeANLS_Response_DEplot2.tiff')

#Non Response condition, Time Biopsy Pre vs. On DE
DE_Volcano(DEGs_res = Riaz.NR.DE,
           p.value = 0.05, FC.value = 1,
           #selectlab_genes = unlist(ipiPROG_NAIVE_TimeANLS_sig),
           plot_title = 'Riaz et al. Non-responders\n(Pre-treatment Time vs. On-treatment Time)',
           pointsize = 2.5, labsize = 4.5, 
           xlimit = c(-6,6),ylimit = c(0,6),
           w = 10, h = 10,
           no_posLab = 45, no_neglab = 26,
           output_file = 'output/Riaz_TimeANLS_NonResponse_DEplot2.tiff')

################################################################################
#--------   Generate Signature from GSEA for ssGSEA Use   ----------------------
#-------------   Plot GSEA table for signatures   ------------------------------
################################################################################
#Generate Signatures from filtered GSEA result ----------
#RespANLS signatures @ Pre Smaples (Top 15 Pathway)
RespANLS.Pre.sig  <- Riaz.Pre.GSEA.FL[1:15,]

#RespANLS signatures @ On Samples
RespANLS.On.sig <- Riaz.On.GSEA.FL[1:15,]

#TimeANLS signatures for Pre/On Samples (Top 15 Pathway in Response, exclude Non-reponse)
TimeANLS.R.NR.inter <- intersect(Riaz.R.GSEA.FL$pathway, Riaz.NR.GSEA.FL$pathway)
Riaz.R.GSEA.FL.unique <- Riaz.R.GSEA.FL[-which(Riaz.R.GSEA.FL$pathway %in% TimeANLS.R.NR.inter),]
TimeANLS.sigs <- Riaz.R.GSEA.FL.unique[1:15,]

#Prepare signature for ssGSEA
Pathway.Sigs <- list('PASS_PRE' = RespANLS.Pre.sig,
                     'PASS_ON' = RespANLS.On.sig,
                     'TimeANLS_Sig' = TimeANLS.sigs)
save(Pathway.Sigs, file='output/Pathway_Singatures.Rdata')

#plot GSEA Table-----------------------------------------
names(reactome) <- gsub('REACTOME_','',names(reactome))
#RespANLS signatures @ Pre Samples (Top 15 Pathway)
Riaz.Pre.DE.rnk <- read_rnk('data/Riaz_PRE_Sample_RespANLS_DE.rnk')
tiff(file='output/RespANLS_Pre_Sig.tiff',width=5,height = 5,unit='in',res=300,compression='lzw')
plotGseaTable(reactome[RespANLS.Pre.sig$pathway], Riaz.Pre.DE.rnk, Riaz.Pre.GSEA, colwidths=c(5, 3.5, 0.5, 0.5, 0.5, 0.5))
dev.off()

#RespANLS signatures @ On Samples (Top 15 Pathway)
Riaz.On.DE.rnk <- read_rnk('data/Riaz_ON_Sample_RespANLS_DE.rnk')
tiff(file='output/RespANLS_ON_Sig.tiff',width=5,height = 5,unit='in',res=300,compression='lzw')
plotGseaTable(reactome[RespANLS.On.sig$pathway], Riaz.On.DE.rnk, Riaz.On.GSEA, colwidths=c(5, 3.5, 0.5, 0.5, 0.5, 0.5))
dev.off()

#TimeANLS signatures @ Response Samples
Riaz.R.DE.rnk <- read_rnk('data/Riaz_Resp_Sample_TimeANLS_DE.rnk')
tiff(file='output/TimeANLS_Sig.tiff',width=5,height = 5,unit='in',res=300,compression='lzw')
plotGseaTable(reactome[TimeANLS.sigs$pathway], Riaz.R.DE.rnk, Riaz.R.GSEA, colwidths=c(5, 3.5, 0.5, 0.5, 0.5, 0.5))
dev.off()

################################################################################
#-------------   Plot Pahtway   ------------------------------
################################################################################
Pathways_Enrichment_plot <- function(selected_pathways_names, rnk, plot_folder, plot_width, plot_height){
  for(i in 1:length(selected_pathways_names)){
    print(i)
    single_pathway <- selected_pathways_names[i]
    p <- plotEnrichment(reactome[[single_pathway]],rnk) + labs(title=single_pathway) + labs(title=single_pathway) + theme(plot.title = element_text(size = 6.5))
    ggsave(p, file=file.path(plot_folder,paste0(single_pathway,'.tiff',sep='')), width = plot_width, height = plot_height, units = "in", dpi=200, compression='lzw')
    print(single_pathway)
  }
}


Pathways_Enrichment_plot(selected_pathways_names = RespANLS.Pre.sig$pathway, 
                         rnk = Riaz.Pre.DE.rnk, 
                         plot_folder = 'data/RespANLS.Pre',
                         plot_width = 5,
                         plot_height = 3)

#RespANLS_ON_Sigs
Pathways_Enrichment_plot(selected_pathways_names = RespANLS.On.sig$pathway, 
                         rnk = Riaz.On.DE.rnk, 
                         plot_folder = 'data/RespANLS.On',
                         plot_width = 5,
                         plot_height = 3)

#TimeANLS_Sigs
Pathways_Enrichment_plot(selected_pathways_names = TimeANLS.sigs$pathway, 
                         rnk = Riaz.R.DE.rnk, 
                         plot_folder = 'data/TimeANLS.sigs',
                         plot_width = 5,
                         plot_height = 3)


