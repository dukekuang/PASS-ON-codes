library(here)
library(splines)
library(maxstat)
library(survival)
library(survminer)
library(RColorBrewer)
source('code/FE_Model_Functions.R') #Feature Engineering & Modeling
source('code/Model_Development/Model_Devep_Fun.R')

prob2odds <- function(prob){
  odds = prob/(1-prob)
  return(odds)
}
########################################
##### Load Sample Information
########################################
#Riaz
load('data/Riaz_data.RData')
#Gide
load('data/Gide_data.RData')
#Lee
load('data/Lee_data.RData')
#MGH
load('data/MGH_ON_data.RData')

########################################
##### Load Sample Log Odds
########################################
load('output/PASS_ON/PASS_ON_SigScore.Rdata')

########################################
#### Riaz Surv Analysis
########################################
Riaz.res <- Riaz_data$Samples[match(rownames(PASS_ON_ss$Riaz), Riaz_data$Samples$Sample),]
Riaz.res$SS <- prob2odds(PASS_ON_ss$Riaz) #Get odds instend of Probability
Riaz.res.surv <- merge(Riaz.res, Riaz_data$Surv.info, by='PatientID')

#The status indicator, normally 0=alive, 1=dead.
Riaz.res.surv$PFS_SOR2 <- ifelse(Riaz.res.surv$PFS_SOR == 0,1,0)
Riaz.res.surv$OS_SOR2 <- ifelse(Riaz.res.surv$OS_SOR == 0,1,0)

Riaz.res.surv$GruopSS <- ifelse(Riaz.res.surv$SS >= mean(Riaz.res.surv$SS),'High_model_score','Low_model_score')
Riaz.res.surv$GruopSS <- as.factor(Riaz.res.surv$GruopSS)

#OS plot
Riaz.OS.fit <- survfit(Surv(OS, OS_SOR2) ~ GruopSS, data = Riaz.res.surv)
tiff(file="output/PASS_ON/Riaz_PASS_ON_OS_odds_mean.tiff",width=7, height=5,unit='in',compression = 'lzw',res=150)
ggsurvplot(Riaz.OS.fit,
           pval = TRUE, 
           pval.method = FALSE,
           risk.table = TRUE, # Add risk table
           risk.table.col = 'black',
           risk.table.y.text = TRUE,
           risk.table.y.text.col = TRUE,
           legend.labs = c("High (n=15)", "Low (n=33)"),
           legend.title = 'Signature Score',
           ggtheme = theme_bw(), # Change ggplot2 theme
           tables.theme = theme_void(),
           palette = c("#E31A1C","#1F78B4"),
           legend = c(0.12,0.20),#"right",
           pval.coord = c(400,0.20),
           pval.method.coord = c(750,0.60),
           title='Overall Survival of Riaz et al. On-treatment Samples'
           )
dev.off()

#PFS plot
Riaz.PFS.fit <- survfit(Surv(PFS, PFS_SOR2) ~ GruopSS, data = Riaz.res.surv)
tiff(file="output/PASS_ON/Riaz_PASS_ON_PFS_odds_mean.tiff",width=7, height=5,unit='in',compression = 'lzw',res=150)
ggsurvplot(Riaz.PFS.fit,
           pval = TRUE, 
           pval.method = FALSE,
           risk.table = TRUE, 
           risk.table.col = 'black',
           risk.table.y.text = TRUE,
           risk.table.y.text.col = TRUE,
           legend.labs = c("High (n=15)", "Low (n=33)"),
           legend.title = 'Signature Score',
           ggtheme = theme_bw(), # Change ggplot2 theme
           tables.theme = theme_void(),
           palette = c("#E31A1C","#1F78B4"), #palette = c("red", "blue")
           legend = c(0.85,0.85),#"right",
           title='Progression Free Survival of Riaz et al. On-treatment Samples',
           pval.coord = c(450,0.75),
           pval.method.coord = c(750,0.60),
           xlim = c(0, 1050)
)
dev.off()

########################################
#### Gide Surv Analysis
########################################
Gide.res <- Gide_data$Samples[match(rownames(PASS_ON_ss$Gide), Gide_data$Samples$Sample),]
Gide.res$SS <- prob2odds(PASS_ON_ss$Gide) #Get odds instend of Probability

#The status indicator, normally 0=alive, 1=dead.
Gide.res$PFS_event <- ifelse(Gide.res$Progressed == 'No',0,1)
Gide.res$OS_event <- ifelse(Gide.res$Last.Followup.Status == 'Alive',0,1)
colnames(Gide.res)[5] <- 'PFS'
colnames(Gide.res)[6] <- 'OS'

Gide.res$GruopSS <- ifelse(Gide.res$SS >= mean(Gide.res$SS),'High_model_score','Low_model_score')
Gide.res$GruopSS <- as.factor(Gide.res$GruopSS)

#OS plot
Gide.OS.fit<-survfit(Surv(OS, OS_event) ~ GruopSS, data = Gide.res)
tiff(file="output/PASS_ON/Gide_PASS_ON_OS_odds_mean.tiff",width=7, height=5,unit='in',compression = 'lzw',res=150)
ggsurvplot(Gide.OS.fit,
           pval = TRUE, 
           pval.method = FALSE,
           risk.table = TRUE, # Add risk table
           risk.table.col = 'black',
           risk.table.y.text = TRUE,
           risk.table.y.text.col = TRUE,
           legend.labs = c("High (n=7)", "Low (n=11)"),
           legend.title = 'Signature Score',
           ggtheme = theme_bw(), # Change ggplot2 theme
           tables.theme = theme_void(),
           palette = c("#E31A1C","#1F78B4"),
           legend = c(0.12,0.20),#"right",
           pval.coord = c(400,0.20),
           title='Overall Survival of Gide et al. On-treatment Samples'
          )
dev.off()

#PFS plot
Gide.PFS.fit<-survfit(Surv(PFS, PFS_event) ~ GruopSS, data = Gide.res)
tiff(file="output/PASS_ON/Gide_PASS_ON_PFS_odds_mean.tiff",width=7, height=5,unit='in',compression = 'lzw',res=150)
ggsurvplot(Gide.PFS.fit,
           pval = TRUE, 
           pval.method = FALSE,
           risk.table = TRUE, 
           risk.table.col = 'black',
           risk.table.y.text = TRUE,
           risk.table.y.text.col = TRUE,
           legend.labs = c("High (n=7)", "Low (n=11)"),
           legend.title = 'Signature Score',
           ggtheme = theme_bw(), # Change ggplot2 theme
           tables.theme = theme_void(),
           palette = c("#E31A1C","#1F78B4"), #palette = c("red", "blue")
           legend = c(0.85,0.85),#"right",
           title='Progression Free Survival of Gide et al. On-treatment Samples',
           pval.coord = c(550,0.85)
           #pval.method.coord = c(900,0.85)
)
dev.off()

########################################
#### Lee Surv Analysis
########################################。
Lee.res <- Lee_data$On_Samples[match(rownames(PASS_ON_ss$Lee), Lee_data$On_Samples$Pt_ID),]
Lee.res$SS <- prob2odds(PASS_ON_ss$Lee)
#Lee.res$SS <- round(Lee.res$SS,2)

which(Lee.res$Patient.ID == 'Patient 24' & Lee.res$On.treatment.Biopsy == 'Responding')
#[1] 7
which(Lee.res$Patient.ID == 'Patient 31' & Lee.res$On.treatment.Biopsy == 'Innate PROG')
#[1] 12 13 14 15
Lee.res[which(Lee.res$Patient.ID == 'Patient 31' & Lee.res$On.treatment.Biopsy == 'Innate PROG'),]

which(Lee.res$Patient.ID == 'Patient 44' & Lee.res$On.treatment.Biopsy == 'Acquired PROG')
#17
which(Lee.res$Patient.ID == 'Patient 46' & Lee.res$On.treatment.Biopsy == 'Acquired PROG')
#[1] 19 20
Lee.res[which(Lee.res$Patient.ID == 'Patient 46' & Lee.res$On.treatment.Biopsy == 'Acquired PROG'),]

which(Lee.res$Patient.ID == 'Patient 48' & Lee.res$On.treatment.Biopsy == 'Innate PROG')
#[1] 22 23
Lee.res[which(Lee.res$Patient.ID == 'Patient 48' & Lee.res$On.treatment.Biopsy == 'Innate PROG'),]

which(Lee.res$Patient.ID == 'Patient 52' & Lee.res$On.treatment.Biopsy == 'Innate PROG')
#[1] 27 28

Lee.res <- Lee.res[-c(7,13,14,15,17,19,22,27,30,31,33,34),]


Lee.res$GruopSS <- ifelse(Lee.res$SS >= mean(Lee.res$SS),'High_model_score','Low_model_score')
#Lee.res$GruopSS <- ifelse(Lee.res$SS >= 0.5,'High_model_score','Low_model_score')
Lee.res$GruopSS <- as.factor(Lee.res$GruopSS)
Lee.res$OS.days <- Lee.res$`Overall.survival..months.`* 30
Lee.res$OS_stauts <- ifelse(Lee.res$Vital.status == "Alive", 0,1) #The status indicator, normally 0=alive, 1=dead.

#OS plot
Lee.OS.fit<-survfit(Surv(OS.days, OS_stauts) ~ GruopSS, data = Lee.res)
tiff(file="output/PASS_ON/Lee_PASS_ON_OS_odds_mean.tiff",width=7, height=5,unit='in',compression = 'lzw',res=150)
ggsurvplot(Lee.OS.fit,
           pval = TRUE, 
           pval.method = FALSE,
           test.for.trend = FALSE,
           risk.table = TRUE, # Add risk table
           risk.table.col = 'black',
           risk.table.y.text = TRUE,
           risk.table.y.text.col = TRUE,
           legend.labs = c("High (n=11)", "Low (n=12)"),
           legend.title = 'Signature Score',
           ggtheme = theme_bw(), # Change ggplot2 theme
           tables.theme = theme_void(),
           palette = c("#E31A1C","#1F78B4"),
           legend = c(0.12,0.20),#"right",
           pval.coord = c(1700,0.15),
           xlim = c(0, 2200),
           title='Overall Survival of Lee et al. On-treatment Samples'
)
dev.off()

########################################
#### MGH Surv Analysis
########################################。
MGH.res <- MGH_ON_data$Samples[match(rownames(PASS_ON_ss$MGH), MGH_ON_data$Samples$Sample),]
MGH.res$SS <- prob2odds(PASS_ON_ss$MGH)

which(MGH.res$patient == 530 & MGH.res$day == 207)
#[1] 4
which(MGH.res$patient == 42 & MGH.res$day == 73)
#[1] 5 6
which(MGH.res$patient == 'IPIPD1001' & MGH.res$day == 434)
#[1]  8 9
which(MGH.res$patient == 'IPIPD1001' & MGH.res$day == 439)
#[1] 10 11 12

which(MGH.res$patient == 51 & MGH.res$day == 589)
#28
which(MGH.res$patient == 530 & MGH.res$day == 207)
#3

which(MGH.res$patient == 62 & MGH.res$day == 45)
#26

#Patient 39
which(MGH.res$patient == 39 & MGH.res$day == 21)
#[1] 14 15 16 17 18
which(MGH.res$patient == 39 & MGH.res$day == 88)
#[1] 19
MGH.res[which(MGH.res$patient == 39 & MGH.res$day == 21),]

which(MGH.res$patient == 148 & MGH.res$day == 209)
#[1] 21 22
MGH.res[which(MGH.res$patient == 148 & MGH.res$day == 209),]

MGH.res <- MGH.res[-c(3,5,6,8,9,10,11,12,15:19,21,26,28),]

MGH.res$GroupSS <- ifelse(MGH.res$SS>=mean(MGH.res$SS),'High_model_score','Low_model_score')
MGH.res$GroupSS <- as.factor(MGH.res$GroupSS)

#OS plot
MGH.OS.fit<-survfit(Surv(os, censure_os) ~ GroupSS, data = MGH.res)
tiff(file="output/PASS_ON/MGH_PASS_ON_OS_odds_mean.tiff",width=7, height=5,unit='in',compression = 'lzw',res=150)
ggsurvplot(MGH.OS.fit,
           pval = TRUE, 
           pval.method = FALSE,
           risk.table = TRUE, # Add risk table
           risk.table.col = 'black',
           risk.table.y.text = TRUE,
           risk.table.y.text.col = TRUE,
           legend.labs = c("High (n=5)", "Low (n=10)"),
           legend.title = 'Signature Score',
           ggtheme = theme_bw(), # Change ggplot2 theme
           tables.theme = theme_void(),
           palette = c("#E31A1C","#1F78B4"),
           legend = c(0.85,0.85),#"right",
           title='Overall Survival of MGH On-treatment Samples',
           pval.coord = c(0.1,0.1),
           pval.method.coord = c(1500,0.85)
)
dev.off()

#OS plot
MGH.PFS.fit<-survfit(Surv(pfs, censure_pfs) ~ GroupSS, data = MGH.res)
tiff(file="output/PASS_ON/MGH_PASS_ON_PFS_odds_mean.tiff",width=7, height=5,unit='in',compression = 'lzw',res=150)
ggsurvplot(MGH.PFS.fit,
           pval = TRUE, 
           pval.method = FALSE,
           risk.table = TRUE, 
           risk.table.col = 'black',
           risk.table.y.text = TRUE,
           risk.table.y.text.col = TRUE,
           legend.labs = c("High (n=5)", "Low (n=10)"),
           legend.title = 'Signature Score',
           ggtheme = theme_bw(), # Change ggplot2 theme
           tables.theme = theme_void(),
           palette = c("#E31A1C","#1F78B4"), #palette = c("red", "blue")
           legend = c(0.85,0.85),#"right",
           title='Progression Free Survival of MGH On-treatment Samples',
           pval.coord = c(600,0.75),
           pval.method.coord = c(800,0.85),
           xlim = c(0, 1350)
           )
dev.off()

########################################
#### All test  Sample Surv Analysis
########################################。
Gide.res.OS.surv2 <- Gide.res[,c('OS','OS_event','GruopSS')]

Lee.res.OS.surv2 <- Lee.res[,c('OS.days', 'OS_stauts','GruopSS')]
colnames(Lee.res.OS.surv2)[1:2] <- c('OS','OS_event')

MGH.res.OS.surv2 <- MGH.res[,c('os','censure_os','GroupSS')]
colnames(MGH.res.OS.surv2) <- c('OS','OS_event','GruopSS')

OS_all_Test_Samples <- rbind(Gide.res.OS.surv2, Lee.res.OS.surv2, MGH.res.OS.surv2)

all.OS.test.sample.fit<-survfit(Surv(OS, OS_event) ~ GruopSS, data = OS_all_Test_Samples)
tiff(file="output/PASS_ON/All_Test_Sample_PASS_ON_OS_odds_mean.tiff",width=7, height=5,unit='in',compression = 'lzw',res=150)
ggsurvplot(all.OS.test.sample.fit,
           pval = TRUE, 
           pval.method = FALSE,
           test.for.trend = FALSE,
           risk.table = TRUE, # Add risk table
           risk.table.col = 'black',
           risk.table.y.text = TRUE,
           risk.table.y.text.col = TRUE,
           #conf.int = TRUE,
           legend.labs = c("High (n=23)", "Low (n=33)"),
           legend.title = 'Signature Score',
           ggtheme = theme_bw(), # Change ggplot2 theme
           tables.theme = theme_void(),
           palette = c("#E31A1C","#1F78B4"),
           legend = c(0.12,0.20),#"right",
           pval.coord = c(1700,0.8),
           #xlim = c(0, 2400),
           title='Overall Survival of all On-treatment Test Samples'
)
dev.off()


#PFS
Gide.res.pfs.surv2 <- Gide.res[,c('PFS','PFS_event','GruopSS')]
MGH.res.pfs.surv2 <- MGH.res[,c('pfs', 'censure_pfs','GroupSS')]
colnames(MGH.res.pfs.surv2) <- c('PFS','PFS_event','GruopSS')
all.test.sample.pfs <- rbind(Gide.res.pfs.surv2,MGH.res.pfs.surv2)
table(all.test.sample.pfs$GruopSS)

all.test.sample.PFS.fit<-survfit(Surv(PFS, PFS_event) ~ GruopSS, data = all.test.sample.pfs)
tiff(file="output/PASS_ON/All_Test_Sample_PASS_ON_PFS_odds_mean.tiff",width=7, height=5,unit='in',compression = 'lzw',res=150)
ggsurvplot(all.test.sample.PFS.fit,
           pval = TRUE, 
           pval.method = FALSE,
           test.for.trend = FALSE,
           risk.table = TRUE, # Add risk table
           risk.table.col = 'black',
           risk.table.y.text = TRUE,
           risk.table.y.text.col = TRUE,
           #conf.int = TRUE,
           legend.labs = c("High (n=12)", "Low (n=21)"),
           legend.title = 'Signature Score',
           ggtheme = theme_bw(), # Change ggplot2 theme
           tables.theme = theme_void(),
           palette = c("#E31A1C","#1F78B4"),
           legend = c(0.85,0.85),#"right",
           pval.coord = c(900,0.7),
           xlim = c(0, 1350),
           title='Progression Free Survival of all On-treatment Test Samples'
)
dev.off()
