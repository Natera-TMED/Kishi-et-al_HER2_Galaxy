# https://docs.google.com/document/d/1xWoVdUrsUKNvnCyHdQppwWT-0PzcxS9R/edit

rm(list=ls())
library(data.table)
library(survival)
library(survminer)
library(ggplot2)
library(reshape2)
library(psych)
library(flextable)
library(gtsummary)
library(vautils)
library(openxlsx)
library(cards)
library(cardx)
library(broom)
library(officer)
library(tidyverse)

source('~/gitlab/translationaloncology/tmed_stat.R')

dataset=fread('Demographics for Natera.csv')
dim(dataset)


names(dataset)
dim(dataset)
dataset[,DFS.event:=`DFS event=1 free=0`]
# HER2 amp (+) was observed in 59 cases (1.6%). 
table(dataset$`HER2 amplification`)
prop.table(table(dataset$`HER2 amplification`))
prop.table(table(dataset[,.(`HER2 amplification`,`Stage I,II,III`)]),margin=2)
dataset$`HER2 amplification`=as.factor(dataset$`HER2 amplification`)
dataset[,histological_type:=`Degree of Histological Differentiation por/sig or not`]
dataset[histological_type=='',histological_type:=NA]
dataset$histological_type=factor(dataset$histological_type,level=c('Differentiated','Por/Sig'))
table(dataset$histological_type)
dataset[,location:=`Location 解析4 colon vs rectum`]
table(dataset$location)
dataset[,left_right:=`Location classification right vs left vs rectum`]
dataset[`Location classification right vs left vs rectum`=='',`Location classification right vs left vs rectum`:=NA]
dataset$`Location classification right vs left vs rectum`=factor(dataset$`Location classification right vs left vs rectum`,level=c('RIGHT','LEFT','RECTUM'))
table(dataset$`Location classification right vs left vs rectum`)
table(dataset[,.(left_right,`HER2 amplification`)])
dataset[left_right=='',left_right:=NA]
dataset$left_right=factor(dataset$left_right,level=c('RIGHT','LEFT'))
# dataset[left_right=='RECTUM',left_right:=NA]
dataset[`MSI category`=='',`MSI category`:=NA]
dataset[RAS=='',RAS:=NA]
dataset[RAS=='Not Done',RAS:=NA]
dataset[RAS=='Unknown',RAS:=NA]
dataset[`BRAF category`=='',`BRAF category`:=NA]
dataset[`ctDNA Result 4W`=='TNP',`ctDNA Result 4W`:=NA]
dataset[`ctDNA Result 4W`=='NC',`ctDNA Result 4W`:=NA]
dataset[`ctDNA Result 4W`=='',`ctDNA Result 4W`:=NA]
dataset[,StageIII:=fifelse(`Stage I,II,III`=='III',TRUE,fifelse(`Stage I,II,III` %in% c('I','II'),FALSE,NA))]
dataset[,MSI:=factor(`MSI category`,levels=c('MSI-HIGH','MSS'))]
dataset[,BRAF:=`BRAF category`]
dataset[,ACT:=`Perioperative treatmentChemotherapy`]
dataset[,`ctDNA status at 4 weeks post surgery`:=`ctDNA Result 4W`]
dataset[,`Histological type`:=histological_type]
dataset[`Histological type`=='Por/Sig',`Histological type`:='Poorly differentiated /Signet ring']
dataset[,HER2_amplification:=`HER2 amplification`]


# dataset[!is.na(DFS)] %>% dim() #3152  178
# dataset[!is.na(`DFS event=1 free=0`)] %>% dim() #3153  178
# dataset[!is.na(`Recurrence free period`)] %>% dim() #3082  178
# dataset[!is.na(`rec=1 rec fee=0`)] %>% dim() #3155  178
# table(dataset[,.(`rec=1 rec fee=0`,`DFS event=1 free=0`)],useNA='always')
# dataset[`rec=1 rec fee=0`==0 & `DFS event=1 free=0`==1,.(`Date of death`,`Cause of death`,Outcome)]
# 
# dataset
# dataset[,RFS:=`Recurrence free period`]
# dataset[,RFS.event:=`rec=1 rec fee=0`]
# View(dataset[(`rec=1 rec fee=0`!=`DFS event=1 free=0`),])
# View(dataset[(`rec=1 rec fee=0`==1 & Outcome=='Dead'),.()])
# dataset[,.(`Date of recurrence`,`Cause of death`,`Date of death`)] %>% head(100)
# dataset[,date_recur:=as.Date(`Date of recurrence`,format='%m/%d/%y')]
# dataset[,cause_of_death:=as.Date(`Cause of death`,format='%m/%d/%y')]
# dataset[,date_death:=as.Date(`Date of death`,format='%m/%d/%y')]
# dataset[,.(date_recur,cause_of_death,date_death)]
# dataset[date_recur<date_death,.(RFS.event,`DFS event=1 free=0`)]
# dataset[is.na(date_recur) & `DFS event=1 free=0`==1,.(RFS.event,date_death,cause_of_death)]
# 
# dataset[`DFS event=1 free=0`==1 & RFS.event==0,.(date_recur,cause_of_death,date_death)]


dt.pat.table.source = dataset[, .(`HER2 amplification`, Sex, Age, `Location classification`=`Location classification right vs left vs rectum`,histological_type,`Stage I,II,III` ,RAS,`BRAF V600E`,`MSI category`,`ctDNA Result 4W`)]
dt.pat.table.source[!`BRAF V600E` %in% c('Negative (wild type)','Positive (mutant)'),`BRAF V600E`:=NA]

reset_gtsummary_theme()
theme_gtsummary_compact()


table <- dt.pat.table.source |>
  tbl_summary(
    by = `HER2 amplification`,
    missing = "no"
  ) |>
  
  add_p(
    test = list(all_categorical() ~ "fisher.test"),
    pvalue_fun = ~style_pvalue(.x, digits = 3)
  ) |>
  bold_p(t = 0.05)
table

table_flextable=as_flex_table(table)
doc <- read_docx() %>%
  body_add_flextable(table_flextable)
print(doc,target='Table1.docx')


###########################
## check consort diagram ##
###########################
dim(dataset) #3607

dataset[!is.na(DFS)] %>% dim() #3152  178
dataset[!is.na(`DFS event=1 free=0`)] %>% dim() #3153  178
dataset[!is.na(RFS)] %>% dim() #3082  178

# follow up
dataset[!is.na(DFS),DFS] %>% summary()
dataset[!is.na(RFS),DFS] %>% summary()

# check DFS/RFS subset
dataset[!is.na(RFS) & is.na(DFS)] %>% dim() ## no patient with RFS that has no DFS

##associaiton of HER
dataset[!is.na(RFS) & !is.na(`rec=1 rec fee=0`) & !is.na(HER2_amplification)] %>% dim() #3082  178
dataset[!is.na(RFS) & !is.na(`rec=1 rec fee=0`),HER2_amplification]  %>% table()
dataset[is.na(RFS) | is.na(`rec=1 rec fee=0`)] %>% dim() #525
RFSdataset=dataset[!is.na(RFS) & !is.na(`rec=1 rec fee=0`)] 
dim(RFSdataset)

## effect of ACT
dataset[!is.na(DFS) & (ACT) %in% c('Yes','No') & StageIII==TRUE] %>% dim() #1520  178
dataset[!is.na(RFS) & (ACT) %in% c('Yes','No') & StageIII==TRUE] %>% dim() #1520  178
dataset[!is.na(DFS) & (ACT) %in% c('Yes','No') & `Stage classification [pStage] (UICC TNM Classification, 8th Edition)` %like% 'III' ] %>% dim() #1520  178
RFSdataset[!StageIII==TRUE]%>% dim()
RFSdataset[(!StageIII==TRUE) | !(ACT %in% c('Yes','No')) ]%>% dim() #1595
dataset[!is.na(RFS) & !is.na(`rec=1 rec fee=0`) & (ACT) %in% c('Yes','No') & StageIII==TRUE] %>% dim() #1487


## copy number with recurrence
dataset[!is.na(RFS) & !is.na(`rec=1 rec fee=0`) & HER2_amplification=='negative'] %>% dim()
dataset[!is.na(DFS) & HER2_amplification=='positive'] %>% dim() #48
dataset[!is.na(RFS) & !is.na(`rec=1 rec fee=0`) & HER2_amplification=='positive'] %>% dim() #46
HERpositiveRFSdataset=dataset[!is.na(RFS) & !is.na(`rec=1 rec fee=0`) & HER2_amplification=='positive']

##pre-operate
dataset[!is.na(DFS) & HER2_amplification=='positive' & `ctDNA Result initial`=='Positive'] %>% dim() # 41
HERpositiveRFSdataset[!is.na(RFS) & !is.na(`rec=1 rec fee=0`) & HER2_amplification=='positive' & `ctDNA Result initial`=='Positive'] %>% dim() # 39
HERpositiveRFSdataset[ `ctDNA Result initial`!='Positive'] %>% dim() # 39

## copy number and ctDNA MRD
dataset[!is.na(DFS) & HER2_amplification=='positive' & !is.na(`ctDNA Result 4W`) ] %>% dim() # 43
HERpositiveRFSdataset[!is.na(`ctDNA Result 4W`) ] %>% dim() # 41
HERpositiveRFSdataset[is.na(`ctDNA Result 4W`) ] %>% dim() # 41




######################
## HER2 copy cutoff ##
######################
pdf('FigureS1.pdf')
ggplot(dataset, aes(x = `copy number`)) +
  geom_histogram(binwidth = 1, boundary = 0.5, color = "white") +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(1, 50, 2)) 
dev.off()