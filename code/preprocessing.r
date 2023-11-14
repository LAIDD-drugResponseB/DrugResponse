library(dplyr)
library(data.table)
library(tidyr)

getwd()

## gene.exp
dirv = '' # input data가 있는 경로 넣어주기

RNA<-fread(paste0(driv,'/OmicsExpressionProteinCodingGenesTPMLogp1.csv'))
RNA %>% dim # 1450 19194

#### for RNA expression data annotation (ARXSPAN_ID(V1) -> COSMIC.ID)
drug.res = read.csv(paste0(dirv,'/sanger-dose-response.csv'))

anno<-drug.res[,c("COSMIC_ID","ARXSPAN_ID")]

RNA_name<-data.frame(RNA$V1)
colnames(RNA_name)<-"ARXSPAN_ID"
name_merge<-merge(anno,RNA_name,by="ARXSPAN_ID")

name_merge<- name_merge %>% unique
name_merge$COSMIC_ID %>% unique %>% length #696

colnames(RNA)<-gsub('V1',"ARXSPAN_ID",colnames(RNA))

which(colnames(RNA)%in%"ARXSPAN_ID") #1개만 바뀌었는지 확인
merge_RNA<-merge(name_merge,RNA,by="ARXSPAN_ID")
merge_RNA[1:10,1:3]
identical(merge_RNA$COSMIC_ID,unique(merge_RNA$COSMIC_ID))
merge_RNA<-merge_RNA[,-c(1)]
merge_RNA %>% colnames %>% head

### gene name (*) remove 
col_exp<-colnames(merge_RNA)
col_exp %>% length # 19194
col_exp<-gsub("\\s*\\(\\d+\\)","",col_exp) 
# \\s 공백 
# \\( ( 
# \\d 숫자 
# \\) ) 
# + 1개이상
# * 0개이상
# gsub(a,b,data) => a를 b로 변경해주는 코드 
# 위의 상황에서는 " (숫자-여러개)" -> ""로 변경 (없애줌)
col_exp %>% head
identical(merge_RNA %>% colnames() %>% length,length(col_exp))
colnames(merge_RNA)<-col_exp
gene.exp %>% dim # 664 19194
merge_RNA %>% dim # 696 19194

gene.exp %>% colnames() %>% head
gene.exp$"COSMIC.ID" %>% head

write.csv(merge_RNA,file=paste0(dirv,"/RNA_express_preprocessing.csv"))
merge_RNA %>% dim #696 19194


###### gene.mut 
# 하나의 컬럼을 여러개로 IS Mutated -> gene 이름으로

genetic<-fread(paste0(dirv,'PANCANCER_Genetic_features_variant_GDSC1.csv')
genetic %>% head
genetic %>% colnames 
#"Cell Line Name" "COSMIC ID" "GDSC Desc1" "GDSC Desc2" "TCGA Desc" 
#"Genetic Feature" "IS Mutated" "Recurrent Gain Loss" "Genes in Segment"  

genetic %>% dim # 289540 9

genetic_2<-fread(paste0(dirv,'/PANCANCER_Genetic_features_variant_GDSC2.csv'))
genetic_2 %>% head
genetic_2 %>% dim #288300 9
'''
[1] "Cell Line Name"      "COSMIC ID"           "GDSC Desc1"         
[4] "GDSC Desc2"          "TCGA Desc"           "Genetic Feature"    
[7] "IS Mutated"          "Recurrent Gain Loss" "Genes in Segment"  

'''
genetic<-rbind(genetic,genetic_2)
genetic$"COSMIC ID" %>% unique %>% length # 939 cosmicID

genetic<-genetic %>% unique 

genetic$'Genetic Feature' %>% head
genetic$'Genetic Feature'<-gsub("(*)_mut", "",genetic$'Genetic Feature')
genetic$'COSMIC ID' %>% unique %>% length #939
genetic %>% head
genetic %>% colnames

merge_RNA %>% colnames %>% head
genetic %>% colnames

genetic<-genetic[,c("COSMIC ID","Genetic Feature","IS Mutated")] %>% data.frame
colnames(genetic)<-c("COSMIC_ID","Genetic_Feature","IS_Mutated")
merged_omics<-merge(genetic,merge_RNA,by="COSMIC_ID")
RNA_name<-merge_RNA [,1] %>% data.frame
colnames(RNA_name)<-"COSMIC_ID" 
merged_omics<-merge(genetic,RNA_name,by="COSMIC_ID")
merged_omics %>% dim
merged_omics<-merged_omics %>% unique

merged_omics$COSMIC_ID %>% unique %>% length # 661 
merged_omics<-merged_omics %>% unique

temp_name<-merged_omics$COSMIC_ID %>% unique
setdiff(gene.mut$COSMIC.ID %>% unique,temp_name)


table(genetic$'COSMIC ID',genetic$'Genetic Feature') %>% dim
melted_data<-genetic[,c(2,6,7)]
melted_data %>% head
melted_data<-melted_data %>% data.frame
colnames(melted_data)<-c("COSMIC.ID","mut","gene")

spread_data <- spread(melted_data, key = "mut", value = "gene") #spread <-> gather
write.csv(spread_data,file=paste0(dirv,"/genetic_data.csv"))

#####gene expression & gene mutation cosmic ID merge
mut_raw<-read.csv(file=paste0(dirv,"/genetic_data.csv"))
mut_raw %>% colnames %>% head
mut_raw$COSMIC.ID %>% head
intersect(mut_raw$COSMIC.ID,merge_RNA$COSMIC_ID) %>% length 
inter_cosmic<-intersect(mut_raw$COSMIC.ID,merge_RNA$COSMIC_ID)
which(mut_raw$COSMIC.ID%in%inter_cosmic) %>% length

merge_RNA$COSMIC.ID %>% head
merge_RNA %>% colnames %>% head

## cell.info (gene.mut 먼저 preprocessing 하기) 
#

cell.info = read.csv(file=paste0(dirv,'/Cell_info_filtered.csv'))
cell.info %>% dim # 664 5
cell.info  %>% colnames

# "Cell.Line.Name" "COSMIC.ID" "GDSC.Desc1" "GDSC.Desc2" "TCGA.Desc"     

genetic %>% colnames
'''[1] "Cell Line Name"      "COSMIC ID"           "GDSC Desc1"         
[4] "GDSC Desc2"          "TCGA Desc"           "Genetic Feature"    
[7] "IS Mutated"          "Recurrent Gain Loss" "Genes in Segment"   '''

cell_info<-genetic[,1:5]
colnames(cell_info)<-c("Cell.Line.Name", "COSMIC.ID", "GDSC.Desc1", "GDSC.Desc2", "TCGA.Desc")

#drug.res  # 변경할 사항 없음

drug.res = read.csv(paste0(dirv,'/sanger-dose-response.csv'))

drug.res %>% head


## drug.info # 변경할 사항 없음

drug.info = read.csv(file=paste0(dirv,'/screened_compounds_rel_8.4.csv'))
