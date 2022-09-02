
# setwd("~/Documents/R/Genbank/bryophytes_2022/")

library(tidyverse)
library(bold)    # API interface to BOLD
library(taxize)  # for NCBI taxonomy lookup
library(seqinr)  # for FASTA output
library(rentrez)
library(readxl)

options(ENTREZ_KEY = "b2ccdcd6619a29d3bf31de74e7cde9a1c209")

# import data / clean ---------------------------------
tax_0<-#read_csv("Data/VocabRequest_ArhtropodSTE.csv") 
  read_xlsx("Data/VocabRequest_ArthropodSTE.xlsx", sheet="All arthropods")
tax<-tax_0  %>% 
  mutate(Species=trimSpace(FinalID)) %>%
  # filter(!is.na(Family))
  filter(TaxonomicLevelCode==40)
taxa<-tax$Species
#taxa<-(droplevels(taxa[1:10]))
sort(taxa)

###############################
# Spellcheck -----------
###########################
# http://viktoriawagner.weebly.com/blog/cleaning-species-names-with-r-ii-taxize

src <- c("EOL", "The International Plant Names Index",
         "Index Fungorum", "ITIS", "Catalogue of Life",
         "Tropicos - Missouri Botanical Garden")

result.long <-gnr_resolve(sci = as.character(tax$Species), data_source_ids = c(1,5,12,165,167), 
                          with_canonical_ranks=T)
# write.table(result.long, file="clipboard", sep="\t", row.names=F)
write_csv(result.long, "Outputs/Arthropod_spellcheck.csv")

result.long %>% select(user_supplied_name, matched_name2) %>% unique() %>%
  group_by(user_supplied_name)  %>% mutate(n=length(matched_name2)) %>% filter(n>1) 

# result.long %>% filter(user_supplied_name == "Tricholepididae")

taxa<-#sort(unique(result.long$matched_name2))
  intersect(taxa, result.long$matched_name2) %>% #Remove those with bad or unknown spellings
  sort() %>%
  unique()

###################################################
# BOLD -------------------------------------------------
###################################################
#Check a few taxa to see what markers are widely used?
# junk<-bold_seqspec("Zygentoma")
# junk$markercode %>% unique() %>% sort()
# junk %>% select(markercode) %>% unique() %>% dput()

#Carabidae markers c("COI-5P", "COI-3P", "", "16S", "18S-3P", "ITS2", "COI-5PNMT1", "ND2", "CYTB", "COII", "ND3", "ND4L", "ND6", "COXIII", "ND1", "ND4", "ND5-0")     
#Formicidae markers c("COI-5P", "", "28S-D2", "CYTB", "LWRHO", "12S", "H3", "16S", "18S", "ITS", "ITS2", "COI-NUMT", "COI-3P", "COI-PSEUDO", "28S-D2-D3", "28S", "H3-NUMT", "CYTB-NUMT", "Wnt1")
#Lycosidae markers c("COI-5P", "18S-3P", "COI-3P", "","16S", "ND6", "ND4L", "ND3", "ND2", "ND1", "COII", "COXIII", "CYTB", "ND5-0", "ND4")
#Pseudoscorpiones markers  c("COI-5P", "", "COI-3P", "ND3",  "COII", "ND4L", "ND2", "ND4", "ND5-0", "ND6", "ND1", "COXIII", "CYTB")
#Dermaptera markers c("COI-5P", "COI-3P", "", "COXIII", "ND4", "COII", "CYTB", "ND1", "ND2", "ND4L", "ND5-0", "ND3", "ND6")
#Archaeognatha markers c("COI-5P", "", "ND4L", "ND3", "ND1", "COII", "ND2", "ND6", "COXIII", "CYTB", "ND4", "ND5-0")
#Zygentoma markers  c("COI-5P", "", "COI-3P", "28S",  "16S", "Wnt1", "ND1", "ND3", "ND6", "COXIII", "ND4L", "COII", "ND4", "CYTB", "ND2", "ND5-0")
bold_df_summary<-  lapply(taxa, function(x){
  print(x)
  mydf=bold_seqspec(x) 
  if(class(mydf)=="data.frame"){ 
    mydf<-mydf %>% 
      mutate(Locus = case_when(#str_detect(seq_primers %>% tolower(),"rbcl")~"rbcL",
                               #str_detect(seq_primers %>% tolower(),"trnl")~"trnL",
                               str_detect(marker_codes %>% tolower(),"its")~"ITS", #Ribosomal
                               str_detect(markercode %>% tolower(),"its")~"ITS", #Ribosomal
                               str_detect(marker_codes %>% tolower(),"coi")~"COI", #Mitochondrial
                               str_detect(marker_codes %>% tolower(),"co1")~"COI", #Mitochondrial
                               str_detect(markercode %>% tolower(),"coi")~"COI",#Mitochondrial
                               str_detect(markercode %>% tolower(),"co1")~"COI",#Mitochondrial
                               
                               #
                               str_detect(markercode %>% tolower(),"18s")~"18S", #Ribosomal
                               str_detect(marker_codes %>% tolower(),"18s")~"18S", #Ribosomal
                               str_detect(markercode %>% tolower(),"16s")~"16S", #Ribosomal but sometimes mitochondrial
                               str_detect(marker_codes %>% tolower(),"16s")~"16S", #Ribosomal but sometimes mitochondrial
                               str_detect(markercode %>% tolower(),"12s")~"12S", #Mitochondrial
                               str_detect(marker_codes %>% tolower(),"12s")~"12S", #Mitochondrial
                               str_detect(markercode %>% tolower(),"cytb")~"CYTB", #Mitochondrial ##SAME AS COI?
                               str_detect(marker_codes %>% tolower(),"cytb")~"CYTB", #Mitochondrial ##SAME AS COI?
                               str_detect(markercode %>% tolower(),"coxiii")~"COI", #Mitochondrial
                               str_detect(marker_codes %>% tolower(),"coxiii")~"COI", #Mitochondrial
                               
                               T~"Other"))%>%
      filter(Locus!="Other")
    xdf<-tibble(FinalID=x,
                BOLD_Records=nrow(mydf),
                BOLD_Records_USA=sum(mydf$country=="United States", na.rm=T),
                BOLD_Records_CA=sum(mydf$province_state=="California", na.rm=T),
                BOLD_Records_SW = sum(mydf$province_state %in% c("California","Arizona", "New Mexico", "Texas", "Nevada", "Utah", "Colorado"), na.rm=T),
                
                #ITS, ribosomal
                BOLD_Records_ITS=sum(mydf$Locus=="ITS", na.rm = T),
                BOLD_Records_ITS_USA=sum(mydf$Locus=="ITS" & mydf$country=="United States", na.rm=T),
                BOLD_Records_ITS_CA=sum(mydf$Locus=="ITS" &mydf$province_state=="California", na.rm=T),
                BOLD_Records_ITS_SW = sum( mydf$Locus=="ITS" &mydf$province_state %in% c("California","Arizona", "New Mexico", "Texas", "Nevada", "Utah", "Colorado"), na.rm=T),
                
                #COI, mitochondrial
                BOLD_Records_COI=sum(mydf$Locus=="COI", na.rm = T),
                BOLD_Records_COI_USA=sum(mydf$Locus=="COI" & mydf$country=="United States", na.rm=T),
                BOLD_Records_COI_CA=sum(mydf$Locus=="COI" &mydf$province_state=="California", na.rm=T),
                BOLD_Records_COI_SW = sum( mydf$Locus=="COI" &mydf$province_state %in% c("California","Arizona", "New Mexico", "Texas", "Nevada", "Utah", "Colorado"), na.rm=T),
                
                #CYTB, mitochondrial
                BOLD_Records_CYTB=sum(mydf$Locus=="CYTB", na.rm = T),
                BOLD_Records_CYTB_USA=sum(mydf$Locus=="CYTB" & mydf$country=="United States", na.rm=T),
                BOLD_Records_CYTB_CA=sum(mydf$Locus=="CYTB" &mydf$province_state=="California", na.rm=T),
                BOLD_Records_CYTB_SW = sum( mydf$Locus=="CYTB" &mydf$province_state %in% c("California","Arizona", "New Mexico", "Texas", "Nevada", "Utah", "Colorado"), na.rm=T),
                
                #12S, ribosomal
                BOLD_Records_12S=sum(mydf$Locus=="12S", na.rm = T),
                BOLD_Records_12S_USA=sum(mydf$Locus=="12S" & mydf$country=="United States", na.rm=T),
                BOLD_Records_12S_CA=sum(mydf$Locus=="12S" &mydf$province_state=="California", na.rm=T),
                BOLD_Records_12S_SW = sum( mydf$Locus=="12S" &mydf$province_state %in% c("California","Arizona", "New Mexico", "Texas", "Nevada", "Utah", "Colorado"), na.rm=T),
                
                #16S, ribosomal
                BOLD_Records_16S=sum(mydf$Locus=="16S", na.rm = T),
                BOLD_Records_16S_USA=sum(mydf$Locus=="16S" & mydf$country=="United States", na.rm=T),
                BOLD_Records_16S_CA=sum(mydf$Locus=="16S" &mydf$province_state=="California", na.rm=T),
                BOLD_Records_16S_SW = sum( mydf$Locus=="16S" &mydf$province_state %in% c("California","Arizona", "New Mexico", "Texas", "Nevada", "Utah", "Colorado"), na.rm=T),
                
                #18S, ribosomal
                BOLD_Records_18S=sum(mydf$Locus=="18S", na.rm = T),
                BOLD_Records_18S_USA=sum(mydf$Locus=="18S" & mydf$country=="United States", na.rm=T),
                BOLD_Records_18S_CA=sum(mydf$Locus=="18S" &mydf$province_state=="California", na.rm=T),
                BOLD_Records_18S_SW = sum( mydf$Locus=="18S" &mydf$province_state %in% c("California","Arizona", "New Mexico", "Texas", "Nevada", "Utah", "Colorado"), na.rm=T)
                
    )
  }
  else  {
    xdf= tibble(FinalID=x,  
                BOLD_Records=0, BOLD_Records_USA=0, BOLD_Records_CA=0, BOLD_Records_SW = 0,
                BOLD_Records_ITS=0, BOLD_Records_ITS_USA=0, BOLD_Records_ITS_CA=0, BOLD_Records_ITS_SW = 0,
                BOLD_Records_COI=0, BOLD_Records_COI_USA=0, BOLD_Records_COI_CA=0, BOLD_Records_COI_SW = 0,
                BOLD_Records_CYTB=0, BOLD_Records_CYTB_USA=0, BOLD_Records_CYTB_CA=0, BOLD_Records_CYTB_SW = 0,
                BOLD_Records_12S=0, BOLD_Records_12S_USA=0, BOLD_Records_12S_CA=0, BOLD_Records_12S_SW = 0,
                BOLD_Records_16S=0, BOLD_Records_16S_USA=0, BOLD_Records_16S_CA=0, BOLD_Records_16S_SW = 0,
                BOLD_Records_18S=0, BOLD_Records_18S_USA=0, BOLD_Records_18S_CA=0, BOLD_Records_18S_SW = 0
                )
  }
  xdf
})  %>% bind_rows()

bold_summary<-bold_df_summary %>%
  rename(QryTax=FinalID) %>%
  inner_join(result.long %>%
               select(QryTax= matched_name2, FinalID=user_supplied_name) %>%
               unique()) %>%
  select(-QryTax) %>%
  pivot_longer(cols=starts_with("BOLD")) %>%
  group_by(FinalID, name) %>%
  summarise(value=sum(value)) %>% 
  ungroup() %>%
  pivot_wider(names_from=name, values_from=value) %>%
  inner_join(tax)


write_csv(bold_summary, file="Outputs/arthro_bold_summary.csv")

###################################################
# GENBANK -------------------------------------------------------------
###################################################

#The function can't handle the whole data set, break into smaller pieces


#mito_search
gb_df_summary<- tibble(taxa = taxa)
gb_df_summary$GB_mito_count<-sapply(gb_df_summary$taxa, function(x){
  term.x = paste0(x,"[PORGN]", " AND 200:3400[SLEN]", " AND mitochondrion[filter]")
  xdf = entrez_search(db="nuccore", term=term.x, FILT=1, api_key="b2ccdcd6619a29d3bf31de74e7cde9a1c209"  )
  print(paste(x, xdf$count))
  xdf$count
})
#mito, full length search
gb_df_summary$GB_mito_fl_count<-sapply(gb_df_summary$taxa, function(x){
  term.x = paste0(x,"[PORGN]", " AND 200:34000000[SLEN]", " AND mitochondrion[filter]")
  xdf = entrez_search(db="nuccore", term=term.x, FILT=1, api_key="b2ccdcd6619a29d3bf31de74e7cde9a1c209"  )
  print(paste(x, xdf$count))
  xdf$count
})
# #chloroplast search
# gb_df_summary$GB_chloro_count<-sapply(gb_df_summary$taxa, function(x){
#   term.x = paste0(x,"[PORGN]", " AND 200:3400[SLEN]", " AND chloroplast[filter]")
#   xdf = entrez_search(db="nuccore", term=term.x, FILT=1, api_key="b2ccdcd6619a29d3bf31de74e7cde9a1c209"  )
#   print(paste(x, xdf$count))
#   xdf$count
# })

#18S search
gb_df_summary$GB_18S_count<-sapply(gb_df_summary$taxa, function(x){
  term.x = paste0(x,"[PORGN]", " AND 200:3400[SLEN]", " AND 18S[All fields]")
  xdf = entrez_search(db="nuccore", term=term.x, FILT=1, api_key="b2ccdcd6619a29d3bf31de74e7cde9a1c209"  )
  print(paste(x, xdf$count))
  xdf$count
})




gb_summary<-gb_df_summary %>%
  rename(QryTax=taxa) %>%
  inner_join(result.long %>%
               select(QryTax= matched_name2, FinalID=user_supplied_name) %>%
               unique()) %>%
  select(-QryTax) %>%
  pivot_longer(cols=starts_with("GB_")) %>%
  group_by(FinalID, name) %>%
  summarise(value=sum(value)) %>% 
  ungroup() %>%
  pivot_wider(names_from=name, values_from=value) %>%
  inner_join(tax)


write_csv(gb_summary, file="Outputs/arthro_gb_summary.csv")

tax_bold_genbank<-tax %>%
  left_join(bold_df_summary %>%
              rename(QryTax=FinalID) %>%
              inner_join(result.long %>%
                           select(QryTax= matched_name2, FinalID=user_supplied_name) %>%
                           unique()) %>%
              select(-QryTax) %>%
              pivot_longer(cols=starts_with("BOLD")) %>%
              group_by(FinalID, name) %>%
              summarise(value=sum(value)) %>% 
              ungroup() %>%
              pivot_wider(names_from=name, values_from=value) ) %>%
  left_join(gb_df_summary %>%
              rename(QryTax=taxa) %>%
              inner_join(result.long %>%
                           select(QryTax= matched_name2, FinalID=user_supplied_name) %>%
                           unique()) %>%
              select(-QryTax) %>%
              pivot_longer(cols=starts_with("GB_")) %>%
              group_by(FinalID, name) %>%
              summarise(value=sum(value)) %>% 
              ungroup() %>%
              pivot_wider(names_from=name, values_from=value) )

write_csv(tax_bold_genbank, file="Outputs/arthros_bold_genbank.csv")



