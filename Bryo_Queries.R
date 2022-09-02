
# setwd("~/Documents/R/Genbank/bryophytes_2022/")

library(tidyverse)
library(bold)    # API interface to BOLD
library(taxize)  # for NCBI taxonomy lookup
library(seqinr)  # for FASTA output
library(rentrez)
library(readxl)

options(ENTREZ_KEY = "b2ccdcd6619a29d3bf31de74e7cde9a1c209")

# import data / clean ---------------------------------
tax_0<-read_csv("Data/VocabRequest_BryophyteSTE.csv") 
tax<-tax_0  %>% 
  mutate(Species=trimSpace(FinalID)) %>%
  filter(!is.na(Family))
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

write_csv(result.long, "Outputs/Bryophyte_spellcheck.csv")

taxa<-sort(unique(result.long$matched_name2))


#Run this code to find out what markers are commonly used in BOLD
# junk<-bold_seqspec("Funariaceae")
# junk %>% select(markercode) %>% unique() %>% dput()
#Bryum c("rbcL", "trnL-F", "ITS2", "ITS", "", "rbcLa", "trnH-psbA", "COI-5P")
#Tortula c("rbcLa", "trnL-F", "ITS2", "rbcL", "ITS", "")
#Timmiella  c("trnL-F", "ITS2", "rbcL")
#Pohlia c("ITS2", "trnL-F", "rbcLa", "", "rbcL", "ITS")
#Pottiaceae c("ITS2", "trnL-F", "rbcLa", "rbcL", "matK", "ITS1", "ITS", "", "atpF-atpH", "trnH-psbA")
###################################################
# BOLD -------------------------------------------------
###################################################

#Look at seq_primers and marker_codes to see if it's rbcL vs trnL vs ITS COI or CO1 or matK
#It's usually better to look at "markercode", but occasionally "marker_codes" has stuff too when markercode is blank
bold_df_summary<-  lapply(taxa, function(x){
  print(x)
  mydf=bold_seqspec(x) 
  if(class(mydf)=="data.frame"){ 
    mydf<-mydf %>% 
      mutate(Locus = case_when( str_detect(markercode %>% tolower(),"rbcl")~"rbcL", #chloroplast
                                str_detect(markercode %>% tolower(),"trnl")~"trnL", #chloroplast
                                str_detect(markercode %>% tolower(),"its")~"ITS", #nuclear ribosomal
                                str_detect(markercode %>% tolower(),"coi")~"COI",#mitochondria
                                str_detect(markercode %>% tolower(),"co1")~"COI", #mitochondria
                                str_detect(marker_codes %>% tolower(),"rbcl")~"rbcL", 
                                str_detect(marker_codes %>% tolower(),"trnl")~"trnL",
                                str_detect(marker_codes %>% tolower(),"its")~"ITS",
                                str_detect(marker_codes %>% tolower(),"coi")~"COI",
                                str_detect(marker_codes %>% tolower(),"co1")~"COI", 
                                
                                
                                
                                str_detect(markercode %>% tolower(),"matk")~"matK", #Chloroplast
                                str_detect(marker_codes %>% tolower(),"matk")~"matK", #Chloroplast
                                
                                T~"Other")) %>%
      filter(Locus!="Other")
    xdf<-tibble(FinalID=x,
                #Any sequence
                BOLD_Records=nrow(mydf),
                BOLD_Records_USA=sum(mydf$country=="United States", na.rm=T),
                BOLD_Records_CA=sum(mydf$province_state=="California", na.rm=T),
                BOLD_Records_SW = sum(mydf$province_state %in% c("California","Arizona", "New Mexico", "Texas", "Nevada", "Utah", "Colorado"), na.rm=T),
                
                #rbcL sequences. rbcL is a chloroplast marker
                BOLD_Records_rbcL=sum(mydf$Locus=="rbcL", na.rm = T),
                BOLD_Records_rbcL_USA=sum(mydf$Locus=="rbcL" & mydf$country=="United States", na.rm=T),
                BOLD_Records_rbcL_CA=sum(mydf$Locus=="rbcL" &mydf$province_state=="California", na.rm=T),
                BOLD_Records_rbcL_SW = sum(mydf$Locus=="rbcL" &mydf$province_state %in% c("California","Arizona", "New Mexico", "Texas", "Nevada", "Utah", "Colorado"), na.rm=T),
                
                #trnL sequences. trnL is a chloroplast marker
                BOLD_Records_trnL=sum(mydf$Locus=="trnL", na.rm = T),
                BOLD_Records_trnL_USA=sum(mydf$Locus=="trnL" & mydf$country=="United States", na.rm=T),
                BOLD_Records_trnL_CA=sum(mydf$Locus=="trnL" &mydf$province_state=="California", na.rm=T),
                BOLD_Records_trnL_SW = sum( mydf$Locus=="trnL" &mydf$province_state %in% c("California","Arizona", "New Mexico", "Texas", "Nevada", "Utah", "Colorado"), na.rm=T),
                
                #ITS sequences. ITS is a ribosomal nuclear marker
                BOLD_Records_ITS=sum(mydf$Locus=="ITS", na.rm = T),
                BOLD_Records_ITS_USA=sum(mydf$Locus=="ITS" & mydf$country=="United States", na.rm=T),
                BOLD_Records_ITS_CA=sum(mydf$Locus=="ITS" &mydf$province_state=="California", na.rm=T),
                BOLD_Records_ITS_SW = sum( mydf$Locus=="ITS" &mydf$province_state %in% c("California","Arizona", "New Mexico", "Texas", "Nevada", "Utah", "Colorado"), na.rm=T),
                
                #COI sequences. COI is a mitochondrial marker
                BOLD_Records_COI=sum(mydf$Locus=="COI", na.rm = T),
                BOLD_Records_COI_USA=sum(mydf$Locus=="COI" & mydf$country=="United States", na.rm=T),
                BOLD_Records_COI_CA=sum(mydf$Locus=="COI" &mydf$province_state=="California", na.rm=T),
                BOLD_Records_COI_SW = sum( mydf$Locus=="COI" &mydf$province_state %in% c("California","Arizona", "New Mexico", "Texas", "Nevada", "Utah", "Colorado"), na.rm=T),
                
                #matK, a chloroplast marker
                BOLD_Records_matK=sum(mydf$Locus=="matK", na.rm = T),
                BOLD_Records_matK_USA=sum(mydf$Locus=="matK" & mydf$country=="United States", na.rm=T),
                BOLD_Records_matK_CA=sum(mydf$Locus=="matK" &mydf$province_state=="California", na.rm=T),
                BOLD_Records_matK_SW = sum( mydf$Locus=="matK" &mydf$province_state %in% c("California","Arizona", "New Mexico", "Texas", "Nevada", "Utah", "Colorado"), na.rm=T)
                
    )
  }
  else  {
    xdf= tibble(FinalID=x,  
                BOLD_Records=0, BOLD_Records_USA=0, BOLD_Records_CA=0, BOLD_Records_SW = 0,
                BOLD_Records_rbcL=0, BOLD_Records_rbcL_USA=0, BOLD_Records_rbcL_CA=0, BOLD_Records_rbcL_SW = 0,
                BOLD_Records_trnL=0, BOLD_Records_trnL_USA=0, BOLD_Records_trnL_CA=0, BOLD_Records_trnL_SW = 0,
                BOLD_Records_matK=0, BOLD_Records_matK_USA=0, BOLD_Records_matK_CA=0, BOLD_Records_matK_SW = 0,
                BOLD_Records_COI=0, BOLD_Records_COI_USA=0, BOLD_Records_COI_CA=0, BOLD_Records_COI_SW = 0
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


write_csv(bold_summary, file="Outputs/bryos_bold_summary.csv")

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
#chloroplast search
gb_df_summary$GB_chloro_count<-sapply(gb_df_summary$taxa, function(x){
  term.x = paste0(x,"[PORGN]", " AND 200:3400[SLEN]", " AND chloroplast[filter]")
  xdf = entrez_search(db="nuccore", term=term.x, FILT=1, api_key="b2ccdcd6619a29d3bf31de74e7cde9a1c209"  )
  print(paste(x, xdf$count))
  xdf$count
})

# junk<-paste0("Bryum[PORGN] AND 200:3400[SLEN] AND 18S[All fields]")
# junkx<-entrez_search(db="nuccore", term=junk, FILT=1, api_key="b2ccdcd6619a29d3bf31de74e7cde9a1c209"  )
# junkx$count
# print(junkx)

#18S search. 18S is ribosomal
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


write_csv(gb_summary, file="Outputs/bryos_gb_summary.csv")
bryos_gb_summary<-read_csv("Outputs/bryos_gb_summary.csv")

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

write_csv(tax_bold_genbank, file="Outputs/bryos_bold_genbank.csv")


##############
#Summary
bryos_bold_genbank<-read_csv("Outputs/bryos_bold_genbank.csv")

bryos_bold_genbank %>%
  group_by(TaxonomicLevelCode) %>%
  tally()


bryos_bold_genbank %>%
  filter(TaxonomicLevelCode==40)%>%
  select(FinalID)


bry_boldplot_dat<-bryos_bold_genbank %>%
  select(FinalID, TaxonomicLevelCode, starts_with("BOLD_"), starts_with("GB_")) %>%
  pivot_longer(cols=c(starts_with("BOLD_"), starts_with("GB_")), names_to="Locus", values_to="NSeq", values_drop_na = F) %>%
  mutate(Geography=case_when(str_detect(Locus,"_CA")~"CA",
                             str_detect(Locus,"_SW")~"SW",
                             str_detect(Locus,"_USA")~"USA",
                             T~"World"),
         Loc=case_when(str_detect(Locus,"COI")~"COI",
                       str_detect(Locus,"ITS")~"ITS",
                       str_detect(Locus,"trnL")~"trnL",
                       str_detect(Locus,"rbcL")~"rbcL",
                       str_detect(Locus,"matK")~"matK",
                       T~"Any"),
         DB=case_when(str_detect(Locus,"BOLD_")~"BOLD",
                      str_detect(Locus,"GB_")~"GenBank",
                      T~"XXXX"),
         TaxLevel=case_when(TaxonomicLevelCode==40~"Family",
                            TaxonomicLevelCode==50~"Genus",
                            TaxonomicLevelCode==60~"Species",T~"X")
  ) %>%
  group_by(TaxLevel,TaxonomicLevelCode, Loc, Geography, DB) %>%
  summarise(SeqsPresent=sum(NSeq>1, na.rm=T)) %>%
  ungroup() %>%
  filter(DB=="BOLD", TaxonomicLevelCode<66) %>%
  pivot_wider(names_from=Geography, values_from = SeqsPresent) %>%
  transmute(TaxLevel,Loc,
            CA,
            SW=SW-CA,
            USA=USA-(SW+CA),
            World=World-(CA+SW+USA)) %>%
  pivot_longer(cols=c("CA","SW","USA","World")) 

tot_taxa<-bryos_bold_genbank %>%
  filter(TaxonomicLevelCode<66) %>%
  mutate(TaxLevel=case_when(TaxonomicLevelCode==40~"Family",
                            TaxonomicLevelCode==50~"Genus",
                            TaxonomicLevelCode==60~"Species",T~"X")) %>%
  group_by(TaxLevel) %>% tally()

bry_boldplot_dat2<-bry_boldplot_dat %>%
  bind_rows(
    bry_boldplot_dat %>%
      group_by(TaxLevel, Loc ) %>%
      summarise(value=sum(value)) %>%
      inner_join(tot_taxa) %>%
      mutate(Absent=n-value) %>%
      transmute(TaxLevel, Loc, name="No sequences", value=Absent)
  )
bry_boldplot_dat2$name<-factor(bry_boldplot_dat2$name, levels=rev(c("CA","SW","USA","World","No sequences")))

bry_boldplot_dat2<-bry_boldplot_dat2 %>% arrange(name)

bry_boldplot<-ggplot(data=bry_boldplot_dat2)+ 
  geom_col(aes(x=Loc,
               y=value, 
               fill=name), position=position_stack(), color="black")+
  facet_wrap(~TaxLevel, scales="free_y")+
  scale_y_continuous()+
  scale_fill_brewer(palette="YlOrRd", name="Sequences in BOLD", labels=rev(c("California","Rest of SW", "Rest of USA", "Rest of World", "No sequences")))+
  ylab("# taxa")+xlab("Locus")+
  theme_bw()
bry_boldplot
ggsave(bry_boldplot, filename="Outputs/bry_boldplot.jpg", dpi=300, height=5, width=6.5)
#####

bry_gbplot_dat<-
  bryos_bold_genbank %>%
  select(FinalID, TaxonomicLevelCode,  starts_with("GB_")) %>%
  pivot_longer(cols= c(starts_with("GB_")), names_to="Locus", values_to="NSeq", values_drop_na = F)  %>%
  mutate(Loc=case_when(str_detect(Locus,"mito_fl")~"Mito\n(full)",
                       str_detect(Locus,"mito_")~"Mito\n(short)",
                       str_detect(Locus,"chloro")~"Chl",
                       str_detect(Locus,"18S")~"18S",
                       T~"Any"),
         DB=case_when(str_detect(Locus,"BOLD_")~"BOLD",
                      str_detect(Locus,"GB_")~"GenBank",
                      T~"XXXX"),
         TaxLevel=case_when(TaxonomicLevelCode==40~"Family",
                            TaxonomicLevelCode==50~"Genus",
                            TaxonomicLevelCode==60~"Species",T~"X")
  ) %>%
  group_by(TaxLevel,TaxonomicLevelCode, Loc, DB) %>%
  summarise(SeqsPresent=sum(NSeq>1, na.rm=T)) %>%
  ungroup() %>%
  filter(DB=="GenBank", TaxonomicLevelCode<66) 


bry_gbplot_dat2<-bry_gbplot_dat %>%
  group_by(TaxLevel, Loc ) %>%
  summarise(Present=sum(SeqsPresent)) %>%
  inner_join(tot_taxa) %>%
  mutate(Absent=n-Present) %>%
  select(-n) %>%
  pivot_longer(cols=c(Present, Absent))

bry_gbplot<-ggplot(data=bry_gbplot_dat2)    +
  geom_col(aes(x=Loc,
               y=value, 
               fill=name), position=position_stack(), color="black")+
  facet_wrap(~TaxLevel, scales="free_y")+
  scale_y_continuous()+
  scale_fill_manual(values=c("#ffffb2","#bd0026"), name="Sequences in GenBank" , labels=(c("No sequences", "Present")))+
  ylab("# taxa")+xlab("Locus")+
  theme_bw()

library(ggpubr)
bold_genbank_plots<-ggarrange(bry_boldplot, bry_gbplot, ncol=1, align="hv")
ggsave(bold_genbank_plots, filename="Outputs/bold_genbank_plots.jpg", dpi=300, width=6.5, height=6)

bryos_bold_genbank %>%
  select(FinalID, Family, Genus, Species, TaxonomicLevelCode,
         starts_with("GB_")) %>%
  mutate(GB_any = (GB_18S_count>0 |  GB_chloro_count  >0 | GB_mito_count  > 0 |  GB_mito_fl_count >0) ) %>%
  filter(!is.na(GB_any)) %>%
  filter(TaxonomicLevelCode<66) %>%
  group_by(TaxonomicLevelCode, GB_any) %>% 
  tally() 


bryos_bold_genbank %>%
  select(FinalID, Family, Genus, Species, TaxonomicLevelCode,
         starts_with("GB_"), starts_with("BOLD_")) %>%
  transmute(FinalID, Family, Genus, Species, TaxonomicLevelCode,
            GB_any = (GB_18S_count>0 |  GB_chloro_count  >0 | GB_mito_count  > 0 |  GB_mito_fl_count >0),
            BOLD_any = BOLD_Records>0) %>%
  filter(!is.na(GB_any)) %>%
  filter(TaxonomicLevelCode<66)
  
  
  