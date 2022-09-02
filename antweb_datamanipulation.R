library(tidyverse)

ants_df<-bind_rows(
  read_delim("Data/AntWebLists/AZ_Ants_speciesListDownload.txt") %>% mutate(State="Arizona"),
  read_delim("Data/AntWebLists/CA_Ants_speciesListDownload.txt") %>% mutate(State="California"),
  read_delim("Data/AntWebLists/BJMX_Ants_speciesListDownload.txt") %>% mutate(State="BajaCalifornia"),
  read_delim("Data/AntWebLists/CO_Ants_speciesListDownload.txt") %>% mutate(State="Colorado"),
  read_delim("Data/AntWebLists/NM_Ants_speciesListDownload.txt") %>% mutate(State="NewMexico"),
  read_delim("Data/AntWebLists/NV_Ants_speciesListDownload.txt") %>% mutate(State="Nevada"),
  read_delim("Data/AntWebLists/SONMX_Ants_speciesListDownload.txt") %>% mutate(State="Sonora"),
  read_delim("Data/AntWebLists/UT_Ants_speciesListDownload.txt") %>% mutate(State="Utah")
) %>%
  rename(Author_Date=`Author Date`)


ants_df %>%
  mutate(Present=1) %>%
  pivot_wider(names_from = State,
              values_from = Present,
              values_fill = 0) %>%
  mutate(GenSp = paste(Genus, Species),
         GenSpSp = case_when(is.na(Subspecies)~NA_character_ ,
                             T~paste(Genus, Species, Subspecies))
  ) %>%
  # as.data.frame() %>% head()
  bind_rows(
    #Genus summary
    ants_df %>%
      mutate(Present=1) %>%
      group_by(Subfamily, Genus, State) %>% tally(name="Taxa") %>%
      pivot_wider(names_from = State,
                  values_from = Taxa,
                  values_fill = 0),
    #Subfamily summary
    ants_df %>%
      mutate(Present=1) %>%
      group_by(Subfamily,  State) %>% tally(name="Taxa") %>%
      pivot_wider(names_from = State,
                  values_from = Taxa,
                  values_fill = 0)
  
  ) %>%
  mutate(TaxonomicLevelCode = case_when(is.na(Genus)~42,
                                        is.na(Species)~50,
                                        is.na(Subspecies)~60,
                                        T~66),
         FinalID = case_when(TaxonomicLevelCode==66~GenSpSp,
                             TaxonomicLevelCode==60~GenSp,
                             TaxonomicLevelCode==50~Genus,
                             TaxonomicLevelCode==42~Subfamily,
                             T~NA_character_),
         ParentFinalID = case_when(TaxonomicLevelCode==66~GenSp,
                                   TaxonomicLevelCode==60~Genus,
                                   TaxonomicLevelCode==50~"tribe",
                                   TaxonomicLevelCode==42~"Formicidae",
                                   T~NA_character_)
         ) %>%
  
  transmute(FinalID, ParentFinalID, 
            Active= -1,
            Phylum="Arthropoda", Subphylum="Hexapoda", 
            Class="Insecta", Order="Hymenoptera", Family="Formicidae",
            Subfamily,
            # Tribe=NA_character_,
            Genus, Species=GenSp, Subspecies=GenSpSp,
            TaxonomicLevelCode, TaxaList=0,LifeStageDefault="A",
            CaliforniaTaxon=0,
            TaxonomicAuthority=Author_Date,
            FinalIDAuthority="CSUMB-WEE", Source="CSUMB-WEE",
            Notes=case_when(California==1~"Found in CA",T~"Found in SW USA"),
            California, Arizona,  BajaCalifornia, Sonora, Nevada,NewMexico, Colorado, Utah
            ) %>%
  # as.data.frame() %>% head()
  write_csv(file="Data/AntWebLists/AntWeb_Consolidated.csv")
