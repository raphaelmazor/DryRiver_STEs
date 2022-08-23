
# setwd("~/Documents/R/Genbank/bryophytes_2022/")

library(tidyverse)
library(bold)    # API interface to BOLD
library(taxize)  # for NCBI taxonomy lookup
library(seqinr)  # for FASTA output
library(rentrez)

options(ENTREZ_KEY = "b2ccdcd6619a29d3bf31de74e7cde9a1c209")

# import data / clean ---------------------------------
tax<-read_csv("Susie_082022/California Bryophyte STE.csv") %>%
  mutate(Species=trimSpace(Species))
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

write_csv(result.long, "Susie_082022/spellcheck.csv")

taxa<-sort(unique(result.long$matched_name2))

###################################################
# BOLD -------------------------------------------------
###################################################

# get specimen info from BOLD using our taxa list ----------

out1<-list()

for(i in 1:length(taxa)) { # i <- 4
  sp<-(taxa[i])
  #foo1<-bold_specimens(taxon = sp, response = T, timeout_ms = 10)
  foo2<-bold_stats(taxon = sp)
  out1[paste(sp)]<- paste0(sp, "_",foo2$total_records)
  print(paste(sp, foo2$total_records))
}
out.bold<-as.data.frame(t(t(out1)))
foo<-as.data.frame(do.call(rbind, strsplit(as.character(out.bold$V1), "_")))
names(foo)<-c("Species","Count")
out.bold<-foo

foo2x<-sapply(taxa, function(x){
  bold_x<-bold_stats(x)
})

bold_stats(taxa[3])
# california only -------------DOES NOT WORK. Returns *all* records from a UC institution, ignores taxa.


out1<-list()

for(i in 1:length(taxa)) { # i <- 4
  sp<-(taxa[i]) # sp<-"Drepanocladus aduncus"
  #foo1<-bold_specimens(taxon = sp, response = T, timeout_ms = 10)
  foo2<-bold_stats(taxon = sp, geo = "California")
  out1[paste(sp)]<- paste0(sp, "_",foo2$total_records)
  print(paste(sp, foo2$total_records))
}

out.bold<-as.data.frame(t(t(out1)))
foo<-as.data.frame(do.call(rbind, strsplit(as.character(out.bold$V1), "_")))
names(foo)<-c("Species","Count")
out.bold.ca<-foo
write.csv(as.matrix(out.bold.ca),"bold.ca.csv" )



###################################################
# GENBANK -------------------------------------------------------------
###################################################

#The function can't handle the whole data set, break into smaller pieces

# mito 
outGB<-list()
for (i in taxa) { # i <- "emiliania huxleyi"
  term1<-paste0(i,"[PORGN]", " AND 200:3400[SLEN]", " AND mitochondrion[filter]")
  foo<-entrez_search(db="nuccore", term=term1, FILT=1, api_key="b2ccdcd6619a29d3bf31de74e7cde9a1c209"  )
  outGB[i]<-foo$count
  print(paste(i, foo$count))}

outGB.mito<-outGB

outGB<-list()
for (i in taxa) { # i <- "emiliania huxleyi"
  term1<-paste0(i,"[PORGN]", " AND 200:34000000[SLEN]", " AND mitochondrion[filter]")
  foo<-entrez_search(db="nuccore", term=term1, FILT=1, api_key="b2ccdcd6619a29d3bf31de74e7cde9a1c209"  )
  outGB[i]<-foo$count
  print(paste(i, foo$count))}

outGB.mito.fl<-outGB

# chl 
outGB<-list()
for (i in taxa) { # i <- "emiliania huxleyi"
  term1<-paste0(i,"[PORGN]", " AND 200:3400[SLEN]", " AND chloroplast[filter]")
  foo<-entrez_search(db="nuccore", term=term1, FILT=1, api_key="b2ccdcd6619a29d3bf31de74e7cde9a1c209"  )
  outGB[i]<-foo$count
  print(paste(i, foo$count))}

outGB.chl<-outGB

# 18S

outGB<-list()
for (i in taxa) { # i <- "emiliania huxleyi"
  term1<-paste0(i,"[PORGN]", " AND 200:3400[SLEN]", " AND 18S[All fields]")
  foo<-entrez_search(db="nuccore", term=term1, FILT=1, api_key="b2ccdcd6619a29d3bf31de74e7cde9a1c209"  )
  outGB[i]<-foo$count
  print(paste(i, foo$count))}

outGB.18S<-outGB



###################################################
# GEOMDB -----------------------
###################################################
library(geomedb)

outgeomedb<-list()
  
  for(i in 1:length(taxa)) { # i<-2
    sp<-droplevels(taxa[i])
    foo1<-do.call(rbind, strsplit((as.character(sp)), " "))
    j<-foo1[1]
    k<-foo1[2]
    if (is.na(k) == T) {query1<-paste("genus =", j) }
    if (is.na(k) == F) {query1<-paste("genus =", j, "AND specificEpithet =", k)}
    data1 <- querySanger(locus = 'CO1', query = query1 )
    outgeomedb[i]<-(paste(j,k, (data1)))
    rm(data1)}


### 
a<-as.data.frame(t(t(outGB.18S)))
b<-as.data.frame(t(t(outGB.mito)))
c<-as.data.frame(t(t(outGB.mito.fl)))
d<-as.data.frame(t(t(outGB.chl)))
e<-as.data.frame(t(t(outgeomedb)))
f<-out.bold

out<-cbind(a,b,c,d,e,f)
names(out) <- c("GB.18S", "GB.mito", "GB.mito.fl","out.GB.chl", "GEOME", "BOLD", "BOLD.ct")

write.csv(as.matrix(out), "out.CA.csv")





