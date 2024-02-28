install.packages("VennDiagram")
library(VennDiagram)

base_proteins <- read.csv("Data/venndiagram.csv")
View(base_proteins)

if (!require(devtools)) install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram")
library(ggVennDiagram)
library(ggplot2)

  citosol <- c("Factor_E.","Acil_CoA", "ZapE")
  membrana <-  c("SubB_ATPs","DUF3068") 
  division <- ("ZapE") 
  pared <- c("IniB", "PE_PGRS")
  metabolismo <- c("Acil_CoA","SubB_ATPs","DUF3068","Diacilglicerol") 
  superficie <- ("PTP")
  secretion <- c("MPB64_MPT64", "PE_PGRS","Diacilglicerol","PTP")


ggVennDiagram(x)

####
install.packages("ggvenn")

library(ggvenn)
library(RColorBrewer)
citosol <- c("Factor_E.","Acil_CoA", "ZapE")
membrana <-  c("SubB_ATPs","DUF3068") 
division <- c("ZapE") 
pared <- c("IniB", "PE_PGRS")
metabolismo <- c("Acil_CoA","SubB_ATPs","DUF3068","Diacilglicerol") 
superficie <- c("PTP")
secretion <- c("MPB64_MPT64", "PE_PGRS","Diacilglicerol","PTP")


x <- list(citosol=citosol, membrana=membrana, division=division, pared=pared, metabolismo=metabolismo, superficie=superficie,secretion=secretion)

ggvenn(x, show_elements = T, label_sep = "\n")


##############################
if (!"venn" %in% installed.packages()) {
  install.packages("venn")
}
if (!"ggplot2" %in% installed.packages()) {
  install.packages("ggplot2")
}
if (!"ggpolypath" %in% installed.packages()) {
  install.packages("ggpolypath")
}

library(venn)
library(ggplot2)
library(ggpolypath)


Proteinas <- read.csv("Data/Libro1.csv")
Proteinas
Proteinas2 <- Proteinas[,-1]
row.names(Proteinas2) <- Proteinas[,1]
#samp2 <- samp[,-1]
#rownames(samp2) <- samp[,1]
Proteinas
t(Proteinas2)

 
venn(Proteinas2, ggplot = T)

#############################


library(tidyverse)
#install.packages("ggvenn")
library(ggvenn)
#> Loading required package: grid

Name <- c("Citoplasma", "Metabolismo", "Extracelular", "Division", "Membrana", "Pared","Superficie")
IniB <- c(0,0,0,0,0,1,0)
FactorE <- c(1,0,0,0,0,0,0)
MPB64_MPT64 <- c(0,0,1,0,0,0,0)
PE_PEGRS <- c(0,0,1,0,0,1,0)
Acil_CoA <- c( 1, 1,0,0,0,0,0)
SubB_ATPs <- c(0,1,0,0,1,0,0)
DUF3068 <- c(0,1,0,0,1,0,0)
DGAT <- c(0,1,1,0,0,0,0)
ZapE <- c(1,0,0,1,0,0,0)
PTP <- c(0,0,1,0,0,0,1)

df1 <- data.frame(IniB, FactorE, MPB64_MPT64, PE_PEGRS, Acil_CoA, SubB_ATPs, DUF3068, DGAT,ZapE, PTP, Name)

df1 %>%
  mutate(across(starts_with("Event"), as.logical)) %>%
  ggplot() +
  geom_venn(aes(A = Event1, B = Event2, C = Event3, D = Event4),
            set_names = Name)


