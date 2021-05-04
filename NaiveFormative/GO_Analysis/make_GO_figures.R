# Ian Jones
# GGPLOT for Niko
# 9/23/20

library(ggplot2)
library(readxl)
library(dplyr)
library(forcats)

# Set working dir


# Read in data
df <- read_excel('for_analysis/David_Annotation_SIPs.xlsx')
naive <- df[df$`Cell Type`=='Naïve',]
form <- df[df$`Cell Type`=='Formative',]

class(naive)


# Make figures
n<-ggplot(naive, aes(x=reorder(`Functional Categories`, `Enrichment Score`), y=`Enrichment Score`)) +
    geom_bar(stat="identity", fill="darkseagreen3", width=0.5)+
    theme_classic() + 
    coord_flip() +
    ggtitle("Naïve terms") +
    xlab("") +
    ylab("Enrichment Score") +
    scale_y_continuous(position = "right") +
    theme(plot.title = element_text(hjust=0, size=14,face="bold"),
          axis.text=element_text(size=12))
n

f<-ggplot(form, aes(x=reorder(`Functional Categories`, `Enrichment Score`), y=`Enrichment Score`)) +
  geom_bar(stat="identity", fill="darkred", width=0.5)+
  theme_classic() + 
  coord_flip() +
  ggtitle("Formative terms") +
  xlab("") +
  ylab("Enrichment Score") +
  scale_y_continuous(position = "right") +
  theme(plot.title = element_text(hjust=0, size=14,face="bold"),
        axis.text=element_text(size=12))
f


df <- read_excel('for_analysis/David_Annotation_Form_Bivalent_Promoter.xlsx')
naive.biv <- df[df$`Cell Type`=='Naive Bivalent',]
naive.active <- df[df$`Cell Type`=='Naive Active',]

class(naive.biv)


# Make figures
n.biv<-ggplot(naive.biv, aes(x=reorder(`Functional Categories`, `Enrichment Score`), y=`Enrichment Score`)) +
    geom_bar(stat="identity", fill="#009E73", width=0.5)+
    theme_classic() + 
    coord_flip() +
    ggtitle("Naïve Bivalent terms") +
    xlab("") +
    ylab("Enrichment Score") +
    scale_y_continuous(position = "right") +
    theme(plot.title = element_text(hjust=0, size=14,face="bold"),
          axis.text=element_text(size=12))
n.biv

ggsave(n.biv,file='plots/David_Annotation_Form_Biv_subset_Naive_Biv_GO.pdf',height=7, width=7)

n.active<-ggplot(naive.active, aes(x=reorder(`Functional Categories`, `Enrichment Score`), y=`Enrichment Score`)) +
  geom_bar(stat="identity", fill="#0072B2", width=0.5)+
  theme_classic() + 
  coord_flip() +
  ggtitle("Naive Active terms") +
  xlab("") +
  ylab("Enrichment Score") +
  scale_y_continuous(position = "right") +
  theme(plot.title = element_text(hjust=0, size=14,face="bold"),
        axis.text=element_text(size=12))
n.active

ggsave(n.active,file='plots/David_Annotation_Form_Biv_subset_Naive_Active_GO.pdf',height=7, width=7)
