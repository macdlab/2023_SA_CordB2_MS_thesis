library(ggplot2)
library(ggplotify)
library(ggpubr)
library(gridExtra)
library(tibble)
library(tidyverse)

#expression plots


#annotation for datasets- blue, brown, greenyellow
#blue, brown, greenyellow

#blue
#read in the dataset
blue1=read.csv("BLUE_7_B2.csv")
#column to row for child_no
blue1=column_to_rownames(blue1, var = "gene_symbol")
#find the median
median(blue1$cord_b2)
#median comes 299
#annotate, H for more than 299 , and L for less than 299
blue1$Range[blue1$cord_b2<299]="L"
blue1$Range[blue1$cord_b2>=299]="H"


#brown
#read in the dataset
brown=read.csv("BROWN_7_B2.csv")
#column to row for child_no
brown=column_to_rownames(brown, var = "gene_symbol")
#find the median
median(brown$cord_b2)
#median comes 299
#annotate, H for more than 299 , and L for less than 299
brown$Range[brown$cord_b2>=299]="H"
brown$Range[brown$cord_b2<299]="L"

#greenyellow
#read in the dataset
greenyellow=read.csv("GREENYELLOW_7_B2.csv")
#column to row for child_no
greenyellow=column_to_rownames(greenyellow, var = "gene_symbol")
#find the median
median(greenyellow$cord_b2)
#median comes 299
#annotate, H for more than 299 , and L for less than 299
greenyellow$Range[greenyellow$cord_b2>=299]="H"
greenyellow$Range[greenyellow$cord_b2<299]="L"

#boxplot

#blue module

NOP14<- ggplot(blue1, aes(x = fct_rev(Range), y = NOP14)) + geom_boxplot(stat = "boxplot", fill = c("lightblue", "pink")) + labs(x = "cord B2 levels")
NOP14

DDX54<- ggplot(blue1, aes(x = fct_rev(Range), y = DDX54)) + geom_boxplot(stat = "boxplot", fill = c("lightblue", "pink")) + labs(x = "cord B2 levels")
DDX54

BOP1<- ggplot(blue1, aes(x = fct_rev(Range), y = BOP1)) + geom_boxplot(stat = "boxplot", fill = c("lightblue", "pink")) + labs(x = "cord B2 levels")
BOP1

NOC2L<- ggplot(blue1, aes(x = fct_rev(Range), y = NOC2L)) + geom_boxplot(stat = "boxplot", fill = c("lightblue", "pink")) + labs(x = "cord B2 levels")
NOC2L

TBL3<- ggplot(blue1, aes(x = fct_rev(Range), y = TBL3)) + geom_boxplot(stat = "boxplot", fill = c("lightblue", "pink")) + labs(x = "cord B2 levels")
TBL3

#multiple boxplots in one graph
library(ggpubr)
figure_blue<-ggarrange(DDX54, NOP14, BOP1, NOC2L, TBL3, labels = c("1", "2", "3", "4", "5"), ncol = 3, nrow = 2)
annotate_figure(figure_blue, right = "H = high levels of B2. L = low levels of B2")


#brown

FTSJ3<- ggplot(brown, aes(x = fct_rev(Range), y = FTSJ3)) + geom_boxplot(stat = "boxplot", fill = c("lightblue", "pink")) + labs(x = "cord B2 levels")
FTSJ3

SDAD1<- ggplot(brown, aes(x = fct_rev(Range), y = SDAD1)) + geom_boxplot(stat = "boxplot", fill = c("lightblue", "pink")) + labs(x = "cord B2 levels")
SDAD1

GTPBP4<- ggplot(brown, aes(x = fct_rev(Range), y = GTPBP4)) + geom_boxplot(stat = "boxplot", fill = c("lightblue", "pink")) + labs(x = "cord B2 levels")
GTPBP4

DDX10<- ggplot(brown, aes(x = fct_rev(Range), y = DDX10)) + geom_boxplot(stat = "boxplot", fill = c("lightblue", "pink")) + labs(x = "cord B2 levels")
DDX10

NOL10<- ggplot(brown, aes(x = fct_rev(Range), y = NOL10)) + geom_boxplot(stat = "boxplot", fill = c("lightblue", "pink")) + labs(x = "cord B2 levels")
NOL10


#multiple boxplots in one graph
library(ggpubr)
figure_brown<-ggarrange(FTSJ3, SDAD1, GTPBP4, DDX10, NOL10, labels = c("1", "2", "3", "4", "5"), ncol = 3, nrow = 2)
annotate_figure(figure_brown, right = "H = high levels of B2. L = low levels of B2")



#greenyellow module

HCK<- ggplot(greenyellow, aes(x = fct_rev(Range), y = HCK)) + geom_boxplot(stat = "boxplot", fill = c("lightblue", "pink")) + labs(x = "cord B2 levels")
HCK

SYK<- ggplot(greenyellow, aes(x = fct_rev(Range), y = SYK)) + geom_boxplot(stat = "boxplot", fill = c("lightblue", "pink")) + labs(x = "cord B2 levels")
SYK

RAC1<- ggplot(greenyellow, aes(x = fct_rev(Range), y = RAC1)) + geom_boxplot(stat = "boxplot", fill = c("lightblue", "pink")) + labs(x = "cord B2 levels")
RAC1

GRB2<- ggplot(greenyellow, aes(x = fct_rev(Range), y = GRB2)) + geom_boxplot(stat = "boxplot", fill = c("lightblue", "pink")) + labs(x = "cord B2 levels")
GRB2

PIK3CB<- ggplot(greenyellow, aes(x = fct_rev(Range), y = PIK3CB)) + geom_boxplot(stat = "boxplot", fill = c("lightblue", "pink")) + labs(x = "cord B2 levels")
PIK3CB

#multiple boxplots in one graph
figure_greenyellow<-ggarrange(PIK3CB, GRB2, RAC1, SYK, HCK, labels = c("1", "2", "3", "4", "5"), ncol = 3, nrow = 2)
annotate_figure(figure_greenyellow, right = "H = high levels of B2. L = low levels of B2")

