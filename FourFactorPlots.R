########################
## A R program to plot two differential condition volcano plots in four quadrants. 
## Usage. Rscript FourFactorPlots.R input1.xlsx input2.xlsx condition_1 condition_2
#########################
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(openxlsx)
library(readxl)
library(ggrepel)

#Arguement passing
args = commandArgs(trailingOnly=TRUE)

# Setting work directory
getwd()
setwd(getwd())

# Loading matrix of highconfidence genes
file1= args[1] # for ex. resN2_F180del_Input_results.csv
file2= args[2] # for ex. resN2_F180del_IP_results.csv


difcondf = read.xlsx(file1, sheet = as.numeric(args[1])) #read.table(file1, header = 1, sep = "\t") #args[1]
difconds = read.xlsx(file2, sheet = as.numeric(args[2])) #read.table(file2, header = 1, sep = "\t") #args[2]

upfirst = nrow(difcondf %>% filter(logFC >= 1 & FDR < 0.05))
dwfirst = nrow(difcondf %>% filter(logFC >= -1 & FDR < 0.05))

upsecond = nrow(difconds %>% filter(logFC >= 1 & FDR < 0.05))
dwsecond = nrow(difconds %>% filter(logFC >= -1 & FDR < 0.05))

list(upfirst, dwfirst, upsecond, dwsecond)

#m1dw = read.table("XV102vsXV33_CD_NewCountWS279_Mat2x2.csv.diff_2.xls", header = 1, sep = ",")
#m1up = read.table("XV139vsXV33_CD_NewCountWS279_Mat3x2.csv.diff_2.xls", header = 1, sep = ",")

#head(difcondf)
#head(difconds)

# Setting differential conditions to add in the dataframe
#dif1 = str_split(file1, "_")[[1]]
dif1_col = args[3] #"IL6_YvsN_UNX_N_F" #paste(dif1[[2]],dif1[[3]], sep = "_")
dif1lfc = paste('logFC',dif1_col,sep="_")

#dif2 = str_split(file2, "_")[[1]]
dif2_col = args[4] #"IL6_YvsN_UNX_Y_F" #paste(dif2[[2]],dif2[[3]], sep = "_")
dif2lfc = paste('logFC',dif2_col,sep="_")

list(dif1_col, dif2_col)




m1dw_1 = tibble(difcondf) %>% select(genes, logFC, FDR) %>%
  mutate(diff_G=dif1_col) #%>% mutate(genes=gsub("cel-","",genes)) 
head(m1dw_1)
#m1dw_1 %>% filter(genes %in% c('hbl-1','lin-14','lin-28','lin-29','lin-41','lin-42','lin-46')) %>% mutate(fold_change_XV102=log)

m1up_1 = tibble(difconds) %>% select(genes, logFC, FDR) %>%
  mutate(diff_G=dif2_col) %>% mutate(genes=gsub("^cel-","",genes))
head(m1up_1)

metadiff = inner_join(m1dw_1, m1up_1, by=c( "genes" = "genes"))


#%>% mutate(dif_statusFDR.x <= 0.05 & FDR.y <= 0.05) %>%
#  mutate(dif_status=case_when(logFC.x >= 0.5 & logFC.y >=0.5 ~ "Up", logFC.x <= -0.5 & logFC.y <= -0.5 ~ "Down")) %>%
#  replace_na(list(dif_status="NS"))
tail(metadiff)
dim(metadiff)
#nonsignficantGenes = metadiff %>% filter(dif_status %in% c("NS"))
#write.table(nonsignficantGenes, file = args[4], sep = "\t", quote = FALSE, row.names = FALSE)
#metadiff %>% filter(genes %in% c('hbl-1','lin-14','lin-28','lin-29','lin-41','lin-42','lin-46'))

# Add colour, size and alpha (transparency) to volcano plot --------------------
cols <- c("Up" = "#ffad73", "Down" = "#26b3ff", "NS" = "gray") 
sizes <- c("Up" = 2, "Down" = 2, "NS" = 1) 
alphas <- c("Up" = 1, "Down" = 1, "NS" = 0.5)

# Setting points for diff genes
upsigenesboth = metadiff %>% filter((logFC.x >= 1 & logFC.y >=1) & (FDR.x <= 0.05 & FDR.y <= 0.05)) %>% add_column(dif_status="Up (both)")
# mutate(dif_status=case_when(logFC.x >= 0.5 & logFC.y >=0.5 & FDR.x <= 0.05 ~ "Up")) #, logFC.x <= -0.5 & logFC.y <= -0.5 ~ "Down"))
dim(upsigenesboth)[1]
#filter(dif_status %in% c("Up"))

dwsigenesboth = metadiff %>% filter((logFC.x <= -1 & logFC.y <= -1) & (FDR.x <= 0.05 & FDR.y <= 0.05)) %>% add_column(dif_status="Down (both)")
#mutate(dif_status=logFC.x <= -0.5 & logFC.y <= -0.5 ~ "Down"))) #filter(dif_status %in% c("Down"))
dwsigenesboth

upinfirst = metadiff %>% filter(logFC.x >= 1 & FDR.x <= 0.05) %>% filter(!genes %in% upsigenesboth$genes) %>% add_column(dif_status=paste0("Up (",dif1_col,")"))
upinfirst

dwinfirst = metadiff %>% filter(logFC.x <= -1 & FDR.x <= 0.05) %>% filter(!genes %in% dwsigenesboth$genes) %>% add_column(dif_status=paste0("Down (",dif1_col,")"))
dwinfirst

#second condition
upinsecond = metadiff %>% filter(logFC.y >= 1 & FDR.y <= 0.05) %>% filter(!genes %in% upsigenesboth$genes) %>% add_column(dif_status=paste0("Up (",dif2_col,")"))
upinsecond

dwinsecond = metadiff %>% filter(logFC.y <= -1 & FDR.y <= 0.05) %>% filter(!genes %in% dwsigenesboth$genes) %>% add_column(dif_status=paste0("Down (",dif2_col,")"))
dwinsecond

nssigene = metadiff %>% filter(FDR.x > 0.05 & FDR.y > 0.05) %>% add_column(dif_status="NS")
nssigene

cldf = rbind(upsigenesboth, upinfirst, upinsecond, dwsigenesboth, dwinfirst, dwinsecond, nssigene)
head(cldf)

metadiffv1 = full_join(metadiff, cldf, by=c( "genes" = "genes")) %>% select(1:7,14) %>% 
  rename('logFC.x' = 'logFC.x.x', 'logFC.y' = 'logFC.y.x' )
head(metadiffv1)

#write.table(metadiffv1, file =args[4], quote = FALSE, sep = "\t", row.names = FALSE)

lmin = min(metadiff$logFC.x)
lmin1 = min(metadiff$logFC.y)
mlimlow = mean(lmin,lmin1)
mlimlow


lmax = max(metadiff$logFC.x)
lmax1 = max(metadiff$logFC.y)
mlimhigh = mean(lmax,lmax1)
minlim = mlimhigh + mlimlow
tlim = (minlim) - (mlimlow)
-tlim

names(metadiffv1)

#png(file=args[3], width = 8, height = 6, res=1200, pointsize = 1, units="in", bg="white")
#png(file="F180del_InputIP.HCGv2.png", width = 8, height = 6, res=1200, pointsize = 1, units="in", bg="white")

str_split(metadiffv1$genes,"_")[[2]]

metadiffv1 %>% mutate(label = gsub("ENSMUSG\\d+\\.\\d+\\_","",genes)) %>% filter(grepl("Up|Down", dif_status)) 

ggplot(data=metadiffv1, aes(x=logFC.x, y=logFC.y, color=dif_status)) + geom_point(alpha=0.1, size=3) +
  geom_point(data = upsigenesboth,shape = 16,size = 5,fill = "firebrick", colour="firebrick", alpha=0.5) +
  geom_point(data = upinfirst,shape = 16,size = 2,fill = "red", colour="red", alpha=0.5) +
  geom_point(data = upinsecond,shape = 16,size = 2,fill = "orange", colour="orange", alpha=0.5) +
  geom_point(data = dwsigenesboth,shape = 16,size = 5,fill = "steelblue", colour="steelblue", alpha=0.5) +
  geom_point(data = dwinfirst,shape = 16,size = 2,fill = "blue", colour="blue", alpha=0.5) +
  geom_point(data = dwinsecond,shape = 16,size = 2,fill = "skyblue", colour="skyblue", alpha=0.5) +
  geom_point(data = nssigene,shape = 16,size = 1,fill = "gray78", colour="gray78", alpha=0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color="gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color="gray") +
  geom_hline(yintercept = -1, linetype = "dashed", color="black") +
  geom_hline(yintercept = 1, linetype = "dashed", color="black") + 
  geom_vline(xintercept = -1, linetype = "dashed", color="black") +
  geom_vline(xintercept = 1, linetype = "dashed", color="black") + 
  xlim((-tlim), mlimhigh) + ylim((-tlim), mlimhigh) + theme_bw() +
  geom_text_repel(
    data = . %>% mutate(label = gsub("ENSMUSG\\d+\\.\\d+\\_","",genes)) %>% filter(grepl("Up|Down", dif_status)),
    aes(label = label),
    box.padding = .2,
    show.legend = FALSE,
    color = "black",
    alpha = 0.9,
    size = 3.5, fontface="bold"
  ) +
  scale_color_manual(name="Expression", breaks = c("Up (both)", paste0("Up (",dif1_col,")"), paste0("Up (",dif2_col,")"), "Down (both)", paste0("Down (",dif1_col,")"),paste0("Down (",dif2_col,")"), "NS"), 
                     values=c("firebrick", "red", "orange","steelblue", 'blue', 'skyblue','gray78'),
                     labels= c(paste0("Up (both)", "[",dim(upsigenesboth)[1],"]"), 
                               paste0("Up (",dif1_col,") ", "[",nrow(upinfirst),"]"), 
                               paste0("Up (",dif2_col,")", "[",nrow(upinsecond),"]"), 
                               paste0("Down (both)", "[",nrow(dwsigenesboth),"]"), 
                               paste0("Down (",dif1_col,")", "[",nrow(dwinfirst),"]"),
                               paste0("Down (",dif2_col,")","[",nrow(dwinsecond),"]"), 
                               "NS FDR > 0.05")) +
  labs(title = paste0("Differential Expressed Genes in ", dif1_col, " and ",dif2_col),
       x = paste0("logFC_", dif1_col),
       y = paste0("logFC_", dif2_col),
       colour = "Expression \nchange") +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  theme(plot.title = element_text(hjust=0.5, face = "bold", size=10)) +
  theme(axis.title.x = element_text(hjust=0.5, vjust=-0.5, face = "bold", size=10), 
        axis.title.y = element_text(hjust=0.5, vjust=-0.5, face = "bold", size=10)) +
  theme(axis.text.x = element_text(size=10, angle = 0, hjust=1, face = "bold")) +
  theme(axis.text.y = element_text(size=10, angle = 0, hjust=1, face = "bold")) +
  theme(legend.title = element_text(face="bold", size=10, vjust=1), legend.text = element_text(face="bold", size=10, angle=0),
        legend.key.size=unit(2, "line")) 


ggsave(file = paste(dif1_col,dif2_col,"FourFactor","png",sep = "."),  plot = last_plot(),  device = "png",  width = 8,  height = 6,  units = "in",  bg = "white",  scale = 1.2,
       dpi = 600)

# Modify legend labels by re-ordering gene_type levels -------------------------
#metadiff <- metadiff %>%
# mutate(dif_status = fct_relevel(dif_status, "Up", "Down", "NS")) 



#metadiff$dif_status = factor(metadiff$dif_status, levels=unique(metadiff$dif_status), ordered=T)

#ggplot(data = metadiff, aes(x=logFC.x, y=logFC.y, color=dif_status)) + geom_point(alpha=.5, size=1.9) + 
#    scale_color_manual(name="Expression", breaks = c("Up", "NS", "Down"), values=c("firebrick", "gray", "steelblue")) +
#  xlim((-tlim), mlimhigh) + ylim((-tlim), mlimhigh) +
#  theme_bw() + theme(text=element_text(size=10, face = "bold"))  +
#  geom_hline(yintercept = 0, linetype = "dashed", color="gray") +
#  geom_vline(xintercept = 0, linetype = "dashed", color="gray") +
#  geom_hline(yintercept = -0.5, linetype = "dashed", color="black") +
#  geom_hline(yintercept = 0.5, linetype = "dashed", color="black") + 
#  geom_vline(xintercept = -0.5, linetype = "dashed", color="black") +
#  geom_vline(xintercept =0.5, linetype = "dashed", color="black") +
#    labs(title = paste0("Differential Expressed miRNAs in ", dif1[[2]], " Input/IP"),
#         x = paste0("log2FC_", dif1_col),
#         y = paste0("log2FC_", dif2_col),
#         colour = "Expression \nchange") +
#    theme_bw() + # Select theme with a white background  
#    theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
#          panel.grid.minor = element_blank(),
#          panel.grid.major = element_blank()) +
#    theme(plot.title = element_text(hjust=0.5, face = "bold", size=10)) +
#    theme(axis.title.x = element_text(hjust=0.5, vjust=-0.5, face = "bold", size=10), 
#          axis.title.y = element_text(hjust=0.5, vjust=-0.5, face = "bold", size=10)) +
#    theme(axis.text.x = element_text(size=10, angle = 0, hjust=1, face = "bold")) +
#    theme(axis.text.y = element_text(size=10, angle = 0, hjust=0, face = "bold")) +
#    theme(legend.title = element_text(face="bold", size=10, vjust=1), legend.text = element_text(face="bold", size=8, angle=0),
#          legend.key.size=unit(.8, "line")) 

#dev.off()
#ggsave("V254I_InputIP.HCGv1.png", width = 8, height = 6, units = "in", bg = "white", dpi = 300, scale = 1)

#geom_point(data = dwsigenesboth,shape = 16,size = 3,fill = "steelblue", colour="steelblue", alpha=0.5) +
#  geom_point(data = upsigenesboth,shape = 16,size = 3,fill = "firebrick", colour="firebrick", alpha=0.5) +
#  geom_point(data = dwinsecond,shape = 16,size = 2,fill = "skyblue", colour="skyblue", alpha=0.5) +
#  geom_point(data = upinfirst,shape = 16,size = 2,fill = "red", colour="red", alpha=0.5) +
#  geom_point(data = nssigene,shape = 16,size = 1,fill = "gray", colour="gray", alpha=0.5) +
#  geom_point(data = dwinfirst,shape = 16,size = 2,fill = "blue", colour="blue", alpha=0.5) +
#  geom_point(data = upinsecond,shape = 16,size = 2,fill = "orange", colour="orange", alpha=0.5) +
