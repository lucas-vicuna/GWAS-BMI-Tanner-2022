pdf("/data/lucas/tanner/figures/manhattans_survival_MF_tanner12.pdf",  
    width=6, height=5.5)

library(graphics)
library(dplyr)
library(data.table)

layout(matrix(c(1,2), nrow = 2, ncol = 1, byrow = TRUE))
par(cex=0.6)
par(mai=c(0.3,0.6,0.3,0.3)) # adjusts thickness of margins (white space)
par(oma=c(0.5,0.3,0.3,0.3)) 
################################################
################################################
################################################

# Boys:

chrs=fread("/data/genetica/Tanner/tanner12_boys_survival.csv",header=T)
colnames(chrs)=c("CHR","SNP","Beta","std_error","Z_val","P_val")
snps=fread("/data/genetica/Tanner/snps_filtrados_tanner12_boys_survival.csv",
           header=T)
chrs = semi_join(chrs,snps,c("CHR","SNP"))

# GOCS annots:

vep=fread("/data/lucas/selection_gocs/vep/gocs_vep.txt",header = T, stringsAsFactors = F)
vep = vep[,c("SNP","Location", "Consequence", "SYMBOL","BIOTYPE")]
vep$Location = gsub('\\-.*', '',vep$Location)
vep$Location = gsub(".*:","",vep$Location)
vep = vep[!duplicated(vep$SNP, fromLast=T), ]
colnames(vep)[2]="POS"

chrs2 = left_join(chrs,vep,"SNP") %>%
  group_by(CHR) %>%
  arrange(CHR,as.numeric(POS)) %>%
  ungroup() %>%
  mutate(SNP_pos = row_number()) 

#chrs2 = chrs2[-log10(chrs2$P_val) > 2,]

# Count SNPs per chr:
count_chrs = chrs2 %>% 
  group_by(CHR) %>% 
  summarise(n = n()) %>%
  mutate(LastSNP = cumsum(n),
         LastSNPplus1 = LastSNP+1,
         FirstSNP = lag(LastSNPplus1, default = first(LastSNPplus1)),
         FirstSNP2 = ifelse(LastSNPplus1!=FirstSNP,FirstSNP-1,1),
         chr_size_div2=n/2,
         middle_val=FirstSNP2+chr_size_div2) %>%
  select(-FirstSNP)

X2 = chrs2[-log10(chrs2$P_val) > 2,]

minPboys = min(X2$P_val)

#create data
x = X2$SNP_pos
y = -log10(X2$P_val)
z = data.frame(x,y)

#cut in segments
my_segments = count_chrs$LastSNP

gr <- cut(z$x, my_segments,labels = FALSE, right = T)
gr[is.na(gr)] <- 0

# create color vector 
z$color <- ifelse(gr %% 2 == 0, "coral1", "deepskyblue4")

plot(z$x, z$y, type="p", cex = 0.8, pch = 16,
     col = z$color,
     lwd=0.1,
     frame.plot = F,
     xaxt = 'n', # removes x labels,
     xlim = c(1,count_chrs$LastSNP[22]),
     ylim = c(2, 8),
     las = 2,
     cex.lab=1, # size of axis labels
     ann = FALSE, # remove axis titles
     mgp = c(3, 0.8, 0)) 

# adjust y axis label size
par(cex.axis= 1.1, tck=-0.03)

abline(h = 7.3, col= "black", lty=3, lwd=0.7)

mtext(side = 1, text = "Chromosome", line = 2.1, cex = 1) 
mtext(side = 2, 
      text = expression(paste("-log"[10]," ",italic("P"),"-value",sep="")),
      line = 1.7, adj = 0.5, cex = 1) 

options(scipen = 999) # unables scientific notation por axes

#Add ticks to X axis
axis(side=1, at = c(1,count_chrs$LastSNP), 
     labels = F, las=1, cex.axis=0.6, font = 2, tck=-0.03, pos = 1.85)

# Add labels between ticks
axis(side=1, count_chrs$middle_val, labels = seq(1,22,1), 
     las=2, cex.axis=1, font = 1, tck= F, lwd=0, mgp = c(3, 0.4, 0))

mtext(side = 3, text = "A", line = 1, adj = -0.115, padj = 0.75, las = 1, cex=1.8) # to adjust distance of text from axis

################################################
################################################
################################################

# Girls:

chrs=fread("/data/genetica/Tanner/tanner12_girls_survival.csv",header=T)
colnames(chrs)=c("CHR","SNP","Beta","std_error","Z_val","P_val")
snps=fread("/data/genetica/Tanner/snps_filtrados_tanner12_girls_survival.csv",
           header=T)
chrs = semi_join(chrs,snps,c("CHR","SNP"))

chrs2 = left_join(chrs,vep,"SNP") %>%
  group_by(CHR) %>%
  arrange(CHR,as.numeric(POS)) %>%
  ungroup() %>%
  mutate(SNP_pos = row_number()) 

# Count SNPs per chr:
count_chrs = chrs2 %>% 
  group_by(CHR) %>% 
  summarise(n = n()) %>%
  mutate(LastSNP = cumsum(n),
         LastSNPplus1 = LastSNP+1,
         FirstSNP = lag(LastSNPplus1, default = first(LastSNPplus1)),
         FirstSNP2 = ifelse(LastSNPplus1!=FirstSNP,FirstSNP-1,1),
         chr_size_div2=n/2,
         middle_val=FirstSNP2+chr_size_div2) %>%
  select(-FirstSNP)

X2 = chrs2[-log10(chrs2$P_val) > 2,]

minPgirls = min(X2$P_val)

#create data
x = X2$SNP_pos
y = -log10(X2$P_val)
z = data.frame(x,y)

#cut in segments
my_segments = count_chrs$LastSNP

gr <- cut(z$x, my_segments,labels = FALSE, right = T)
gr[is.na(gr)] <- 0

# create color vector
z$color <- ifelse(gr %% 2 == 0, "coral1", "deepskyblue4")

plot(z$x, z$y, type="p", cex = 0.8, pch = 16,
     col = z$color,
     lwd=0.1,
     frame.plot = F,
     xaxt = 'n', # removes x labels,
     xlim = c(1,count_chrs$LastSNP[22]),
     ylim = c(2, 8),
     las = 2,
     cex.lab=1, # size of axis labels
     ann = FALSE, # remove axis titles
     mgp = c(3, 0.8, 0)) 

# adjust y axis label size
par(cex.axis= 1.1, tck=-0.03)

abline(h = 7.3, col= "black", lty=3, lwd=0.7)

mtext(side = 1, text = "Chromosome", line = 2.1, cex = 1) 
mtext(side = 2, 
      text = expression(paste("-log"[10]," ",italic("P"),"-value",sep="")),
      line = 1.7, adj = 0.5, cex = 1) 

options(scipen = 999) # unables scientific notation por axes

#Add ticks to X axis
axis(side=1, at = c(1,count_chrs$LastSNP), 
     labels = F, las=1, cex.axis=0.6, font = 2, tck=-0.03, pos = 1.85)

# Add labels between ticks
axis(side=1, count_chrs$middle_val, labels = seq(1,22,1), 
     las=2, cex.axis=1, font = 1, tck= F, lwd=0, mgp = c(3, 0.4, 0))

mtext(side = 3, text = "B", line = 1, adj = -0.115, padj = 0.75, las = 1, cex=1.8) # to adjust distance of text from axis

dev.off()

