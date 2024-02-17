pdf("/data/lucas/tanner/figures/manhattans_survival_MF_tanner34.pdf",  
    width=6, height=5.5)


library(graphics)
library(dplyr)
library(data.table)

layout(matrix(c(1,2), nrow = 2, ncol = 1, byrow = TRUE))
par(cex=0.6) # size of chart
par(mai=c(0.3,0.6,0.3,0.3)) # adjusts thickness of margins (white space)
par(oma=c(0.5,0.3,0.3,0.3)) 

################################################
################################################
################################################

# Boys:

chrs=fread("/data/genetica/Tanner/tanner34_boys_survival.csv",
           header=T, sep = ",")
colnames(chrs)=c("CHR","SNP","Beta","std_error","Z_val","P_val")
snps=fread("/data/genetica/Tanner/snps_filtrados_tanner34_boys_survival.csv",
           header=T, sep = ",")
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

# create color vector with 1 == red, and 2 == deepskyblue2
z$color <- ifelse(gr %% 2 == 0, "coral1", "deepskyblue4")

plot(z$x, z$y, type="p", cex = 0.8, pch = 16,
     col = z$color,
     lwd=0.1,
     frame.plot = F,
     xaxt = 'n', # removes x labels,
     xlim = c(1,count_chrs$LastSNP[22]),
     ylim = c(2, 12),
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
mtext(side = 3, text = expression(italic("LARS2")), 
      line = 1, adj = 0.23, padj = 10.7, las = 1, cex=0.7, col="coral1") 
mtext(side = 3, text = expression(italic("LIMD1")), 
      line = 1, adj = 0.23, padj = 14, las = 1, cex=0.7, col="coral1") 
mtext(side = 3, text = expression(italic("FAM83B")), 
      line = 1, adj = 0.45, padj = 11.7, las = 1, cex=0.7, col="deepskyblue4") 
mtext(side = 3, text = expression(italic("ZNF320")), 
      line = 1, adj = 0.99, padj = 13.0, las = 1, cex=0.7, col="coral1") 

################################################
################################################
################################################

# Girls:

chrs=fread("/data/genetica/Tanner/tanner34_girls_survival.csv",
           header=T, sep = ",")
colnames(chrs)=c("CHR","SNP","Beta","std_error","Z_val","P_val")
snps=fread("/data/genetica/Tanner/snps_filtrados_tanner34_girls_survival.csv",
           header=T, sep = ",")
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

# create color vector with 1 == red, and 2 == deepskyblue2
z$color <- ifelse(gr %% 2 == 0, "coral1", "deepskyblue4")

plot(z$x, z$y, type="p", cex = 0.8, pch = 16,
     col = z$color,
     lwd=0.1,
     #     xlab= "physical position",
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
mtext(side = 3, text = expression(italic("AL157359.3")), 
      line = 1, adj = 0.99, padj = 4.5, las = 1, cex=0.7, col="coral1") 

dev.off()

