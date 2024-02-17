library(data.table)

# Plot global ancestry

setwd("/data/genetica/tomas_globalanc_gocs.mpu.aym.pel95.ibs.yri")

anc = fread("ga.chi_mpu_aym_pel_ibs_yri.k4.tsv", header = F)
colnames(anc) = c("POP","sampleID","AFR","AYM","EUR","MAP")

# Select pops:

chi = anc %>% filter(POP=="chi") %>% arrange(MAP)
ibs = anc %>% filter(POP=="ibs") %>% arrange(EUR)
aym = anc %>% filter(POP=="aym") %>% arrange(AYM)
mpu = anc %>% filter(POP=="mpu") %>% arrange(MAP)
yri = anc %>% filter(POP=="yri") %>% arrange(AFR)

pops = rbind(chi,aym,mpu,ibs,yri)

pops$index <- as.numeric(row.names(pops))
pops = pops[order(pops$index), ]
pops1 = pops[,3:6]

# Estimate global ancestries:

#eur_anc = mean(pops$EUR) #
#mpu_anc = mean(pops$MAP) # 
#aym_anc = mean(pops$AYM) # 

pdf("/data/lucas/tanner/figures/global_ancestry.pdf",  
    width=8, height=5)

#Draw outside plot area
par(xpd = TRUE)

barplot(t(as.matrix(pops1)),
        space = c(0),
        mgp = c(3, 0.6, 0),
        col=c("green4","orangered3","royalblue2","tan1"),
        border=NA,
        xaxt = 'n', # removes x labels,
        ann = FALSE)

segments(x0=1,x1=900,y0=-0.02,y1=-0.02,lwd=1,col="black") 
segments(x0=908,x1=973,y0=-0.02,y1=-0.02,lwd=1,col="black") 
segments(x0=979,x1=986,y0=-0.02,y1=-0.02,lwd=1,col="black") 
segments(x0=993,x1=1090,y0=-0.02,y1=-0.02,lwd=1,col="black") 
segments(x0=1100,x1=1203,y0=-0.02,y1=-0.02,lwd=1,col="black") 

text(x=450, y=-0.1, label="CHI", srt=90, cex=0.7)
text(x=935, y=-0.105, label="AYM", srt=90, cex=0.7)
text(x=982, y=-0.11, label="MAP", srt=90, cex=0.7)
text(x=1040, y=-0.1, label="IBS", srt=90, cex=0.7)
text(x=1150, y=-0.1, label="YRI", srt=90, cex=0.7)

mtext(side = 2, text = "Global ancestry proportions (k = 4)", 
      line = 1.7, adj = 0.65)

mtext(side = 1, text ="Subjects", line = 2.5) 

dev.off()




