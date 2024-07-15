#----------------------- PACKAGES AND FUNCTIONS --------------------------------

library(r2vcftools)
library(LEA)
library(tess3r) 
library(seqinr)
library(ade4)
library(adegenet)
library(vcfR)
library(SNPRelate)
library(dartR)
library(devtools)
library(usethis)
library(poppr)
library(hierfstat)
library(reshape2)
library(ggplot2)
library(dplyr)
library(readxl)
library(SNPfiltR)
library(r2vcftools)

VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}

GenDiv <- function(vcf){
  HE <- Query(vcf, type="het")
  
  HE$HO <- (HE$N_SITES-HE$O.HOM.)/(HE$N_SITES) ## Observed heterozygosity (HO)
  error_ho <- qt(0.975,df=length(HE$HO)-1)*sd(HE$HO)/sqrt(length(HE$HO))
  ho <- mean(HE$HO)
  left_ho <- ho-error_ho
  right_ho <- ho+error_ho
  HO.df <- data.frame(He_O=ho, low_He_O=left_ho, up_He_O=right_ho)
  
  HE$HE <- (HE$N_SITES-HE$E.HOM.)/(HE$N_SITES) ## Expected heterozygosity (HE)
  error_he <- qt(0.975,df=length(HE$HE)-1)*sd(HE$HE)/sqrt(length(HE$HE))
  he <- mean(HE$HE)
  left_he <- he-error_he
  right_he <- he+error_he
  HE.df <- data.frame(He_E=he, low_He_E=left_he, up_He_E=right_he)
  
  error_f <- qt(0.975,df=length(HE$F)-1)*sd(HE$F)/sqrt(length(HE$F))
  f <- mean(HE$F)
  left_f <- f-error_f
  right_f <- f+error_f
  F.df <- data.frame(F=f, low_F=left_f, up_F=right_f)
  
  PI <- Query(vcf, type="site-pi")
  error_pi <- qt(0.975,df=length(PI$PI)-1)*sd(PI$PI)/sqrt(length(PI$PI))
  pi <- mean(PI$PI)
  left_pi <- pi-error_pi
  right_pi <- pi+error_pi
  PI.df <- data.frame(PI=pi, low_PI=left_pi, up_PI=right_pi)
  
  print(paste0("OUTPUT NUCLEOTIDE DIVERGENCE STATISTICS")) ## DIVERGENCE STATISTICS
  print(paste0("Observed heterozygosity (He_O)"))
  print(paste0("Expected heterozygosity (He_E)"))
  print(paste0("Coefficient of inbreeding (F)"))
  print(paste0("Measures nucleotide divergency on a per-site basis (PI)"))
  print(paste0("The 95% confidence interval - Lower limit (low_)"))
  print(paste0("The 95% confidence interval - Upper limit (up_)"))
  RES <- cbind(HO.df,  HE.df,  F.df,  PI.df)
  return(RES)
}

PlotK <- function(snmfproject){
  
  Mk1 <- mean(cross.entropy(snmfproject, K=1))
  SDk1 <- sd(cross.entropy(snmfproject, K=1))
  Mk2 <- mean(cross.entropy(snmfproject, K=2))
  SDk2 <- sd(cross.entropy(snmfproject, K=2))
  Mk3 <- mean(cross.entropy(snmfproject, K=3))
  SDk3 <- sd(cross.entropy(snmfproject, K=3))
  Mk4 <- mean(cross.entropy(snmfproject, K=4))
  SDk4 <- sd(cross.entropy(snmfproject, K=4))
  Mk5 <- mean(cross.entropy(snmfproject, K=5))
  SDk5 <- sd(cross.entropy(snmfproject, K=5))
  Mk6 <- mean(cross.entropy(snmfproject, K=6))
  SDk6 <- sd(cross.entropy(snmfproject, K=6))
  Mk7 <- mean(cross.entropy(snmfproject, K=7))
  SDk7 <- sd(cross.entropy(snmfproject, K=7))
  Mk8 <- mean(cross.entropy(snmfproject, K=8))
  SDk8 <- sd(cross.entropy(snmfproject, K=8))
  Mk9 <- mean(cross.entropy(snmfproject, K=9))
  SDk9 <- sd(cross.entropy(snmfproject, K=9))
  Mk10 <- mean(cross.entropy(snmfproject, K=10))
  SDk10 <- sd(cross.entropy(snmfproject, K=10))
  
  CE <- data.frame(K=c(1:10), Mean = c(Mk1,Mk2,Mk3,Mk4,Mk5,Mk6,Mk7,Mk8,Mk9,Mk10),
                   SD = c(SDk1,SDk2,SDk3,SDk4,SDk5,SDk6,SDk7,SDk8,SDk9,SDk10))
  
  library(ggplot2)
  
  ggplot(CE, aes(x=K, y=Mean)) + 
    geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.2)+
    geom_line() + 
    geom_point(size=4, shape=21, fill="red", color="darkred") + xlab("Number of ancestral populations") + ylab("Cross-entropy")+
    theme_bw() +
    theme(axis.text.x=element_text(angle=0, vjust=0.5), axis.text = element_text(size=15, face ="bold" , color = "black"), axis.title.x = element_text(size=15, face="bold", color="black"),axis.title.y = element_text(size=15, face="bold", color="black")) +
    scale_x_continuous(breaks = seq(0,10, by=2))
}

Best.run <- function(nrep, optimalK, p1, p2, p3, p4){
  ce1 = LEA::cross.entropy(p1, K = optimalK) # get the cross-entropy of each run for optimal K
  ce2 = LEA::cross.entropy(p2, K = optimalK) # get the cross-entropy of each run for optimal K
  ce3 = LEA::cross.entropy(p3, K = optimalK) # get the cross-entropy of each run for optimal K
  ce4 = LEA::cross.entropy(p4, K = optimalK) # get the cross-entropy of each run for optimal K
  AllProjects <- rbind(ce1, ce2, ce3, ce4)
  rownames(AllProjects) <- NULL
  AllProjects <- as.data.frame(AllProjects)
  AllProjects$Project <- c(rep(1, nrep), rep(2,nrep), rep(3, nrep), rep(4, nrep))
  Best_project <- AllProjects[AllProjects[, 1]==min(AllProjects[, 1]), 2]
  Best_runs <- AllProjects[AllProjects$Project==Best_project, ]
  Best_runs$Nrun <- 1:nrow(Best_runs)
  Best_run <- Best_runs[Best_runs[, 1]==min(Best_runs[, 1]), 3]
  print(paste0("Best run is: ", "project = ", Best_project, ", run = ", Best_run)) 
}

barplotK <- function(Qfile, Pop, Run_B){
  Q.matrix_1 <-LEA::Q(Qfile, K = Pop, run = Run_B)
  Q.matrix <- as.qmatrix(Q.matrix_1)
  barplot(Q.matrix, xlab = "Sampled individuals",
          ylab = "Ancestry coefficients",
          main = "", cex.axis = 1.5, cex.lab = 1)
}

fst = function(project, run = 1, K, ploidy = 2){
  library(LEA)
  l = dim(G(project, K = K, run = run))[1]
  q = apply(Q(project, K = K, run = run), MARGIN = 2, mean)
  if (ploidy == 2) {
    G1.t = G(project, K = K, run = run)[seq(2,l,by = 3),]
    G2.t = G(project, K = K, run = run)[seq(3,l,by = 3),]
    freq = G1.t/2 + G2.t}
  else {
    freq = G(project, K = K, run = run)[seq(2,l,by = 2),]}
  H.s = apply(freq*(1-freq), MARGIN = 1, FUN = function(x) sum(q*x))
  P.t = apply(freq, MARGIN = 1, FUN = function(x) sum(q*x))
  H.t = P.t*(1-P.t)
  return(1-H.s/H.t)
}

GIF <- function(project, run, K, fst.values){
  n = dim(Q(project, K, run))[1]
  fst.values[fst.values<0] = 0.000001
  z.scores = sqrt(fst.values*(n-K)/(1-fst.values))
  lambda = median(z.scores^2)/qchisq(1/2, df = K-1)
  print(lambda)
  return(lambda)
}

candidates <- function(alpha=0.05, adj.p.values){ ## alpha is significance level
  L = length(adj.p.values) ## Number of loci
  w = which(sort(adj.p.values) < alpha * (1:L) / L)
  candidates = order(adj.p.values)[w]
  print(length(candidates))
  return(candidates)
}

ManPlot <- function(adj.p.values, candidates, title){
  plot(-log10(adj.p.values), main=title, xlab = "Locus", cex = .7, col = "grey")
  points(candidates, -log10(adj.p.values)[candidates], pch = 19, cex = .7, col = "red")
}


#Load VCF
vcf <- read.vcfR("v.vcf") #new SNP with more SNPs retained
vcf <- filter_biallelic(vcf) #1156 SNPs, 0.027% of all input SNPs contained more than 2 alleles 

#convert VCF into a genlight
vcf.gl <- vcfR2genlight(vcf) 
pop(vcf.gl) <- as.factor(c("PGR", "PGR", "BHE", "BHE", "BHE", "BHE", "BHE", "CHC", "CHC",
                           "CHC", "CHC", "CHC", "DUF", "DUF",
                           "DUF", "DUR", "DUR", "DUR",
                           "DUR", "GIS", "GIS", "GIS", "GMN", "GMN", "GMN",
                           "GMN", "IVO", "KNA", "KNA", "KNA", "KNA",
                           "PMN", "PGR", "PGR", "PKW", "PKW",
                           "PKW", "PKW", "PKW", "TAO", "TAO", "TAO",
                           "TAO", "TAO", "TGA", "TGA", "TGA", "TGA", "TGA", "TEA", "TEA", "WLG",
                           "WLG", "WLG", "WLG", "WLG", "WTP", "WTP", "WTP"))

#convert genlight to genind
vcf.gi <- gl2gi(vcf.gl)

#PCA 
PCA <- tab(vcf.gi, freq=TRUE, NA.method="mean")
PCA.data <- dudi.pca(PCA, scannf= FALSE, center=TRUE, scale=FALSE) 
PCA.data

percent <- PCA.data$eig/sum(PCA.data$eig)*100 
percent

ind.coords <- as.data.frame(PCA.data$li)
colnames(ind.coords) <- c("Axis1", "Axis2")

ind.coords$Ind <- indNames(vcf.gi) 
ind.coords$pop <- vcf.gi$pop 

centroid = aggregate(cbind(Axis1, Axis2) ~ pop, data = ind.coords, FUN = mean)
ind.coords <- left_join(ind.coords, centroid, by = "pop", suffix = c("", ".cen"))

#custom x and y labels 
xlab = paste("PC 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("PC 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

#Create colour pallette 
cc <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
        "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
        "#c49c94", "#f7b6d2")

ggtheme= theme(axis.text.y = element_text(colour="black", size=15),
               axis.text.x = element_text(colour="black", size=15),
               axis.title = element_text(colour="black", size=15),
               panel.border = element_rect(colour="black", fill=NA, size=1),
               panel.background = element_blank(),
               plot.title = element_text(hjust=0.5, size=15),
               legend.text = element_text(size = 14), 
               legend.title = element_text(size = 15, face = 'bold'))



ggplot(data = ind.coords, aes(x = Axis1, y = Axis2), frame=TRUE) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_point(aes(fill = pop), shape = 21, size = 5, show.legend = T)+
  scale_fill_manual(values = cc)+
  scale_colour_manual(values = cc)+
  labs(x = xlab, y = ylab)+
  labs(fill='Population')+
  ggtheme


##PCA for only North Island, South Island, and Australia

pop(vcf.gi) <- as.factor(c("North Island, NZ", "North Island, NZ", "South Island, NZ", "South Island, NZ", "South Island, NZ", "South Island, NZ",
                           "South Island, NZ", "South Island, NZ", "South Island, NZ",
                           "South Island, NZ", "South Island, NZ", "South Island, NZ", "South Island, NZ", "South Island, NZ",
                           "South Island, NZ", "South Island, NZ", "South Island, NZ", "South Island, NZ",
                           "South Island, NZ", "North Island, NZ", "North Island, NZ", "North Island, NZ", "South Island, NZ", "South Island, NZ", "South Island, NZ",
                           "South Island, NZ", "South Island, NZ", "Australia", "Australia", "Australia", "Australia",
                           "North Island, NZ", "North Island, NZ", "North Island, NZ", "North Island, NZ", "North Island, NZ",
                           "North Island, NZ", "North Island, NZ", "North Island, NZ", "North Island, NZ", "North Island, NZ", "North Island, NZ",
                           "North Island, NZ", "North Island, NZ", "North Island, NZ", "North Island, NZ", "North Island, NZ", "North Island, NZ",
                           "North Island, NZ", "North Island, NZ", "North Island, NZ", "North Island, NZ",
                           "North Island, NZ", "North Island, NZ", "North Island, NZ", "North Island, NZ", "North Island, NZ", "North Island, NZ", "North Island, NZ"))

#PCA 
PCA <- tab(vcf.gi, freq=TRUE, NA.method="mean")
PCA.data <- dudi.pca(PCA, scannf= FALSE, center=TRUE, scale=FALSE) 
PCA.data

percent <- PCA.data$eig/sum(PCA.data$eig)*100 
percent

ind.coords <- as.data.frame(PCA.data$li)
colnames(ind.coords) <- c("Axis1", "Axis2")

ind.coords$Ind <- indNames(vcf.gi) 
ind.coords$pop <- vcf.gi$pop 

centroid = aggregate(cbind(Axis1, Axis2) ~ pop, data = ind.coords, FUN = mean)
ind.coords <- left_join(ind.coords, centroid, by = "pop", suffix = c("", ".cen"))

#custom x and y labels 
xlab = paste("PC 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("PC 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

#Create colour pallette 
cc <- c("#1f77b4", "#ff7f0e", "#2ca02c")

ggtheme= theme(axis.text.y = element_text(colour="black", size=15),
               axis.text.x = element_text(colour="black", size=15),
               axis.title = element_text(colour="black", size=15),
               panel.border = element_rect(colour="black", fill=NA, size=1),
               panel.background = element_blank(),
               plot.title = element_text(hjust=0.5, size=15),
               legend.text = element_text(size = 14), 
               legend.title = element_text(size = 15, face = 'bold'))



ggplot(data = ind.coords, aes(x = Axis1, y = Axis2), frame=TRUE) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_point(aes(fill = pop), shape = 21, size = 5, show.legend = T)+
  scale_fill_manual(values = cc)+
  scale_colour_manual(values = cc)+
  labs(x = xlab, y = ylab)+
  labs(fill='Population')+
  ggtheme


#----Original Admixtures (Population Assignment using sNMF to estimate K) ------
snps_fil_ldF <-  vcfLink("vicina.vcf")


#make the .geno file 
gl2geno(vcf.gl, outfile = "gl_geno", outpath = getwd(), verbose = NULL)

# Run SNMF (LEA) using different alpha
set.seed(1234)
project_snmf1 = snmf("gl_geno.geno", K=1:10, rep =10, entropy=T, project = "new", CPU=4) #CPU=4 uses 4 CPUs 

#load snmf project
project1 = load.snmfProject("gl_geno.snmfProject")

##Cross-Entropy plot
par(mar=c(5,6,4,1)+.1)
plot(project1, col = "red", pch=19)

optimal_K = 5

#extract cross entropy of all runs where K=2
ce = cross.entropy(project1, K = optimal_K)
ce

# Find the run with the lowest cross-entropy
lowest.ce = which.min(ce)
lowest.ce

#Making the Q-matrix
qmatrix = as.data.frame(Q(project1, K = optimal_K, run = lowest.ce)) 
head(qmatrix)

#Label column names of qmatrix
ncol(qmatrix)
cluster_names = c()
for (i in 1:ncol(qmatrix)){
  cluster_names[i] = paste("Cluster", i)
}
cluster_names
colnames(qmatrix) = cluster_names
head(qmatrix)

#load in the meta stuff
popIds <- read_excel("meta.xlsx")

# Add individual IDs
Ind <- popIds$`SampleID`
qmatrix$Ind = Ind

# Add site IDs
Site <- popIds$Location
qmatrix$Site = Site

# Convert dataframe to long format
qlong = melt(qmatrix, id.vars=c("Ind","Site"))
head(qlong)

# Adjust facet labels (if you want to change the order). 
levels(qlong$Site)
facet.labs = c("TEA", "PGR", "PKW", "TGA", "GIS",
               "TAO", "PMN", "WLG", "WTP", "BHE",
               "GMN", "CHC", "DUF", "DUR", "IVO", "KNA")
levels(qlong$Site) = facet.labs
levels(qlong$Site)

# Define colour palette
pal = colorRampPalette(c("cyan","green1", "deeppink", "yellow","blue1")) 
cols = pal(length(unique(qlong$variable)))

#order levels before plotting - you can change the order of these and it will change the order of how they are plotted in the bar graph :)
qlong$Site <- ordered(qlong$Site, levels = c("TEA", "PGR", "PKW", "TGA", "GIS",
                                             "TAO", "PMN", "WLG", "WTP", "BHE",
                                             "GMN", "CHC", "DUF", "DUR", "IVO", "KNA")) 
qlong$Site

#plot  
ggplot(data=qlong, aes(x=Ind, y=value, fill = variable)) +
  geom_bar(stat="identity")+
  geom_col(color = "gray", size = 0.1) +
  scale_fill_manual(values = cols)+
  theme_minimal() + labs(x = "Individuals", title = "K=5", y = "Admixture proportion") +
  facet_grid(~Site, scales="free", space="free")+
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    panel.spacing.x = unit(0, "lines"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_text(colour="black", size=9),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    legend.position = "none"
  )

#---------------------- Making neutral VCF ---------------------------------

snps_fil_ldF@meta$PopID_snmf <- popIds

#Save new vcf file and pop ID file
Save(snps_fil_ldF, "NEUTRAL.vcf")   #produces a 'meta file'.

#These next steps are all to work out the neutral SNPs only:

#Compute the FST statistics using best run
best1 = which.min(ce)
FST = fst(project1, best1, optimal_K) #you need at least 2 populations for a population-based test, so K>1.

#Compute the GIF
lambda <- GIF(project1, best1, optimal_K, fst.values=FST) 

#Compute adjusted p-values from the combined z-scores and plot histogram of p-values
n = dim(Q(project1, best1, optimal_K))[1]
z.scores = sqrt(FST*(n-optimal_K)/(1-FST)) #NaNs produced 
adj.p.values = pchisq(z.scores^2/lambda, df = optimal_K-1, lower = FALSE)
hist(adj.p.values, col = "red")

#Test different lambda values and plot histogram of p-values
adj.p.values = pchisq(z.scores^2/1.5, df = optimal_K-1, lower = FALSE)
hist(adj.p.values, col = "green")

#Candidate loci
#FDR control: Benjamini-Hochberg at level q
C_fst <- candidates(alpha=0.05, adj.p.values)

#Manhatan plot
ManPlot(adj.p.values, C_fst,"Fst")

#Exclude candidate FST outlier
snps_fil_ldF_candidate <- Subset(snps_fil_ldF, sites=C_fst)
snps_fil_ldF_candidate@site_id ## These are all the candidate SNPs
Chrom(snps_fil_ldF_candidate)

candidates_fst <- snps_fil_ldF_candidate@site_id
All_snp <- snps_fil_ldF@site_id
N_snp <- All_snp[!(All_snp %in% candidates_fst)] #Exclude all candidate loci

snps_fil_ldF_neutral <- Subset(snps_fil_ldF, sites=N_snp)
VCFsummary(snps_fil_ldF_neutral)

length(N_snp)
length(snps_fil_ldF@site_id)-length(C_fst) 
length(snps_fil_ldF_neutral@site_id)

#Save neutral snp dataset
Save(snps_fil_ldF_neutral, "NEUTRAL_SNPS.vcf")
neutral.vcf <- read.vcfR("NEUTRAL_SNPS.vcf") 

#COMPARE VCF's
neutral.vcf ##the neutral vcf has 21,031 variants
vcf #the normal vcf has 21,354 variants 

#323 snps were removed

#####THE FOLLOWING ANALYSIS IS JUST A REPETITION OF THE PREVIOUS WITH THE NEUTRAL DATASET

vcf <- read.vcfR("NEUTRAL_new.vcf") 
vcf <- filter_biallelic(vcf)

#######PCA WITH THE NEUTRAL DATASET########

#convert VCF into a genlight 
vcf.gl <- vcfR2genlight(vcf) 
pop(vcf.gl) <- as.factor(c("PGR", "PGR", "BHE", "BHE", "BHE", "BHE", "BHE", "CHC", "CHC",
                           "CHC", "CHC", "CHC", "DUF", "DUF",
                           "DUF", "DUR", "DUR", "DUR",
                           "DUR", "GIS", "GIS", "GIS", "GMN", "GMN", "GMN",
                           "GMN", "IVO", "KNA", "KNA", "KNA", "KNA",
                           "PMN", "PGR", "PGR", "PKW", "PKW",
                           "PKW", "PKW", "PKW", "TAO", "TAO", "TAO",
                           "TAO", "TAO", "TGA", "TGA", "TGA", "TGA", "TGA", "TEA", "TEA", "WLG",
                           "WLG", "WLG", "WLG", "WLG", "WTP", "WTP", "WTP"))

#convert genlight to genind
vcf.gi <- gl2gi(vcf.gl)
#PCA 
PCA <- tab(vcf.gi, freq=TRUE, NA.method="mean")
PCA.data <- dudi.pca(PCA, scannf= FALSE, center=TRUE, scale=FALSE) 
PCA.data

percent <- PCA.data$eig/sum(PCA.data$eig)*100 
percent

ind.coords <- as.data.frame(PCA.data$li)
colnames(ind.coords) <- c("Axis1", "Axis2")

ind.coords$Ind <- indNames(vcf.gi) 
ind.coords$pop <- vcf.gi$pop 

centroid = aggregate(cbind(Axis1, Axis2) ~ pop, data = ind.coords, FUN = mean)
ind.coords <- left_join(ind.coords, centroid, by = "pop", suffix = c("", ".cen"))

#custom x and y labels 
xlab = paste("PC 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("PC 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

cc <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
        "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
        "#c49c94", "#f7b6d2")

ggtheme= theme(axis.text.y = element_text(colour="black", size=15),
               axis.text.x = element_text(colour="black", size=15),
               axis.title = element_text(colour="black", size=15),
               panel.border = element_rect(colour="black", fill=NA, size=1),
               panel.background = element_blank(),
               plot.title = element_text(hjust=0.5, size=15),
               legend.text = element_text(size = 14), 
               legend.title = element_text(size = 15, face = 'bold'))



ggplot(data = ind.coords, aes(x = Axis1, y = Axis2), frame=TRUE) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_point(aes(fill = pop), shape = 21, size = 5, show.legend = T)+
  scale_fill_manual(values = cc)+
  scale_colour_manual(values = cc)+
  labs(x = xlab, y = ylab)+
  labs(fill='Population')+
  ggtheme

##PCA for only North Island, South Island, and Australia

pop(vcf.gl) <- as.factor(c("North Island, NZ", "North Island, NZ", "South Island, NZ", "South Island, NZ", "South Island, NZ", "South Island, NZ",
                           "South Island, NZ", "South Island, NZ", "South Island, NZ",
                           "South Island, NZ", "South Island, NZ", "South Island, NZ", "South Island, NZ", "South Island, NZ",
                           "South Island, NZ", "South Island, NZ", "South Island, NZ", "South Island, NZ",
                           "South Island, NZ", "North Island, NZ", "North Island, NZ", "North Island, NZ", "South Island, NZ", "South Island, NZ", "South Island, NZ",
                           "South Island, NZ", "South Island, NZ", "Australia", "Australia", "Australia", "Australia",
                           "North Island, NZ", "North Island, NZ", "North Island, NZ", "North Island, NZ", "North Island, NZ",
                           "North Island, NZ", "North Island, NZ", "North Island, NZ", "North Island, NZ", "North Island, NZ", "North Island, NZ",
                           "North Island, NZ", "North Island, NZ", "North Island, NZ", "North Island, NZ", "North Island, NZ", "North Island, NZ",
                           "North Island, NZ", "North Island, NZ", "North Island, NZ", "North Island, NZ",
                           "North Island, NZ", "North Island, NZ", "North Island, NZ", "North Island, NZ", "North Island, NZ", "North Island, NZ", "North Island, NZ"))

#convert genlight to genind
vcf.gi <- gl2gi(vcf.gl)


#PCA 
PCA <- tab(vcf.gi, freq=TRUE, NA.method="mean")
PCA.data <- dudi.pca(PCA, scannf= FALSE, center=TRUE, scale=FALSE) 
PCA.data

percent <- PCA.data$eig/sum(PCA.data$eig)*100 
percent

ind.coords <- as.data.frame(PCA.data$li)
colnames(ind.coords) <- c("Axis1", "Axis2")

ind.coords$Ind <- indNames(vcf.gi) 
ind.coords$pop <- vcf.gi$pop 

centroid = aggregate(cbind(Axis1, Axis2) ~ pop, data = ind.coords, FUN = mean)
ind.coords <- left_join(ind.coords, centroid, by = "pop", suffix = c("", ".cen"))

#custom x and y labels 
xlab = paste("PC 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("PC 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

#Create colour pallette 
cc <- c("#1f77b4", "#ff7f0e", "#2ca02c")

ggtheme= theme(axis.text.y = element_text(colour="black", size=15),
               axis.text.x = element_text(colour="black", size=15),
               axis.title = element_text(colour="black", size=15),
               panel.border = element_rect(colour="black", fill=NA, size=1),
               panel.background = element_blank(),
               plot.title = element_text(hjust=0.5, size=15),
               legend.text = element_text(size = 14), 
               legend.title = element_text(size = 15, face = 'bold'))



ggplot(data = ind.coords, aes(x = Axis1, y = Axis2), frame=TRUE) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_point(aes(fill = pop), shape = 21, size = 5, show.legend = T)+
  scale_fill_manual(values = cc)+
  scale_colour_manual(values = cc)+
  labs(x = xlab, y = ylab)+
  labs(fill='Population')+
  ggtheme



#----NEUTRAL DATASET Admixtures (Population Assignment using sNMF to estimate K) ------

snps_fil_ldF <-  vcfLink("NEUTRAL_SNPS.vcf")

#make the .geno file 
gl2geno(vcf.gl, outfile = "gl_geno", outpath = getwd(), verbose = NULL)

# Run SNMF (LEA) using different alpha
set.seed(1234)
project_snmf1 = snmf("gl_geno.geno", K=1:10, rep =10, entropy=T, project = "new", CPU=4) #CPU=4 uses 4 CPUs 

#load snmf project
project1 = load.snmfProject("gl_geno.snmfProject")

##Cross-Entropy plot
par(mar=c(5,6,4,1)+.1)
plot(project1, col = "red", pch=19)

optimal_K = 5

#extract cross entropy of all runs where K=2
ce = cross.entropy(project1, K = optimal_K)
ce

# Find the run with the lowest cross-entropy
lowest.ce = which.min(ce)
lowest.ce

#Making the Q-matrix
qmatrix = as.data.frame(Q(project1, K = optimal_K, run = lowest.ce)) 
head(qmatrix)

#Label column names of qmatrix
ncol(qmatrix)
cluster_names = c()
for (i in 1:ncol(qmatrix)){
  cluster_names[i] = paste("Cluster", i)
}
cluster_names
colnames(qmatrix) = cluster_names
head(qmatrix)

#load in the meta stuff
popIds <- read_excel("meta.xlsx")

# Add individual IDs
Ind <- popIds$`SampleID`
qmatrix$Ind = Ind

# Add site IDs
Site <- popIds$Location
qmatrix$Site = Site

# Convert dataframe to long format
qlong = melt(qmatrix, id.vars=c("Ind","Site"))
head(qlong)

# Adjust facet labels (if you want to change the order)
levels(qlong$Site)
facet.labs = c("TEA", "PGR", "PKW", "TGA", "GIS",
               "TAO", "PMN", "WLG", "WTP", "BHE",
               "GMN", "CHC", "DUF", "DUR", "IVO", "KNA")
levels(qlong$Site) = facet.labs
levels(qlong$Site)

# Define colour palette
pal = colorRampPalette(c("cyan","green1", "deeppink", "yellow","blue1")) 
cols = pal(length(unique(qlong$variable)))

#order levels before plotting - you can change the order of these and it will change the order of how they are plotted in the bar graph :)
qlong$Site <- ordered(qlong$Site, levels = c("TEA", "PGR", "PKW", "TGA", "GIS",
                                             "TAO", "PMN", "WLG", "WTP", "BHE",
                                             "GMN", "CHC", "DUF", "DUR", "IVO", "KNA")) 
qlong$Site

#plot  
ggplot(data=qlong, aes(x=Ind, y=value, fill = variable)) +
  geom_bar(stat="identity")+
  geom_col(color = "gray", size = 0.1) +
  scale_fill_manual(values = cols)+
  theme_minimal() + labs(x = "Individuals", title = "K=5", y = "Admixture proportion") +
  facet_grid(~Site, scales="free", space="free")+
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    panel.spacing.x = unit(0, "lines"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_text(colour="black", size=9),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    legend.position = "none"
  )

#######BASIC STATS, FSTS, Ho, He, ETC., all with neutral vcf##########

vcf.gl <- gl.drop.ind(vcf.gl, "K76173", recalc = TRUE, mono.rm = FALSE, verbose = NULL) #remove PMN pop due to only having 1 individual (heterozygosity calc won't work otherwise)
vcf.gl <- gl.drop.ind(vcf.gl, "K76126", recalc = TRUE, mono.rm = FALSE, verbose = NULL) #remove IVO pop due to only having 1 individual
vcf.gl <- gl.drop.ind(vcf.gl, "K76245", recalc = TRUE, mono.rm = FALSE, verbose = NULL) #remove TEA pop due to only having 2 individuals
vcf.gl <- gl.drop.ind(vcf.gl, "K76246", recalc = TRUE, mono.rm = FALSE, verbose = NULL) #remove TEA pop due to only having 2 individuals

vcf.gl #check there is four less individuals (i.e., n=55)

#Redo the populations removing PMN, IVO, and TEA
pop(vcf.gl) <- as.factor(c("PGR", "PGR", "BHE", "BHE", "BHE", "BHE", "BHE", "CHC", "CHC",
                           "CHC", "CHC", "CHC", "DUF", "DUF",
                           "DUF", "DUR", "DUR", "DUR",
                           "DUR", "GIS", "GIS", "GIS", "GMN", "GMN", "GMN",
                           "GMN", "KNA", "KNA", "KNA", "KNA",
                           "PGR", "PGR", "PKW", "PKW",
                           "PKW", "PKW", "PKW", "TAO", "TAO", "TAO",
                           "TAO", "TAO", "TGA", "TGA", "TGA", "TGA", "TGA", "WLG",
                           "WLG", "WLG", "WLG", "WLG", "WTP", "WTP", "WTP"))

#convert genlight to genind
vcf.gi <- gl2gi(vcf.gl)

#calculate basic stats using hierfstat (using genind object)
basic_data = basic.stats(vcf.gi, diploid = TRUE)
basic_data

#Calculate mean observed heterozygosity per site 
Ho <- apply(basic_data$Ho, MARGIN = 2, FUN = mean, na.rm = TRUE) %>% 
  round(digits = 3)
Ho

#mean expected heterozygosity per site
He = apply(basic_data$Hs, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 3)
He

##Can plot heterozygosity###
het = data.frame("pop" = names(Ho), Ho = Ho, He = He) %>%
  melt(id.vars = "pop")
custom.theme = theme(
  axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, face = "bold"),
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size = 12),
  axis.title.x = element_blank(),
  axis.line.y = element_line(size = 0.5),
  legend.title = element_blank(),
  legend.text = element_text(size = 12),
  panel.grid = element_blank(),
  panel.background = element_blank(),
  plot.title = element_text(hjust = 0.5, size = 15, face="bold")
)
hetlab.o = expression(italic("H")[o])
hetlab.e = expression(italic("H")[e])

ggplot(data = het, aes(x = pop, y = value, fill = variable))+
  geom_bar(stat = "identity", position = position_dodge(width = 0.6), colour = "black")+
  scale_y_continuous(expand = c(0,0), limits = c(0,0.30))+
  scale_fill_manual(values = c("royalblue", "#bdbdbd"), labels = c(hetlab.o, hetlab.e))+
  ylab("Heterozygosity")+
  ggtitle("C. vicina")+
  custom.theme

#Calculate mean FIS per site.
FIS = apply(basic_data$Fis, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 3)
FIS

##Pairwise Fst using heirfstat##
#subset data to reduce computation time 
gen_sub <- popsub(vcf.gi, c("PKW","PGR", "TGA", "GIS", "TAO",
                            "WLG", "WTP", "BHE", "GMN", "CHC",
                            "CHC", "DUF", "DUR", "KNA"))
pw_fst = genet.dist(gen_sub, method = "WC84") %>% round(digits = 3) #this takes like 3 hours to run
pw_fst

#calculate pw bootstrapping values
boot.ppfst(dat=gen_sub,nboot=100,quant=c(0.025,0.975),diploid=TRUE)


#visualise pairwise Fsts from heirfstat
#Determine order of labels 
lab_order <- c("PKW","PGR", "TGA", "GIS", "TAO",
               "WLG", "WTP", "BHE", "GMN", "CHC",
               "CHC", "DUF", "DUR", "KNA")

#Change order of rows and cols
fst_mat <- as.matrix(pw_fst)
fst_mat1 <- fst_mat[lab_order, ]
fst_mat2 <- fst_mat1[, lab_order]

#Create a data frame
ind <- which(upper.tri(fst_mat2), arr.ind=TRUE)
fst.df <- data.frame(Site1=dimnames(fst_mat2) [[2]] [ind[,2]], Site2=dimnames(fst_mat2)[[1]][ind[,1]], Fst= fst_mat2[ ind ])

#Keep the order of the levels in the data frame for plotting
fst.df$Site1 = factor(fst.df$Site1, levels = unique(fst.df$Site1))
fst.df$Site2 = factor(fst.df$Site2, levels = unique(fst.df$Site2))

#Convert minus values to zero
fst.df$Fst[fst.df$Fst < 0] = 0

#Print data frame summary
fst.df %>% str

#print Fst italic label
fst.label = expression(italic("F")[ST])

#Extract middle Fst value for gradient argument
mid=max(fst.df$Fst)/2

library(scales)

#Plot heatmap 
ggplot(data = fst.df, aes(x = Site1, y = Site2, fill = Fst))+
  geom_tile(colour = "black")+
  geom_text(aes(label = Fst), color="black", size = 4.5)+
  scale_fill_gradient2(low = "lightgoldenrod1", mid = "tomato1", high = "orangered3", midpoint = mid, name = fst.label, limits = c(0, max(fst.df$Fst)), breaks = c(0.25, 0.5, 0.75, 1))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), position = "right")+
  theme(axis.text.y = element_text(colour = "black", size = 12), axis.text.x= element_text(colour="black", size=12),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 10)
  )+scale_x_discrete(labels = wrap_format(10))


#different way to do Fst calculations with p-value (using dartR package)
#new FST

#Load VCF
vcf <- read.vcfR("vicina.vcf") 
vcf <- filter_biallelic(vcf)
#convert VCF into a genlight
vcf.gl <- vcfR2genlight(vcf) 

vcf.gl <- gl.drop.ind(vcf.gl, "K76173", recalc = TRUE, mono.rm = FALSE, verbose = NULL) #remove PMN pop due to only having 1 individual (heterozygosity calc won't work otherwise)
vcf.gl <- gl.drop.ind(vcf.gl, "K76126", recalc = TRUE, mono.rm = FALSE, verbose = NULL) #remove IVO pop due to only having 1 individual
vcf.gl <- gl.drop.ind(vcf.gl, "K76245", recalc = TRUE, mono.rm = FALSE, verbose = NULL) #remove TEA pop due to only having 2 individuals
vcf.gl <- gl.drop.ind(vcf.gl, "K76246", recalc = TRUE, mono.rm = FALSE, verbose = NULL) #remove TEA pop due to only having 2 individuals

vcf.gl #check there is four less individuals (i.e., n=55)

#Redo the populations removing PMN, IVO, and TEA
pop(vcf.gl) <- as.factor(c("PGR", "PGR", "BHE", "BHE", "BHE", "BHE", "BHE", "CHC", "CHC",
                           "CHC", "CHC", "CHC", "DUF", "DUF",
                           "DUF", "DUR", "DUR", "DUR",
                           "DUR", "GIS", "GIS", "GIS", "GMN", "GMN", "GMN",
                           "GMN", "KNA", "KNA", "KNA", "KNA",
                           "PGR", "PGR", "PKW", "PKW",
                           "PKW", "PKW", "PKW", "TAO", "TAO", "TAO",
                           "TAO", "TAO", "TGA", "TGA", "TGA", "TGA", "TGA", "WLG",
                           "WLG", "WLG", "WLG", "WLG", "WTP", "WTP", "WTP"))


#do pairwise FST and calculate p values 
fst <- gl.fst.pop(vcf.gl, nboots = 1000, percent = 95, nclusters = 1, verbose = NULL)
fst

#heatmap

lab_order <- c("PKW","PGR", "TGA", "GIS", "TAO",
               "WLG", "WTP", "BHE", "GMN", "CHC",
               "CHC", "DUF", "DUR", "KNA")

#Change order of rows and cols
fst_mat <- as.matrix(fst$Fsts)
fst_mat
#Create a data frame
ind <- which(lower.tri(fst_mat), arr.ind=TRUE)
fst.df <- data.frame(Site1=dimnames(fst_mat) [[2]] [ind[,2]], Site2=dimnames(fst_mat)[[1]][ind[,1]], Fst= fst_mat[ ind ])


#Keep the order of the levels in the data frame for plotting
fst.df$Site1 = factor(fst.df$Site1, levels = unique(fst.df$Site1))
fst.df$Site2 = factor(fst.df$Site2, levels = unique(fst.df$Site2))

#Convert minus values to zero
fst.df$Fst[fst.df$Fst < 0] = 0

#Print data frame summary
fst.df %>% str

#print Fst italic label
fst.label = expression(italic("F")[ST])

#Extract middle Fst value for gradient argument
mid=max(fst.df$Fst)/2


#Plot heatmap 
ggplot(data = fst.df, aes(x = Site1, y = Site2, fill = Fst))+
  geom_tile(colour = "black")+
  geom_text(aes(label = scales::number(Fst, accuracy = 0.001)), color="black", size = 4.5)+
  scale_fill_gradient2(low = "lightgoldenrod1", mid = "tomato1", high = "orangered3", midpoint = mid, name = fst.label, limits = c(0, max(fst.df$Fst)), breaks = c(0.25, 0.5, 0.75, 1))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), position = "right")+
  theme(axis.text.y = element_text(colour = "black", size = 12), axis.text.x= element_text(colour="black", size=12),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 10)
  )+scale_x_discrete(labels = wrap_format(10))






