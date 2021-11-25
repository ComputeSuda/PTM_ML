#Load package
library(bio3d)
library(readxl)
library(substring)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(stringr)
library(dplyr)

data_only_binding <- read_excel("data.xlsx")  #input dataset
#residue level feature
ki_feature <- data.frame(Uniprot=character(), protein=character(), PDB=character(), res=character(),
                         swiss_resid=character(),functionality=numeric(), anm_top3_bcc=numeric(),
                         gnm_top3_bcc=numeric(), anm_LF_bcc=numeric(), gnm_LF_bcc=numeric(),anm_LTIF_bcc=numeric(),
                         gnm_LTIF_bcc=numeric(), anm_HF_bcc=numeric(),gnm_HF_bcc=numeric(), anm_all_bcc=numeric(),
                         gnm_all_bcc=numeric(), gnm_top3_mode=numeric(), gnm_slow_mode=numeric(), gnm_all_mode=numeric(),
                         b_factor=numeric(), betweenness=numeric(), closeness=numeric(), degree=numeric(),
                         cluster=numeric(), diversity=numeric(), eccentricity=numeric(), strength=numeric(),
                         page_rank=numeric(), mean_sp_bs=numeric(), anm_stiffness=numeric(), ACC=numeric(),
                         Phi=numeric(), PSi=numeric(), anm_effector=numeric(), anm_sensor=numeric(),
                         gnm_effector=numeric(), gnm_sensor=numeric(), hitting_time_col=numeric(), commute_time=numeric(),
                         Iij_betweenness=numeric(), Iij_closeness=numeric(), Iij_degree=numeric(), Iij_cluster=numeric(),
                         Iij_diversity=numeric(), Iij_eccentricity=numeric(), Iij_strength=numeric(),
                         Iij_page_rank=numeric(), hitting_time_row=numeric(), anm_sq=numeric(), gnm_sq=numeric(),
                         stringsAsFactors=FALSE)
for (i in 1:84) {
  setwd("D:/subject/dataset/PSP/regulatory_site/kinase_renumber_pdb")
  pdb <- read.pdb2(paste(as.character(data_only_binding[i,3]), "pdb", sep = "."))  #read pdb file
  pdb_sequence <- as.data.frame(pdbseq(pdb))
  pdb_sequence[, 1] <- as.character(pdb_sequence[, 1])
  colnames(pdb_sequence) <- "res"
  swiss_start <- as.numeric(strsplit(as.character(data_only_binding[i, 13]), split = "-")[[1]][1])
  swiss_end <- as.numeric(strsplit(as.character(data_only_binding[i, 13]), split = "-")[[1]][2])
  pdb_sequence$swiss_resid <- swiss_start:swiss_end  #collect pdb information
  pdb_sequence$Uniprot <- as.character(data_only_binding[i,1])
  pdb_sequence$protein <- as.character(data_only_binding[i,2])
  pdb_sequence$PDB <- as.character(data_only_binding[i,3])
  pdb_sequence$functionality <- "-"
  function_site <- as.numeric(strsplit(as.character(data_only_binding[i, 14]), split = ",")[[1]])
  bind_site <- as.numeric(strsplit(as.character(data_only_binding[i, 16]), split = ",")[[1]])
  non_fun_ptm <- as.numeric(strsplit(as.character(data_only_binding[i, 15]), split = ",")[[1]])
  non_ptm <- as.numeric(strsplit(as.character(data_only_binding[i, 20]), split = ",")[[1]])
  rownames(pdb_sequence) <- NULL
  pdb_sequence[non_ptm - swiss_start + 1, 6] <- "NF"
  
  if(!is.na(bind_site)){
    pdb_sequence[bind_site - swiss_start + 1, 6] <- "bind"
  }
  if(!is.na(function_site)){
    pdb_sequence[function_site - swiss_start + 1, 6] <- "PTM"
  }
  if(!is.na(non_fun_ptm)){
    pdb_sequence[non_fun_ptm - swiss_start + 1, 6] <- "PTM"
  }
  
  ptm <- c(function_site, non_fun_ptm)
  P_B <- intersect(ptm, bind_site)
  if(length(P_B) > 0){
    pdb_sequence[P_B - swiss_start + 1, 6] <- "P/B"
  }
  
  pdb_sequence <- pdb_sequence[ ,c(3:5,1,2,6)]
  df <- data.frame(anm_top3_bcc=NA,gnm_top3_bcc=NA,anm_LF_bcc=NA,gnm_LF_bcc=NA,anm_LTIF_bcc=NA,
                   gnm_LTIF_bcc=NA,anm_HF_bcc=NA,gnm_HF_bcc=NA,anm_all_bcc=NA,gnm_all_bcc=NA,
                   gnm_top3_mode=NA,gnm_slow_mode=NA,gnm_all_mode=NA,b_factor=NA,betweenness=NA,
                   closeness=NA,degree=NA,cluster=NA,diversity=NA,eccentricity=NA,strength=NA,
                   page_rank=NA,mean_sp_bs=NA,anm_stiffness=NA,ACC=NA,Phi=NA,PSi=NA,anm_effector=NA,
                   anm_sensor=NA,gnm_effector=NA,gnm_sensor=NA,hitting_time_col=NA,commute_time=NA,
                   Iij_betweenness=NA,Iij_closeness=NA,Iij_degree=NA,Iij_cluster=NA,Iij_diversity=NA,
                   Iij_eccentricity=NA,Iij_strength=NA,Iij_page_rank=NA,hitting_time_row=NA,anm_sq=NA,
                   gnm_sq=NA,stringsAsFactors=FALSE)
  pdb_sequence <- cbind(pdb_sequence, df)
  b_s <- bind_site - swiss_start + 1
  non_func_site <- non_ptm - swiss_start + 1
  #collect feature
  if(length(b_s) > 1){
    setwd("D:/subject/dataset/PSP/regulatory_site/PSP_pdb/output_2")
    anm_top3_bcc <- abs(as.matrix(read.table(paste(pdb_sequence[i, 3], "_anm_cc_top3.txt", sep = ""))))
    anm_top3_bcc_index <- rowMeans(anm_top3_bcc[ , b_s])
    
    gnm_top3_bcc <- abs(as.matrix(read.table(paste(pdb_sequence[i, 3], "_gnm_cc_top3.txt", sep = ""))))
    gnm_top3_bcc_index <- rowMeans(gnm_top3_bcc[ , b_s])
    
    anm_LF_bcc <- abs(as.matrix(read.table(paste(pdb_sequence[i, 3], "_anm_cc_LF.txt", sep = ""))))
    anm_LF_bcc_index <- rowMeans(anm_LF_bcc[ , b_s])
    
    gnm_LF_bcc <- abs(as.matrix(read.table(paste(pdb_sequence[i, 3], "_gnm_cc_LF.txt", sep = ""))))
    gnm_LF_bcc_index <- rowMeans(gnm_LF_bcc[ , b_s])
    
    anm_LTIF_bcc <- abs(as.matrix(read.table(paste(pdb_sequence[i, 3], "_anm_cc_LTIF.txt", sep = ""))))
    anm_LTIF_bcc_index <- rowMeans(anm_LTIF_bcc[ , b_s])
    
    gnm_LTIF_bcc <- abs(as.matrix(read.table(paste(pdb_sequence[i, 3], "_gnm_cc_LTIF.txt", sep = ""))))
    gnm_LTIF_bcc_index <- rowMeans(gnm_LTIF_bcc[ , b_s])
    
    anm_HF_bcc <- abs(as.matrix(read.table(paste(pdb_sequence[i, 3], "_anm_cc_HF.txt", sep = ""))))
    anm_HF_bcc_index <- rowMeans(anm_HF_bcc[ , b_s])
    
    gnm_HF_bcc <- abs(as.matrix(read.table(paste(pdb_sequence[i, 3], "_gnm_cc_HF.txt", sep = ""))))
    gnm_HF_bcc_index <- rowMeans(gnm_HF_bcc[ , b_s])
    
    setwd("D:/subject/dataset/PSP/regulatory_site/output_2")
    anm_all_bcc <- abs(as.matrix(read.table(paste(pdb_sequence[i, 3], "_anm_cc_all_mode.txt", sep = ""))))
    anm_all_bcc_index <- rowMeans(anm_all_bcc[ , b_s])
    
    gnm_all_bcc <- abs(as.matrix(read.table(paste(pdb_sequence[i, 3], "_gnm_cc_all_mode.txt", sep = ""))))
    gnm_all_bcc_index <- rowMeans(gnm_all_bcc[ , b_s])
    
    setwd("D:/subject/dataset/PSP/regulatory_site/output_1")
    gnm_top3_mode <- abs(as.matrix(read.table(paste(pdb_sequence[i, 3], "_gnm_mode_1_to_3.txt", sep = ""))))
    gnm_top3_mode_index <- rowMeans(gnm_top3_mode)
    
    gnm_slow_mode <- abs(as.matrix(read.table(paste(pdb_sequence[i, 3], "_gnm_slow_mode.txt", sep = ""))))
    gnm_slow_mode_index <- rowMeans(gnm_slow_mode)
    
    gnm_all_mode <- abs(as.matrix(read.table(paste(pdb_sequence[i, 3], "_gnm_all_mode.txt", sep = ""))))
    gnm_all_mode_index <- rowMeans(gnm_all_mode)
    
    anm_stiffness <- abs(as.matrix(read.table(paste(pdb_sequence[i, 3], "_anm_stiffness.txt", sep = ""))))
    anm_stiffness_index <- rowMeans(anm_stiffness[ , b_s])
    
    anm_effector <- abs(as.matrix(read.table(paste(pdb_sequence[i, 3], "_anm_all_prs_effector.txt", sep = ""))))
    anm_effector_index <- anm_effector[ , 1]
    
    anm_sensor <- abs(as.matrix(read.table(paste(pdb_sequence[i, 3], "_anm_all_prs_sensor.txt", sep = ""))))
    anm_sensor_index <- anm_sensor[ , 1]
    
    gnm_effector <- abs(as.matrix(read.table(paste(pdb_sequence[i, 3], "_gnm_all_prs_effector.txt", sep = ""))))
    gnm_effector_index <- gnm_effector[ , 1]
    
    gnm_sensor <- abs(as.matrix(read.table(paste(pdb_sequence[i, 3], "_gnm_all_prs_sensor.txt", sep = ""))))
    gnm_sensor_index <- gnm_sensor[ , 1]
    
    setwd("D:/subject/dataset/PSP/regulatory_site/network parameter")
    b_factor <- read.table(paste(pdb_sequence[i, 3], "_b-factor.txt", sep = ""))
    b_factor_index <- b_factor[ , 1]
    
    betweenness <- read.table(paste(pdb_sequence[i, 3], "_betweenness.txt", sep = ""))
    betweenness_index <- betweenness[ , 1]
    
    closeness <- read.table(paste(pdb_sequence[i, 3], "_closeness.txt", sep = ""))
    closeness_index <- closeness[ , 1]
    
    degree <- read.table(paste(pdb_sequence[i, 3], "_degree.txt", sep = ""))
    degree_index <- degree[ , 1]
    
    cluster <- read.table(paste(pdb_sequence[i, 3], "_cluster.txt", sep = ""))
    cluster_index <- cluster[ , 1]
    
    diversity <- read.table(paste(pdb_sequence[i, 3], "_diversity.txt", sep = ""))
    diversity_index <- diversity[ , 1]
    
    eccentricity <- read.table(paste(pdb_sequence[i, 3], "_eccentricity.txt", sep = ""))
    eccentricity_index <- as.numeric(eccentricity[ , 1])
    
    strength <- read.table(paste(pdb_sequence[i, 3], "_strength.txt", sep = ""))
    strength_index <- strength[ , 1]
    
    page_rank <- read.table(paste(pdb_sequence[i, 3], "_page_rank.txt", sep = ""))
    page_rank_index <- page_rank[ , 1]
    
    mean_sp_bs <- as.matrix(read.table(paste(pdb_sequence[i, 3], "_shortest_path.txt", sep = "")))
    mean_sp_bs_index <- rowMeans(mean_sp_bs[ , b_s])
    
    ACC <- read.table(paste(pdb_sequence[i, 3], "_ACC.txt", sep = ""))
    ACC_index <- as.numeric(ACC[ , 1])
    
    Phi <- read.table(paste(pdb_sequence[i, 3], "_Phi.txt", sep = ""))
    Phi_index <- Phi[ , 1]
    
    PSi <- read.table(paste(pdb_sequence[i, 3], "_PSi.txt", sep = ""))
    PSi_index <- PSi[ , 1]
    
    setwd("D:/subject/Python script/hit_commute/output")
    hitting_time <- as.matrix(read.table(paste(pdb_sequence[i, 3], "_HitTimes.txt", sep = "")))
    hitting_time_index_col <- colMeans(hitting_time[b_s, ])
    hitting_time_index_row <- rowMeans(hitting_time[, b_s])    
    
    commute_time <- as.matrix(read.table(paste(pdb_sequence[i, 3], "_CommuteTimes.txt", sep = "")))
    commute_time_index <- colMeans(commute_time[b_s, ])
    
    setwd("D:/R project/protein_network/feature")
    Iij_betweenness <- read.table(paste(pdb_sequence[i, 3], "_betweenness.txt", sep = ""))
    Iij_betweenness_index <- Iij_betweenness[ , 1]
    
    Iij_closeness <- read.table(paste(pdb_sequence[i, 3], "_closeness.txt", sep = ""))
    Iij_closeness_index <- Iij_closeness[ , 1]
    
    Iij_degree <- read.table(paste(pdb_sequence[i, 3], "_degree.txt", sep = ""))
    Iij_degree_index <- Iij_degree[ , 1]
    
    Iij_cluster <- read.table(paste(pdb_sequence[i, 3], "_cluster.txt", sep = ""))
    Iij_cluster_index <- Iij_cluster[ , 1]
    
    Iij_diversity <- read.table(paste(pdb_sequence[i, 3], "_diversity.txt", sep = ""))
    Iij_diversity_index <- Iij_diversity[ , 1]
    
    Iij_eccentricity <- read.table(paste(pdb_sequence[i, 3], "_eccentricity.txt", sep = ""))
    Iij_eccentricity_index <- Iij_eccentricity[ , 1]
    
    Iij_strength <- read.table(paste(pdb_sequence[i, 3], "_strength.txt", sep = ""))
    Iij_strength_index <- Iij_strength[ , 1]
    
    Iij_page_rank <- read.table(paste(pdb_sequence[i, 3], "_page_rank.txt", sep = ""))
    Iij_page_rank_index <- Iij_page_rank[ , 1]
    
    setwd("D:/subject/dataset/PSP/regulatory_site/kinase_renumber_pdb/output")
    anm_sq <- read.table(paste(pdb_sequence[i, 3], "_anm_all_sq.txt", sep = ""))
    gnm_sq <- read.table(paste(pdb_sequence[i, 3], "_gnm_all_sq.txt", sep = ""))
    
    pdb_sequence[ , 7:50] <- data.frame(anm_top3_bcc_index,gnm_top3_bcc_index,anm_LF_bcc_index,gnm_LF_bcc_index,
                                        anm_LTIF_bcc_index,gnm_LTIF_bcc_index,anm_HF_bcc_index,gnm_HF_bcc_index,
                                        anm_all_bcc_index,gnm_all_bcc_index,gnm_top3_mode_index,gnm_slow_mode_index,
                                        gnm_all_mode_index,b_factor_index,betweenness_index,closeness_index,degree_index,
                                        cluster_index,diversity_index,eccentricity_index,strength_index,page_rank_index,
                                        mean_sp_bs_index,anm_stiffness_index,ACC_index,Phi_index,PSi_index,
                                        anm_effector_index,anm_sensor_index,gnm_effector_index,gnm_sensor_index,
                                        hitting_time_index_col,commute_time_index,Iij_betweenness_index,Iij_closeness_index,
                                        Iij_degree_index,Iij_cluster_index,Iij_diversity_index,Iij_eccentricity_index,
                                        Iij_strength_index,Iij_page_rank_index,hitting_time_index_row,anm_sq,gnm_sq)
    ki_feature <- rbind(ki_feature, pdb_sequence)
  }
}
ki_feature_1 <- ki_feature[-grep("-", ki_feature$functionality), ]  #delee "-" data
write.csv(ki_feature_1, file = "data_res.csv")  #write feature file

#residue violin plot
data_res_<- read.csv("data_res.csv")  #read feature data
data_res[grep("NFP", data_res$functionality), 6] <- "FP"  #NFP change to FP
data_res[grep("FP", data_res$functionality), 6] <- "PTM"  #FP change to PTM
data_res[grep("Bind", data_res$functionality), 6] <- "Orthosteric"  #Bind change to Orthosteric
data_res$functionality <- factor(data_res$functionality, levels = c("PTM", "Orthosteric", "Other"))
my_comparisons <- list(c("PTM", "Orthosteric"), c("Orthosteric", "Other"), c("PTM", "Other"))
setwd("D:/subject/dataset/PSP/regulatory_site/kinase/violin/2021_1_26")
for (i in 14:77) {
  a <- data_res[ , c(6, i)]
  a[ , 2] <- log(a[ , 2])
  p <- ggviolin(a, x = "functionality", y = colnames(a)[2], fill = "functionality",
                add = "boxplot", add.params = list(fill="white"), palette = c("red", "blue", "green"))+
    stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)+
    xlab("")+
    guides(fill=FALSE)
  p
  ggsave(filename = paste(colnames(a)[2], "png", sep = "."), plot = p, dpi = 600, width = 7, height = 6)
  #dev.off()
  i <- i + 1
}


#pocket level feature
data_poc <- data.frame()
for (i in 1:84) {
  wd <- paste("D:/subject/dataset/PSP/regulatory_site/kinase/fpocket", paste(data_ki[i, 3], "out", sep = "_"), sep = "/")
  setwd(wd)  #set work path
  
  #read the pocket information of the protein
  poc_info <- read.table(paste(data_ki[i, 3], "info_1.txt", sep = "_"), header = FALSE,
                         sep = "\t", fill = TRUE, stringsAsFactors=FALSE)
  for (j in 1:as.numeric(dim(poc_info)[1]/22)) {
    poc_info_1 <- as.data.frame(poc_info[c(((j-1) * 22 + 1):(j * 22)), ], stringsAsFactors=FALSE)
    rownames(poc_info_1) <- NULL
    
    df <- data.frame(FP_number=NA, NP_number=NA, Score=NA, Druggability_Score=NA, Total_SASA=NA, Polar_SASA=NA, Apolar_SASA=NA, volume=NA,
                     Flexibility=NA, anm_top3_bcc_m=NA, anm_top3_bcc_s=NA, anm_LF_bcc_m=NA, anm_LF_bcc_s=NA, anm_LTIF_bcc_m=NA,
                     anm_LTIF_bcc_s=NA, anm_HF_bcc_m=NA, anm_HF_bcc_s=NA, gnm_top3_bcc_m=NA, gnm_top3_bcc_s=NA, gnm_LF_bcc_m=NA,
                     gnm_LF_bcc_s=NA, gnm_LTIF_bcc_m=NA, gnm_LTIF_bcc_s=NA, gnm_HF_bcc_m=NA, gnm_HF_bcc_s=NA, anm_all_sq=NA,
                     gnm_all_sq=NA, anm_prs_sensor=NA, anm_prs_effector=NA, anm_prs_binding_row=NA, anm_prs_binding_col=NA,
                     gnm_prs_sensor=NA, gnm_prs_effector=NA, gnm_prs_binding_row=NA, gnm_prs_binding_col=NA, hit_all=NA,
                     hit_binding_row=NA, hit_binding_col=NA, commute_all=NA, commute_binding=NA, ACC=NA, b_factor=NA,
                     betweenness=NA, closeness=NA, degree=NA, cluster=NA, diversity=NA, eccentricity=NA, strength=NA,
                     anm_stiffness=NA, page_rank=NA, shortest_path=NA, Iij_betweenness=NA, Iij_closeness=NA, Iij_degree=NA,
                     Iij_cluster=NA, Iij_diversity=NA, Iij_eccentricity=NA, Iij_strength=NA, Iij_page_rank=NA, cons_mean=NA,
                     cons_sum=NA,mutinfo_mean=NA,mutinfo_sum=NA,omes_mean=NA,omes_sum=NA,sca_mean=NA,sca_sum=NA,
                     di_mean=NA,di_sum=NA,functionality=1,anm_sq_mean=NA,anm_sq_sum=NA,gnm_sq_mean=NA,gnm_sq_sum=NA,
                     anm_sensor_mean=NA,anm_sensor_sum=NA,anm_effector_mean=NA,anm_effector_sum=NA,anm_all_bcc_m=NA,
                     anm_all_bcc_s=NA,gnm_all_bcc_m=NA, gnm_all_bcc_s=NA,sp_m=NA,sp_s=NA,stringsAsFactors=FALSE)
    df[1, 1] <- strsplit(poc_info_1[21, 1], split = " ")[[1]][4]
    df[1, 2] <- strsplit(poc_info_1[22, 1], split = " ")[[1]][4]
    df[1, 3] <- strsplit(poc_info_1[2, 1], split = " ")[[1]][5]
    df[1, 4] <- strsplit(poc_info_1[3, 1], split = " ")[[1]][6]
    df[1, 5] <- strsplit(poc_info_1[5, 1], split = " ")[[1]][6]
    df[1, 6] <- strsplit(poc_info_1[6, 1], split = " ")[[1]][6]
    df[1, 7] <- strsplit(poc_info_1[7, 1], split = " ")[[1]][6]
    df[1, 8] <- strsplit(poc_info_1[8, 1], split = " ")[[1]][5]
    df[1, 9] <- strsplit(poc_info_1[20, 1], split = " ")[[1]][5]
    wd <- paste("D:/subject/dataset/PSP/regulatory_site/kinase/fpocket",
                paste(data_ki[i, 3],"out/pockets", sep = "_"), sep = "/")
    setwd(wd)  #set work path
    poc_pdb <- read.pdb(paste("pocket", paste(j, "atm.pdb", sep = "_"), sep = ""))  # read pdb file
    swiss_start <- as.numeric(strsplit(as.character(data_only_binding[i, 13]), split = "-")[[1]][1])
    swiss_end <- as.numeric(strsplit(as.character(data_only_binding[i, 13]), split = "-")[[1]][2])
    poc_site <- sort(as.numeric(unique(poc_pdb$atom$resno)))  #pocket residues
    poc_site_1 <- poc_site[poc_site  > swiss_start & poc_site < swiss_end]
    poc_site_2 <- poc_site_1 - swiss_start + 1
    bind_site <- as.numeric(strsplit(as.character(data_ki[i, 16]), split = ",")[[1]]) - swiss_start + 1
    
	#calculate pocket feature
    setwd("D:/subject/dataset/PSP/regulatory_site/kinase_renumber_pdb/output")
    anm_top3_cc <- abs(read.table(paste(as.character(data_ki[i, 3]), "_anm_cc_top3.txt", sep = "")))
    df[1, 10] <- mean(rowMeans((anm_top3_cc[poc_site_2, bind_site])))
    df[1, 11] <- sum(anm_top3_cc[poc_site_2, bind_site])
    
    anm_LF_cc <- abs(read.table(paste(as.character(data_ki[i, 3]), "_anm_cc_LF.txt", sep = "")))
    df[1, 12] <- mean(rowMeans((anm_LF_cc[poc_site_2, bind_site])))
    df[1, 13] <- sum(anm_LF_cc[poc_site_2, bind_site])
    
    anm_LTIF_cc <- abs(read.table(paste(as.character(data_ki[i, 3]), "_anm_cc_LTIF.txt", sep = "")))
    df[1, 14] <- mean(rowMeans((anm_LTIF_cc[poc_site_2, bind_site])))
    df[1, 15] <- sum(anm_LTIF_cc[poc_site_2, bind_site])
    
    anm_HF_cc <- abs(read.table(paste(as.character(data_ki[i, 3]), "_anm_cc_HF.txt", sep = "")))
    df[1, 16] <- mean(rowMeans((anm_HF_cc[poc_site_2, bind_site])))
    df[1, 17] <- sum(anm_HF_cc[poc_site_2, bind_site])
    
    anm_all_cc <- abs(read.table(paste(as.character(data_ki[i, 3]), "_anm_all_cc.txt", sep = "")))
    df[1, 80] <- mean(rowMeans((anm_all_cc[poc_site_2, bind_site])))
    df[1, 81] <- sum(anm_all_cc[poc_site_2, bind_site])
    
    gnm_top3_cc <- abs(read.table(paste(as.character(data_ki[i, 3]), "_gnm_cc_top3.txt", sep = "")))
    df[1, 18] <- mean(rowMeans((gnm_top3_cc[poc_site_2, bind_site])))
    df[1, 19] <- sum(gnm_top3_cc[poc_site_2, bind_site])
    
    gnm_LF_cc <- abs(read.table(paste(as.character(data_ki[i, 3]), "_gnm_cc_LF.txt", sep = "")))
    df[1, 20] <- mean(rowMeans((gnm_LF_cc[poc_site_2, bind_site])))
    df[1, 21] <- sum(gnm_LF_cc[poc_site_2, bind_site])
    
    gnm_LTIF_cc <- abs(read.table(paste(as.character(data_ki[i, 3]), "_gnm_cc_LTIF.txt", sep = "")))
    df[1, 22] <- mean(rowMeans((gnm_LTIF_cc[poc_site_2, bind_site])))
    df[1, 23] <- sum(gnm_LTIF_cc[poc_site_2, bind_site])
    
    gnm_HF_cc <- abs(read.table(paste(as.character(data_ki[i, 3]), "_gnm_cc_HF.txt", sep = "")))
    df[1, 24] <- mean(rowMeans((gnm_HF_cc[poc_site_2, bind_site])))
    df[1, 25] <- sum(gnm_HF_cc[poc_site_2, bind_site])
    
    gnm_all_cc <- abs(read.table(paste(as.character(data_ki[i, 3]), "_gnm_all_cc.txt", sep = "")))
    df[1, 82] <- mean(rowMeans((gnm_all_cc[poc_site_2, bind_site])))
    df[1, 83] <- sum(gnm_all_cc[poc_site_2, bind_site])
    
    anm_sq <- read.table(paste(as.character(data_ki[i, 3]), "_anm_all_sq.txt", sep = ""))
    df[1, 26] <- mean(anm_sq[poc_site_2, ])
    
    gnm_sq <- read.table(paste(as.character(data_ki[i, 3]), "_gnm_all_sq.txt", sep = ""))
    df[1, 27] <-mean(gnm_sq[poc_site_2, ])
    
    anm_sensor <- read.table(paste(as.character(data_ki[i, 3]), "_anm_all_prs_sensor.txt", sep = ""))
    df[1, 28] <- mean(anm_sensor[poc_site_2, ])
    
    anm_effector <- read.table(paste(as.character(data_ki[i, 3]), "_anm_all_prs_effector.txt", sep = ""))
    df[1, 29] <- mean(anm_effector[poc_site_2, ])
    
    anm_prs <- read.table(paste(as.character(data_ki[i, 3]), "_anm_prs.txt", sep = ""))
    df[1, 30] <- mean(rowMeans(anm_prs[bind_site, poc_site_2]))
    df[1, 31] <- mean(rowMeans(anm_prs[poc_site_2, bind_site]))
    
    gnm_sensor <- read.table(paste(as.character(data_ki[i, 3]), "_gnm_all_prs_sensor.txt", sep = ""))
    df[1, 32] <- mean(gnm_sensor[poc_site_2, ])
    
    gnm_effector <- read.table(paste(as.character(data_ki[i, 3]), "_gnm_all_prs_effector.txt", sep = ""))
    df[1, 33] <- mean(gnm_effector[poc_site_2, ])
    
    gnm_prs <- read.table(paste(as.character(data_ki[i, 3]), "_gnm_prs.txt", sep = ""))
    df[1, 34] <- mean(rowMeans(gnm_prs[bind_site, poc_site_2]))
    df[1, 35] <- mean(rowMeans(gnm_prs[poc_site_2, bind_site]))
    
    anm_stiffness <- read.table(paste(as.character(data_ki[i, 3]), "_anm_stiffness.txt", sep = ""))
    df[1, 50] <- mean(rowMeans(anm_stiffness[poc_site_2, bind_site]))
    
    setwd("D:/subject/Python script/hit_commute/output")
    hitting_time <- read.table(paste(as.character(data_ki[i, 3]), "_HitTimes.txt", sep = ""))
    df[1, 36] <- mean(rowMeans(hitting_time[poc_site_2, poc_site_2]))
    df[1, 37] <- mean(rowMeans(hitting_time[bind_site, poc_site_2]))
    df[1, 38] <- mean(rowMeans(hitting_time[poc_site_2, bind_site]))
    
    commute_time <- read.table(paste(as.character(data_ki[i, 3]), "_CommuteTimes.txt", sep = ""))
    df[1, 39] <- mean(rowMeans(commute_time[poc_site_2, poc_site_2]))
    df[1, 40] <- mean(rowMeans(commute_time[poc_site_2, bind_site]))
    
    setwd("D:/subject/数据集/PSP/regulatory_site/network parameter")
    ACC <- read.table(paste(as.character(data_ki[i, 3]), "_ACC.txt", sep = ""))
    df[1, 41] <- mean(ACC[poc_site_2, ])
    
    b_factor <- read.table(paste(as.character(data_ki[i, 3]), "_b-factor.txt", sep = ""))
    df[1, 42] <- mean(b_factor[poc_site_2, ])
    
    betweenness <- read.table(paste(as.character(data_ki[i, 3]), "_betweenness.txt", sep = ""))
    df[1, 43] <- mean(betweenness[poc_site_2, ])
    
    closeness <- read.table(paste(as.character(data_ki[i, 3]), "_closeness.txt", sep = ""))
    df[1, 44] <- mean(closeness[poc_site_2, ])
    
    degree <- read.table(paste(as.character(data_ki[i, 3]), "_degree.txt", sep = ""))
    df[1, 45] <- mean(degree[poc_site_2, ])
    
    cluster <- read.table(paste(as.character(data_ki[i, 3]), "_cluster.txt", sep = ""))
    df[1, 46] <- mean(cluster[poc_site_2, ])
    
    diversity <- read.table(paste(as.character(data_ki[i, 3]), "_diversity.txt", sep = ""))
    df[1, 47] <- mean(diversity[poc_site_2, ])
    
    eccentricity <- read.table(paste(as.character(data_ki[i, 3]), "_eccentricity.txt", sep = ""))
    df[1, 48] <- mean(eccentricity[poc_site_2, ])
    
    strength <- read.table(paste(as.character(data_ki[i, 3]), "_strength.txt", sep = ""))
    df[1, 49] <- mean(strength[poc_site_2, ])
    
    page_rank <- read.table(paste(as.character(data_ki[i, 3]), "_page_rank.txt", sep = ""))
    df[1, 51] <- mean(page_rank[poc_site_2, ])
    
    shortest_path <- read.table(paste(as.character(data_ki[i, 3]), "_shortest_path.txt", sep = ""))
    df[1, 52] <- mean(rowMeans(shortest_path[poc_site_2, bind_site]))
    df[1, 84] <- mean(rowMeans(shortest_path[poc_site_2, bind_site]))
    df[1, 85] <- sum(shortest_path[poc_site_2, bind_site])
    
    setwd("D:/R project/protein_network/feature")
    Iij_betweenness <- read.table(paste(as.character(data_ki[i, 3]), "_betweenness.txt", sep = ""))
    df[1, 53] <- mean(Iij_betweenness[poc_site_2, ])
    
    Iij_closeness <- read.table(paste(as.character(data_ki[i, 3]), "_closeness.txt", sep = ""))
    df[1, 54] <- mean(Iij_closeness[poc_site_2, ])
    
    Iij_degree <- read.table(paste(as.character(data_ki[i, 3]), "_degree.txt", sep = ""))
    df[1, 55] <- mean(Iij_degree[poc_site_2, ])
    
    Iij_cluster <- read.table(paste(as.character(data_ki[i, 3]), "_cluster.txt", sep = ""))
    df[1, 56] <- mean(Iij_cluster[poc_site_2, ], na.rm = TRUE)
    
    Iij_diversity <- read.table(paste(as.character(data_ki[i, 3]), "_diversity.txt", sep = ""))
    df[1, 57] <- mean(Iij_diversity[poc_site_2, ], na.rm = TRUE)
    
    Iij_eccentricity <- read.table(paste(as.character(data_ki[i, 3]), "_eccentricity.txt", sep = ""))
    df[1, 58] <- mean(Iij_eccentricity[poc_site_2, ], na.rm = TRUE)
    
    Iij_strength <- read.table(paste(as.character(data_ki[i, 3]), "_strength.txt", sep = ""))
    df[1, 59] <- mean(Iij_strength[poc_site_2, ], na.rm = TRUE)
    
    Iij_page_rank <- read.table(paste(as.character(data_ki[i, 3]), "_page_rank.txt", sep = ""))
    df[1, 60] <- mean(Iij_page_rank[poc_site_2, ], na.rm = TRUE)
    
    setwd("D:/subject/数据集/PSP/regulatory_site/kinase/sequence_feature")
    conservation <- read.table(paste(as.character(data_ki[i, 1]), "_entropy.txt", sep = ""))
    df[1, 61] <- mean(conservation[poc_site_2, 1])
    df[1, 62] <- sum(conservation[poc_site_2, 1])
    
    mutinfo <- read.table(paste(as.character(data_ki[i, 1]), "_mutinfo.txt", sep = ""))
    mutinfo_nor <- (mutinfo - min(mutinfo))/(max(mutinfo) - min(mutinfo)) 
    df[1, 63] <- mean(rowMeans(mutinfo_nor[poc_site_2, bind_site]))
    df[1, 64] <- sum(mutinfo_nor[poc_site_2, bind_site])
    
    omes <- read.table(paste(as.character(data_ki[i, 1]), "_omes.txt", sep = ""))
    omes_nor <- (omes - min(omes))/(max(omes) - min(omes))  
    df[1, 65] <- mean(rowMeans(omes_nor[poc_site_2, bind_site]))
    df[1, 66] <- sum(omes_nor[poc_site_2, bind_site])
    
    sca <- read.table(paste(as.character(data_ki[i, 1]), "_sca.txt", sep = ""))
    sca_nor <- (sca - min(sca))/(max(sca) - min(sca))  
    df[1, 67] <- mean(rowMeans(sca_nor[poc_site_2, bind_site]))
    df[1, 68] <- sum(sca_nor[poc_site_2, bind_site])
    
    di <- read.table(paste(as.character(data_ki[i, 1]), "_di.txt", sep = ""))
    di_nor <- (di - min(di))/(max(di) - min(di))  
    df[1, 69] <- mean(rowMeans(di_nor[poc_site_2, bind_site]))
    df[1, 70] <- sum(di_nor[poc_site_2, bind_site])
    
    if(j == 1){
      df[1, 71] <- 2
    } else if(df$FP_number == "NA" & df$NP_number == "NA") {
      df[1, 71] <- 0
    }
    
    setwd("D:/subject/数据集/PSP/regulatory_site/output_1")
    anm_sq <- read.table(paste(as.character(data_ki[i, 3]), "_anm_all_sq.txt", sep = ""))
    df[1, 72] <- mean(anm_sq[poc_site_2, 1])
    df[1, 73] <- sum(anm_sq[poc_site_2, 1])
    
    gnm_sq <- read.table(paste(as.character(data_ki[i, 3]), "_gnm_all_sq.txt", sep = ""))
    df[1, 74] <- mean(gnm_sq[poc_site_2, 1])
    df[1, 75] <- sum(gnm_sq[poc_site_2, 1])
    
    anm_sensor <- read.table(paste(as.character(data_ki[i, 3]), "_anm_all_prs_sensor.txt", sep = ""))
    df[1, 76] <- mean(anm_sensor[poc_site_2, 1])
    df[1, 77] <- sum(anm_sensor[poc_site_2, 1])
    
    anm_effector <- read.table(paste(as.character(data_ki[i, 3]), "_anm_all_prs_effector.txt", sep = ""))
    df[1, 78] <- mean(anm_effector[poc_site_2, 1])
    df[1, 79] <- sum(anm_effector[poc_site_2, 1])
    
    rownames(df) <- paste(as.character(data_ki[i, 3]), paste("pocket", j, sep = ""), sep = "_")
    a <- rbind(a, df)
  }
}
data_poc <- data_poc[ , c(1,2,71,3:70,72:85)]  #modify column position

#pocket violin plot
data_poc[grep(0, data_poc$functionality), 3] <- "Non-PTM Pocket"
data_poc[grep(1, data_poc$functionality), 3] <- "PTM Pocket"
data_poc[grep(2, data_poc$functionality), 3] <- "Orthosteric Pocket"
data_poc$functionality <- factor(data_poc$functionality, levels = c("PTM Pocket", "Orthosteric Pocket", "Non-PTM Pocket"))
my_comparisons <- list(c("PTM Pocket", "Orthosteric Pocket"), c("Orthosteric Pocket", "Non-PTM Pocket"), c("PTM Pocket", "Non-PTM Pocket"))
setwd("D:/subject/数据集/PSP/regulatory_site/kinase/violin/2021_1_25")
for (i in 13:79) {
  a <- data_poc[ , c(3, i)]
  #a[ , 2] <- as.numeric(a[ , 2])
  a[ , 2] <- log(a[ , 2])
  p <- ggviolin(a, x = "functionality", y = colnames(a)[2], fill = "functionality",
                palette = c("red", "blue", "green"), add = "boxplot", add.params = list(fill="white"))+
    stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)
  p
  ggsave(filename = paste(colnames(a)[2], "png", sep = "."), plot = p, dpi = 600, width = 7, height = 6)
  dev.off()
}







