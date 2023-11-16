####### KS test & Q1-3 loop   ##### 

n <- 100
tas_1 <- matrix(ncol=20, nrow=n)
tas_2 <- matrix(ncol=20, nrow=n)
tas_3 <- matrix(ncol=20, nrow=n)
tas_4 <- matrix(ncol=20, nrow=n)
tas_5 <- matrix(ncol=20, nrow=n)
tas_6 <- matrix(ncol=20, nrow=n)
tas_7 <- matrix(ncol=20, nrow=n)
tas_8 <- matrix(ncol=20, nrow=n)

pr_1 <- matrix(ncol=20, nrow=n)
pr_2 <- matrix(ncol=20, nrow=n)
pr_3 <- matrix(ncol=20, nrow=n)
pr_4 <- matrix(ncol=20, nrow=n)
pr_5 <- matrix(ncol=20, nrow=n)
pr_6 <- matrix(ncol=20, nrow=n)
pr_7 <- matrix(ncol=20, nrow=n)
pr_8 <- matrix(ncol=20, nrow=n)

snd_1 <- matrix(ncol=20, nrow=n)
snd_2 <- matrix(ncol=20, nrow=n)
snd_3 <- matrix(ncol=20, nrow=n)
snd_4 <- matrix(ncol=20, nrow=n)
snd_5 <- matrix(ncol=20, nrow=n)
snd_6 <- matrix(ncol=20, nrow=n)
snd_7 <- matrix(ncol=20, nrow=n)
snd_8 <- matrix(ncol=20, nrow=n)

mrsos_1 <- matrix(ncol=20, nrow=n)
mrsos_2 <- matrix(ncol=20, nrow=n)
mrsos_3 <- matrix(ncol=20, nrow=n)
mrsos_4 <- matrix(ncol=20, nrow=n)
mrsos_5 <- matrix(ncol=20, nrow=n)
mrsos_6 <- matrix(ncol=20, nrow=n)
mrsos_7 <- matrix(ncol=20, nrow=n)
mrsos_8 <- matrix(ncol=20, nrow=n)

cveg_1 <- matrix(ncol=20, nrow=n)
cveg_2 <- matrix(ncol=20, nrow=n)
cveg_3 <- matrix(ncol=20, nrow=n)
cveg_4 <- matrix(ncol=20, nrow=n)
cveg_5 <- matrix(ncol=20, nrow=n)
cveg_6 <- matrix(ncol=20, nrow=n)
cveg_7 <- matrix(ncol=20, nrow=n)
cveg_8 <- matrix(ncol=20, nrow=n)

csoil_1 <- matrix(ncol=20, nrow=n)
csoil_2 <- matrix(ncol=20, nrow=n)
csoil_3 <- matrix(ncol=20, nrow=n)
csoil_4 <- matrix(ncol=20, nrow=n)
csoil_5 <- matrix(ncol=20, nrow=n)
csoil_6 <- matrix(ncol=20, nrow=n)
csoil_7 <- matrix(ncol=20, nrow=n)
csoil_8 <- matrix(ncol=20, nrow=n)

npp_1 <- matrix(ncol=20, nrow=n)
npp_2 <- matrix(ncol=20, nrow=n)
npp_3 <- matrix(ncol=20, nrow=n)
npp_4 <- matrix(ncol=20, nrow=n)
npp_5 <- matrix(ncol=20, nrow=n)
npp_6 <- matrix(ncol=20, nrow=n)
npp_7 <- matrix(ncol=20, nrow=n)
npp_8 <- matrix(ncol=20, nrow=n)

rh_1 <- matrix(ncol=20, nrow=n)
rh_2 <- matrix(ncol=20, nrow=n)
rh_3 <- matrix(ncol=20, nrow=n)
rh_4 <- matrix(ncol=20, nrow=n)
rh_5 <- matrix(ncol=20, nrow=n)
rh_6 <- matrix(ncol=20, nrow=n)
rh_7 <- matrix(ncol=20, nrow=n)
rh_8 <- matrix(ncol=20, nrow=n)


for (i in 1:n) {     # 100  samples
  ##################################################################################################################################
  #tas model_1_ KS test stats 
  tas_1_ks_2020 <- as.numeric(ks.test(tas_2020_1_INTERACT_sites, sample(tas_2020_1_INTERACT_domain$Value, 496))[c(1:2)])
  tas_1_ks_2020_war <- as.numeric(ks.test(tas_2020_1_INTERACT_sites_war, sample(tas_2020_1_INTERACT_domain$Value, 496))[c(1:2)])
  tas_1_ks_2100 <- as.numeric(ks.test(tas_2100_1_INTERACT_sites, sample(tas_2100_1_INTERACT_domain$Value, 496))[c(1:2)])
  tas_1_ks_2020_2100 <- as.numeric(ks.test(sample(tas_2020_1_INTERACT_domain$Value, 496), sample(tas_2100_1_INTERACT_domain$Value, 496))[c(1:2)])
  #tas model_1_ INTERACT medians
  tas_1_I_25_2020 <- as.numeric(quantile(tas_2020_1_INTERACT_sites,probs = c(0.25),na.rm = T))
  tas_1_I_50_2020 <- as.numeric(median(tas_2020_1_INTERACT_sites,na.rm = T))
  tas_1_I_75_2020 <- as.numeric(quantile(tas_2020_1_INTERACT_sites,probs = c(0.75),na.rm = T))
  #tas model_1_ INTERACT (without Russia) medians
  tas_1_W_25_2020 <- as.numeric(quantile(tas_2020_1_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  tas_1_W_50_2020 <- as.numeric(median(tas_2020_1_INTERACT_sites_war,na.rm = T))
  tas_1_W_75_2020 <- as.numeric(quantile(tas_2020_1_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #tas model_1_ Domain medians & 50% CI 2020
  tas_1_D_25_2020 <- as.numeric(quantile(sample(tas_2020_1_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  tas_1_D_50_2020 <- as.numeric(median(sample(tas_2020_1_INTERACT_domain$Value, 496),na.rm = T))
  tas_1_D_75_2020 <- as.numeric(quantile(sample(tas_2020_1_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #tas model_1_ Domain medians & 50% CI 2100
  tas_1_D_25_2100 <- as.numeric(quantile(sample(tas_2100_1_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  tas_1_D_50_2100 <- as.numeric(median(sample(tas_2100_1_INTERACT_domain$Value, 496),na.rm = T))
  tas_1_D_75_2100 <- as.numeric(quantile(sample(tas_2100_1_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  tas_1[i,] <- round(c(tas_1_ks_2020[1],tas_1_ks_2020[2],tas_1_ks_2020_war[1],tas_1_ks_2020_war[2],
                       tas_1_ks_2100[1],tas_1_ks_2100[2],tas_1_ks_2020_2100[1],tas_1_ks_2020_2100[2],
                       tas_1_I_25_2020,tas_1_I_50_2020,tas_1_I_75_2020,tas_1_W_25_2020,tas_1_W_50_2020,tas_1_W_75_2020,
                       tas_1_D_25_2020,tas_1_D_50_2020,tas_1_D_75_2020,tas_1_D_25_2100,tas_1_D_50_2100,tas_1_D_75_2100),4)
  ##################################################################################################################################
  #tas model_2_ KS test stats 
  tas_2_ks_2020 <- as.numeric(ks.test(tas_2020_2_INTERACT_sites, sample(tas_2020_2_INTERACT_domain$Value, 496))[c(1:2)])
  tas_2_ks_2020_war <- as.numeric(ks.test(tas_2020_2_INTERACT_sites_war, sample(tas_2020_2_INTERACT_domain$Value, 496))[c(1:2)])
  tas_2_ks_2100 <- as.numeric(ks.test(tas_2100_2_INTERACT_sites, sample(tas_2100_2_INTERACT_domain$Value, 496))[c(1:2)])
  tas_2_ks_2020_2100 <- as.numeric(ks.test(sample(tas_2020_2_INTERACT_domain$Value, 496), sample(tas_2100_2_INTERACT_domain$Value, 496))[c(1:2)])
  #tas model_2_ INTERACT medians
  tas_2_I_25_2020 <- as.numeric(quantile(tas_2020_2_INTERACT_sites,probs = c(0.25),na.rm = T))
  tas_2_I_50_2020 <- as.numeric(median(tas_2020_2_INTERACT_sites,na.rm = T))
  tas_2_I_75_2020 <- as.numeric(quantile(tas_2020_2_INTERACT_sites,probs = c(0.75),na.rm = T))
  #tas model_2_ INTERACT (without Russia) medians
  tas_2_W_25_2020 <- as.numeric(quantile(tas_2020_2_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  tas_2_W_50_2020 <- as.numeric(median(tas_2020_2_INTERACT_sites_war,na.rm = T))
  tas_2_W_75_2020 <- as.numeric(quantile(tas_2020_2_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #tas model_2_ Domain medians & 50% CI 2020
  tas_2_D_25_2020 <- as.numeric(quantile(sample(tas_2020_2_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  tas_2_D_50_2020 <- as.numeric(median(sample(tas_2020_2_INTERACT_domain$Value, 496),na.rm = T))
  tas_2_D_75_2020 <- as.numeric(quantile(sample(tas_2020_2_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #tas model_2_ Domain medians & 50% CI 2100
  tas_2_D_25_2100 <- as.numeric(quantile(sample(tas_2100_2_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  tas_2_D_50_2100 <- as.numeric(median(sample(tas_2100_2_INTERACT_domain$Value, 496),na.rm = T))
  tas_2_D_75_2100 <- as.numeric(quantile(sample(tas_2100_2_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  tas_2[i,] <- round(c(tas_2_ks_2020[1],tas_2_ks_2020[2],tas_2_ks_2020_war[1],tas_2_ks_2020_war[2],
                       tas_2_ks_2100[1],tas_2_ks_2100[2],tas_2_ks_2020_2100[1],tas_2_ks_2020_2100[2],
                       tas_2_I_25_2020,tas_2_I_50_2020,tas_2_I_75_2020,tas_2_W_25_2020,tas_2_W_50_2020,tas_2_W_75_2020,
                       tas_2_D_25_2020,tas_2_D_50_2020,tas_2_D_75_2020,tas_2_D_25_2100,tas_2_D_50_2100,tas_2_D_75_2100),4)
  ##############################################################################################################################
  #tas model_3_ KS test stats 
  tas_3_ks_2020 <- as.numeric(ks.test(tas_2020_3_INTERACT_sites, sample(tas_2020_3_INTERACT_domain$Value, 496))[c(1:2)])
  tas_3_ks_2020_war <- as.numeric(ks.test(tas_2020_3_INTERACT_sites_war, sample(tas_2020_3_INTERACT_domain$Value, 496))[c(1:2)])
  tas_3_ks_2100 <- as.numeric(ks.test(tas_2100_3_INTERACT_sites, sample(tas_2100_3_INTERACT_domain$Value, 496))[c(1:2)])
  tas_3_ks_2020_2100 <- as.numeric(ks.test(sample(tas_2020_3_INTERACT_domain$Value, 496), sample(tas_2100_3_INTERACT_domain$Value, 496))[c(1:2)])
  #tas model_3_ INTERACT medians
  tas_3_I_25_2020 <- as.numeric(quantile(tas_2020_3_INTERACT_sites,probs = c(0.25),na.rm = T))
  tas_3_I_50_2020 <- as.numeric(median(tas_2020_3_INTERACT_sites,na.rm = T))
  tas_3_I_75_2020 <- as.numeric(quantile(tas_2020_3_INTERACT_sites,probs = c(0.75),na.rm = T))
  #tas model_3_ INTERACT (without Russia) medians
  tas_3_W_25_2020 <- as.numeric(quantile(tas_2020_3_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  tas_3_W_50_2020 <- as.numeric(median(tas_2020_3_INTERACT_sites_war,na.rm = T))
  tas_3_W_75_2020 <- as.numeric(quantile(tas_2020_3_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #tas model_3_ Domain medians & 50% CI 2020
  tas_3_D_25_2020 <- as.numeric(quantile(sample(tas_2020_3_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  tas_3_D_50_2020 <- as.numeric(median(sample(tas_2020_3_INTERACT_domain$Value, 496),na.rm = T))
  tas_3_D_75_2020 <- as.numeric(quantile(sample(tas_2020_3_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #tas model_3_ Domain medians & 50% CI 2100
  tas_3_D_25_2100 <- as.numeric(quantile(sample(tas_2100_3_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  tas_3_D_50_2100 <- as.numeric(median(sample(tas_2100_3_INTERACT_domain$Value, 496),na.rm = T))
  tas_3_D_75_2100 <- as.numeric(quantile(sample(tas_2100_3_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  tas_3[i,] <- round(c(tas_3_ks_2020[1],tas_3_ks_2020[2],tas_3_ks_2020_war[1],tas_3_ks_2020_war[2],
                       tas_3_ks_2100[1],tas_3_ks_2100[2],tas_3_ks_2020_2100[1],tas_3_ks_2020_2100[2],
                       tas_3_I_25_2020,tas_3_I_50_2020,tas_3_I_75_2020,tas_3_W_25_2020,tas_3_W_50_2020,tas_3_W_75_2020,
                       tas_3_D_25_2020,tas_3_D_50_2020,tas_3_D_75_2020,tas_3_D_25_2100,tas_3_D_50_2100,tas_3_D_75_2100),4)
  ##############################################################################################################################
  #tas model_4_ KS test stats 
  tas_4_ks_2020 <- as.numeric(ks.test(tas_2020_4_INTERACT_sites, sample(tas_2020_4_INTERACT_domain$Value, 496))[c(1:2)])
  tas_4_ks_2020_war <- as.numeric(ks.test(tas_2020_4_INTERACT_sites_war, sample(tas_2020_4_INTERACT_domain$Value, 496))[c(1:2)])
  tas_4_ks_2100 <- as.numeric(ks.test(tas_2100_4_INTERACT_sites, sample(tas_2100_4_INTERACT_domain$Value, 496))[c(1:2)])
  tas_4_ks_2020_2100 <- as.numeric(ks.test(sample(tas_2020_4_INTERACT_domain$Value, 496), sample(tas_2100_4_INTERACT_domain$Value, 496))[c(1:2)])
  #tas model_4_ INTERACT medians
  tas_4_I_25_2020 <- as.numeric(quantile(tas_2020_4_INTERACT_sites,probs = c(0.25),na.rm = T))
  tas_4_I_50_2020 <- as.numeric(median(tas_2020_4_INTERACT_sites,na.rm = T))
  tas_4_I_75_2020 <- as.numeric(quantile(tas_2020_4_INTERACT_sites,probs = c(0.75),na.rm = T))
  #tas model_4_ INTERACT (without Russia) medians
  tas_4_W_25_2020 <- as.numeric(quantile(tas_2020_4_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  tas_4_W_50_2020 <- as.numeric(median(tas_2020_4_INTERACT_sites_war,na.rm = T))
  tas_4_W_75_2020 <- as.numeric(quantile(tas_2020_4_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #tas model_4_ Domain medians & 50% CI 2020
  tas_4_D_25_2020 <- as.numeric(quantile(sample(tas_2020_4_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  tas_4_D_50_2020 <- as.numeric(median(sample(tas_2020_4_INTERACT_domain$Value, 496),na.rm = T))
  tas_4_D_75_2020 <- as.numeric(quantile(sample(tas_2020_4_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #tas model_4_ Domain medians & 50% CI 2100
  tas_4_D_25_2100 <- as.numeric(quantile(sample(tas_2100_4_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  tas_4_D_50_2100 <- as.numeric(median(sample(tas_2100_4_INTERACT_domain$Value, 496),na.rm = T))
  tas_4_D_75_2100 <- as.numeric(quantile(sample(tas_2100_4_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  tas_4[i,] <- round(c(tas_4_ks_2020[1],tas_4_ks_2020[2],tas_4_ks_2020_war[1],tas_4_ks_2020_war[2],
                       tas_4_ks_2100[1],tas_4_ks_2100[2],tas_4_ks_2020_2100[1],tas_4_ks_2020_2100[2],
                       tas_4_I_25_2020,tas_4_I_50_2020,tas_4_I_75_2020,tas_4_W_25_2020,tas_4_W_50_2020,tas_4_W_75_2020,
                       tas_4_D_25_2020,tas_4_D_50_2020,tas_4_D_75_2020,tas_4_D_25_2100,tas_4_D_50_2100,tas_4_D_75_2100),4)
  ##############################################################################################################################
  #tas model_5_ KS test stats 
  tas_5_ks_2020 <- as.numeric(ks.test(tas_2020_5_INTERACT_sites, sample(tas_2020_5_INTERACT_domain$Value, 496))[c(1:2)])
  tas_5_ks_2020_war <- as.numeric(ks.test(tas_2020_5_INTERACT_sites_war, sample(tas_2020_5_INTERACT_domain$Value, 496))[c(1:2)])
  tas_5_ks_2100 <- as.numeric(ks.test(tas_2100_5_INTERACT_sites, sample(tas_2100_5_INTERACT_domain$Value, 496))[c(1:2)])
  tas_5_ks_2020_2100 <- as.numeric(ks.test(sample(tas_2020_5_INTERACT_domain$Value, 496), sample(tas_2100_5_INTERACT_domain$Value, 496))[c(1:2)])
  #tas model_5_ INTERACT medians
  tas_5_I_25_2020 <- as.numeric(quantile(tas_2020_5_INTERACT_sites,probs = c(0.25),na.rm = T))
  tas_5_I_50_2020 <- as.numeric(median(tas_2020_5_INTERACT_sites,na.rm = T))
  tas_5_I_75_2020 <- as.numeric(quantile(tas_2020_5_INTERACT_sites,probs = c(0.75),na.rm = T))
  #tas model_5_ INTERACT (without Russia) medians
  tas_5_W_25_2020 <- as.numeric(quantile(tas_2020_5_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  tas_5_W_50_2020 <- as.numeric(median(tas_2020_5_INTERACT_sites_war,na.rm = T))
  tas_5_W_75_2020 <- as.numeric(quantile(tas_2020_5_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #tas model_5_ Domain medians & 50% CI 2020
  tas_5_D_25_2020 <- as.numeric(quantile(sample(tas_2020_5_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  tas_5_D_50_2020 <- as.numeric(median(sample(tas_2020_5_INTERACT_domain$Value, 496),na.rm = T))
  tas_5_D_75_2020 <- as.numeric(quantile(sample(tas_2020_5_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #tas model_5_ Domain medians & 50% CI 2100
  tas_5_D_25_2100 <- as.numeric(quantile(sample(tas_2100_5_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  tas_5_D_50_2100 <- as.numeric(median(sample(tas_2100_5_INTERACT_domain$Value, 496),na.rm = T))
  tas_5_D_75_2100 <- as.numeric(quantile(sample(tas_2100_5_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  tas_5[i,] <- round(c(tas_5_ks_2020[1],tas_5_ks_2020[2],tas_5_ks_2020_war[1],tas_5_ks_2020_war[2],
                       tas_5_ks_2100[1],tas_5_ks_2100[2],tas_5_ks_2020_2100[1],tas_5_ks_2020_2100[2],
                       tas_5_I_25_2020,tas_5_I_50_2020,tas_5_I_75_2020,tas_5_W_25_2020,tas_5_W_50_2020,tas_5_W_75_2020,
                       tas_5_D_25_2020,tas_5_D_50_2020,tas_5_D_75_2020,tas_5_D_25_2100,tas_5_D_50_2100,tas_5_D_75_2100),4)
  ##############################################################################################################################
  #tas model_6_ KS test stats 
  tas_6_ks_2020 <- as.numeric(ks.test(tas_2020_6_INTERACT_sites, sample(tas_2020_6_INTERACT_domain$Value, 496))[c(1:2)])
  tas_6_ks_2020_war <- as.numeric(ks.test(tas_2020_6_INTERACT_sites_war, sample(tas_2020_6_INTERACT_domain$Value, 496))[c(1:2)])
  tas_6_ks_2100 <- as.numeric(ks.test(tas_2100_6_INTERACT_sites, sample(tas_2100_6_INTERACT_domain$Value, 496))[c(1:2)])
  tas_6_ks_2020_2100 <- as.numeric(ks.test(sample(tas_2020_6_INTERACT_domain$Value, 496), sample(tas_2100_6_INTERACT_domain$Value, 496))[c(1:2)])
  #tas model_6_ INTERACT medians
  tas_6_I_25_2020 <- as.numeric(quantile(tas_2020_6_INTERACT_sites,probs = c(0.25),na.rm = T))
  tas_6_I_50_2020 <- as.numeric(median(tas_2020_6_INTERACT_sites,na.rm = T))
  tas_6_I_75_2020 <- as.numeric(quantile(tas_2020_6_INTERACT_sites,probs = c(0.75),na.rm = T))
  #tas model_6_ INTERACT (without Russia) medians
  tas_6_W_25_2020 <- as.numeric(quantile(tas_2020_6_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  tas_6_W_50_2020 <- as.numeric(median(tas_2020_6_INTERACT_sites_war,na.rm = T))
  tas_6_W_75_2020 <- as.numeric(quantile(tas_2020_6_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #tas model_6_ Domain medians & 50% CI 2020
  tas_6_D_25_2020 <- as.numeric(quantile(sample(tas_2020_6_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  tas_6_D_50_2020 <- as.numeric(median(sample(tas_2020_6_INTERACT_domain$Value, 496),na.rm = T))
  tas_6_D_75_2020 <- as.numeric(quantile(sample(tas_2020_6_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #tas model_6_ Domain medians & 50% CI 2100
  tas_6_D_25_2100 <- as.numeric(quantile(sample(tas_2100_6_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  tas_6_D_50_2100 <- as.numeric(median(sample(tas_2100_6_INTERACT_domain$Value, 496),na.rm = T))
  tas_6_D_75_2100 <- as.numeric(quantile(sample(tas_2100_6_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  tas_6[i,] <- round(c(tas_6_ks_2020[1],tas_6_ks_2020[2],tas_6_ks_2020_war[1],tas_6_ks_2020_war[2],
                       tas_6_ks_2100[1],tas_6_ks_2100[2],tas_6_ks_2020_2100[1],tas_6_ks_2020_2100[2],
                       tas_6_I_25_2020,tas_6_I_50_2020,tas_6_I_75_2020,tas_6_W_25_2020,tas_6_W_50_2020,tas_6_W_75_2020,
                       tas_6_D_25_2020,tas_6_D_50_2020,tas_6_D_75_2020,tas_6_D_25_2100,tas_6_D_50_2100,tas_6_D_75_2100),4)
  ##############################################################################################################################
  #tas model_7_ KS test stats 
  tas_7_ks_2020 <- as.numeric(ks.test(tas_2020_7_INTERACT_sites, sample(tas_2020_7_INTERACT_domain$Value, 496))[c(1:2)])
  tas_7_ks_2020_war <- as.numeric(ks.test(tas_2020_7_INTERACT_sites_war, sample(tas_2020_7_INTERACT_domain$Value, 496))[c(1:2)])
  tas_7_ks_2100 <- as.numeric(ks.test(tas_2100_7_INTERACT_sites, sample(tas_2100_7_INTERACT_domain$Value, 496))[c(1:2)])
  tas_7_ks_2020_2100 <- as.numeric(ks.test(sample(tas_2020_7_INTERACT_domain$Value, 496), sample(tas_2100_7_INTERACT_domain$Value, 496))[c(1:2)])
  #tas model_7_ INTERACT medians
  tas_7_I_25_2020 <- as.numeric(quantile(tas_2020_7_INTERACT_sites,probs = c(0.25),na.rm = T))
  tas_7_I_50_2020 <- as.numeric(median(tas_2020_7_INTERACT_sites,na.rm = T))
  tas_7_I_75_2020 <- as.numeric(quantile(tas_2020_7_INTERACT_sites,probs = c(0.75),na.rm = T))
  #tas model_7_ INTERACT (without Russia) medians
  tas_7_W_25_2020 <- as.numeric(quantile(tas_2020_7_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  tas_7_W_50_2020 <- as.numeric(median(tas_2020_7_INTERACT_sites_war,na.rm = T))
  tas_7_W_75_2020 <- as.numeric(quantile(tas_2020_7_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #tas model_7_ Domain medians & 50% CI 2020
  tas_7_D_25_2020 <- as.numeric(quantile(sample(tas_2020_7_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  tas_7_D_50_2020 <- as.numeric(median(sample(tas_2020_7_INTERACT_domain$Value, 496),na.rm = T))
  tas_7_D_75_2020 <- as.numeric(quantile(sample(tas_2020_7_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #tas model_7_ Domain medians & 50% CI 2100
  tas_7_D_25_2100 <- as.numeric(quantile(sample(tas_2100_7_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  tas_7_D_50_2100 <- as.numeric(median(sample(tas_2100_7_INTERACT_domain$Value, 496),na.rm = T))
  tas_7_D_75_2100 <- as.numeric(quantile(sample(tas_2100_7_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  tas_7[i,] <- round(c(tas_7_ks_2020[1],tas_7_ks_2020[2],tas_7_ks_2020_war[1],tas_7_ks_2020_war[2],
                       tas_7_ks_2100[1],tas_7_ks_2100[2],tas_7_ks_2020_2100[1],tas_7_ks_2020_2100[2],
                       tas_7_I_25_2020,tas_7_I_50_2020,tas_7_I_75_2020,tas_7_W_25_2020,tas_7_W_50_2020,tas_7_W_75_2020,
                       tas_7_D_25_2020,tas_7_D_50_2020,tas_7_D_75_2020,tas_7_D_25_2100,tas_7_D_50_2100,tas_7_D_75_2100),4)
  
  ##################################################################################################################################
  #tas model_8_ KS test stats 
  tas_8_ks_2020 <- as.numeric(ks.test(tas_2020_8_INTERACT_sites, sample(tas_2020_8_INTERACT_domain$Value, 496))[c(1:2)])
  tas_8_ks_2020_war <- as.numeric(ks.test(tas_2020_8_INTERACT_sites_war, sample(tas_2020_8_INTERACT_domain$Value, 496))[c(1:2)])
  tas_8_ks_2100 <- as.numeric(ks.test(tas_2100_8_INTERACT_sites, sample(tas_2100_8_INTERACT_domain$Value, 496))[c(1:2)])
  tas_8_ks_2020_2100 <- as.numeric(ks.test(sample(tas_2020_8_INTERACT_domain$Value, 496), sample(tas_2100_8_INTERACT_domain$Value, 496))[c(1:2)])
  #tas model_8_ INTERACT medians
  tas_8_I_25_2020 <- as.numeric(quantile(tas_2020_8_INTERACT_sites,probs = c(0.25),na.rm = T))
  tas_8_I_50_2020 <- as.numeric(median(tas_2020_8_INTERACT_sites,na.rm = T))
  tas_8_I_75_2020 <- as.numeric(quantile(tas_2020_8_INTERACT_sites,probs = c(0.75),na.rm = T))
  #tas model_8_ INTERACT (without Russia) medians
  tas_8_W_25_2020 <- as.numeric(quantile(tas_2020_8_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  tas_8_W_50_2020 <- as.numeric(median(tas_2020_8_INTERACT_sites_war,na.rm = T))
  tas_8_W_75_2020 <- as.numeric(quantile(tas_2020_8_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #tas model_8_ Domain medians & 50% CI 2020
  tas_8_D_25_2020 <- as.numeric(quantile(sample(tas_2020_8_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  tas_8_D_50_2020 <- as.numeric(median(sample(tas_2020_8_INTERACT_domain$Value, 496),na.rm = T))
  tas_8_D_75_2020 <- as.numeric(quantile(sample(tas_2020_8_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #tas model_8_ Domain medians & 50% CI 2100
  tas_8_D_25_2100 <- as.numeric(quantile(sample(tas_2100_8_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  tas_8_D_50_2100 <- as.numeric(median(sample(tas_2100_8_INTERACT_domain$Value, 496),na.rm = T))
  tas_8_D_75_2100 <- as.numeric(quantile(sample(tas_2100_8_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  tas_8[i,] <- round(c(tas_8_ks_2020[1],tas_8_ks_2020[2],tas_8_ks_2020_war[1],tas_8_ks_2020_war[2],
                       tas_8_ks_2100[1],tas_8_ks_2100[2],tas_8_ks_2020_2100[1],tas_8_ks_2020_2100[2],
                       tas_8_I_25_2020,tas_8_I_50_2020,tas_8_I_75_2020,tas_8_W_25_2020,tas_8_W_50_2020,tas_8_W_75_2020,
                       tas_8_D_25_2020,tas_8_D_50_2020,tas_8_D_75_2020,tas_8_D_25_2100,tas_8_D_50_2100,tas_8_D_75_2100),4)
  
  ##################################################################################################################################
  #pr model_1_ KS test stats 
  pr_1_ks_2020 <- as.numeric(ks.test(pr_2020_1_INTERACT_sites, sample(pr_2020_1_INTERACT_domain$Value, 496))[c(1:2)])
  pr_1_ks_2020_war <- as.numeric(ks.test(pr_2020_1_INTERACT_sites_war, sample(pr_2020_1_INTERACT_domain$Value, 496))[c(1:2)])
  pr_1_ks_2100 <- as.numeric(ks.test(pr_2100_1_INTERACT_sites, sample(pr_2100_1_INTERACT_domain$Value, 496))[c(1:2)])
  pr_1_ks_2020_2100 <- as.numeric(ks.test(sample(pr_2020_1_INTERACT_domain$Value, 496), sample(pr_2100_1_INTERACT_domain$Value, 496))[c(1:2)])
  #pr model_1_ INTERACT medians
  pr_1_I_25_2020 <- as.numeric(quantile(pr_2020_1_INTERACT_sites,probs = c(0.25),na.rm = T))
  pr_1_I_50_2020 <- as.numeric(median(pr_2020_1_INTERACT_sites,na.rm = T))
  pr_1_I_75_2020 <- as.numeric(quantile(pr_2020_1_INTERACT_sites,probs = c(0.75),na.rm = T))
  #pr model_1_ INTERACT (without Russia) medians
  pr_1_W_25_2020 <- as.numeric(quantile(pr_2020_1_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  pr_1_W_50_2020 <- as.numeric(median(pr_2020_1_INTERACT_sites_war,na.rm = T))
  pr_1_W_75_2020 <- as.numeric(quantile(pr_2020_1_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #pr model_1_ Domain medians & 50% CI 2020
  pr_1_D_25_2020 <- as.numeric(quantile(sample(pr_2020_1_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  pr_1_D_50_2020 <- as.numeric(median(sample(pr_2020_1_INTERACT_domain$Value, 496),na.rm = T))
  pr_1_D_75_2020 <- as.numeric(quantile(sample(pr_2020_1_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #pr model_1_ Domain medians & 50% CI 2100
  pr_1_D_25_2100 <- as.numeric(quantile(sample(pr_2100_1_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  pr_1_D_50_2100 <- as.numeric(median(sample(pr_2100_1_INTERACT_domain$Value, 496),na.rm = T))
  pr_1_D_75_2100 <- as.numeric(quantile(sample(pr_2100_1_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  pr_1[i,] <- round(c(pr_1_ks_2020[1],pr_1_ks_2020[2],pr_1_ks_2020_war[1],pr_1_ks_2020_war[2],
                      pr_1_ks_2100[1],pr_1_ks_2100[2],pr_1_ks_2020_2100[1],pr_1_ks_2020_2100[2],
                      pr_1_I_25_2020,pr_1_I_50_2020,pr_1_I_75_2020,pr_1_W_25_2020,pr_1_W_50_2020,pr_1_W_75_2020,
                      pr_1_D_25_2020,pr_1_D_50_2020,pr_1_D_75_2020,pr_1_D_25_2100,pr_1_D_50_2100,pr_1_D_75_2100),4)
  ##################################################################################################################################
  #pr model_2_ KS test stats 
  pr_2_ks_2020 <- as.numeric(ks.test(pr_2020_2_INTERACT_sites, sample(pr_2020_2_INTERACT_domain$Value, 496))[c(1:2)])
  pr_2_ks_2020_war <- as.numeric(ks.test(pr_2020_2_INTERACT_sites_war, sample(pr_2020_2_INTERACT_domain$Value, 496))[c(1:2)])
  pr_2_ks_2100 <- as.numeric(ks.test(pr_2100_2_INTERACT_sites, sample(pr_2100_2_INTERACT_domain$Value, 496))[c(1:2)])
  pr_2_ks_2020_2100 <- as.numeric(ks.test(sample(pr_2020_2_INTERACT_domain$Value, 496), sample(pr_2100_2_INTERACT_domain$Value, 496))[c(1:2)])
  #pr model_2_ INTERACT medians
  pr_2_I_25_2020 <- as.numeric(quantile(pr_2020_2_INTERACT_sites,probs = c(0.25),na.rm = T))
  pr_2_I_50_2020 <- as.numeric(median(pr_2020_2_INTERACT_sites,na.rm = T))
  pr_2_I_75_2020 <- as.numeric(quantile(pr_2020_2_INTERACT_sites,probs = c(0.75),na.rm = T))
  #pr model_2_ INTERACT (without Russia) medians
  pr_2_W_25_2020 <- as.numeric(quantile(pr_2020_2_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  pr_2_W_50_2020 <- as.numeric(median(pr_2020_2_INTERACT_sites_war,na.rm = T))
  pr_2_W_75_2020 <- as.numeric(quantile(pr_2020_2_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #pr model_2_ Domain medians & 50% CI 2020
  pr_2_D_25_2020 <- as.numeric(quantile(sample(pr_2020_2_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  pr_2_D_50_2020 <- as.numeric(median(sample(pr_2020_2_INTERACT_domain$Value, 496),na.rm = T))
  pr_2_D_75_2020 <- as.numeric(quantile(sample(pr_2020_2_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #pr model_2_ Domain medians & 50% CI 2100
  pr_2_D_25_2100 <- as.numeric(quantile(sample(pr_2100_2_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  pr_2_D_50_2100 <- as.numeric(median(sample(pr_2100_2_INTERACT_domain$Value, 496),na.rm = T))
  pr_2_D_75_2100 <- as.numeric(quantile(sample(pr_2100_2_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  pr_2[i,] <- round(c(pr_2_ks_2020[1],pr_2_ks_2020[2],pr_2_ks_2020_war[1],pr_2_ks_2020_war[2],
                      pr_2_ks_2100[1],pr_2_ks_2100[2],pr_2_ks_2020_2100[1],pr_2_ks_2020_2100[2],
                      pr_2_I_25_2020,pr_2_I_50_2020,pr_2_I_75_2020,pr_2_W_25_2020,pr_2_W_50_2020,pr_2_W_75_2020,
                      pr_2_D_25_2020,pr_2_D_50_2020,pr_2_D_75_2020,pr_2_D_25_2100,pr_2_D_50_2100,pr_2_D_75_2100),4)
  ##############################################################################################################################
  #pr model_3_ KS test stats 
  pr_3_ks_2020 <- as.numeric(ks.test(pr_2020_3_INTERACT_sites, sample(pr_2020_3_INTERACT_domain$Value, 496))[c(1:2)])
  pr_3_ks_2020_war <- as.numeric(ks.test(pr_2020_3_INTERACT_sites_war, sample(pr_2020_3_INTERACT_domain$Value, 496))[c(1:2)])
  pr_3_ks_2100 <- as.numeric(ks.test(pr_2100_3_INTERACT_sites, sample(pr_2100_3_INTERACT_domain$Value, 496))[c(1:2)])
  pr_3_ks_2020_2100 <- as.numeric(ks.test(sample(pr_2020_3_INTERACT_domain$Value, 496), sample(pr_2100_3_INTERACT_domain$Value, 496))[c(1:2)])
  #pr model_3_ INTERACT medians
  pr_3_I_25_2020 <- as.numeric(quantile(pr_2020_3_INTERACT_sites,probs = c(0.25),na.rm = T))
  pr_3_I_50_2020 <- as.numeric(median(pr_2020_3_INTERACT_sites,na.rm = T))
  pr_3_I_75_2020 <- as.numeric(quantile(pr_2020_3_INTERACT_sites,probs = c(0.75),na.rm = T))
  #pr model_3_ INTERACT (without Russia) medians
  pr_3_W_25_2020 <- as.numeric(quantile(pr_2020_3_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  pr_3_W_50_2020 <- as.numeric(median(pr_2020_3_INTERACT_sites_war,na.rm = T))
  pr_3_W_75_2020 <- as.numeric(quantile(pr_2020_3_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #pr model_3_ Domain medians & 50% CI 2020
  pr_3_D_25_2020 <- as.numeric(quantile(sample(pr_2020_3_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  pr_3_D_50_2020 <- as.numeric(median(sample(pr_2020_3_INTERACT_domain$Value, 496),na.rm = T))
  pr_3_D_75_2020 <- as.numeric(quantile(sample(pr_2020_3_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #pr model_3_ Domain medians & 50% CI 2100
  pr_3_D_25_2100 <- as.numeric(quantile(sample(pr_2100_3_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  pr_3_D_50_2100 <- as.numeric(median(sample(pr_2100_3_INTERACT_domain$Value, 496),na.rm = T))
  pr_3_D_75_2100 <- as.numeric(quantile(sample(pr_2100_3_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  pr_3[i,] <- round(c(pr_3_ks_2020[1],pr_3_ks_2020[2],pr_3_ks_2020_war[1],pr_3_ks_2020_war[2],
                      pr_3_ks_2100[1],pr_3_ks_2100[2],pr_3_ks_2020_2100[1],pr_3_ks_2020_2100[2],
                      pr_3_I_25_2020,pr_3_I_50_2020,pr_3_I_75_2020,pr_3_W_25_2020,pr_3_W_50_2020,pr_3_W_75_2020,
                      pr_3_D_25_2020,pr_3_D_50_2020,pr_3_D_75_2020,pr_3_D_25_2100,pr_3_D_50_2100,pr_3_D_75_2100),4)
  ##############################################################################################################################
  #pr model_4_ KS test stats 
  pr_4_ks_2020 <- as.numeric(ks.test(pr_2020_4_INTERACT_sites, sample(pr_2020_4_INTERACT_domain$Value, 496))[c(1:2)])
  pr_4_ks_2020_war <- as.numeric(ks.test(pr_2020_4_INTERACT_sites_war, sample(pr_2020_4_INTERACT_domain$Value, 496))[c(1:2)])
  pr_4_ks_2100 <- as.numeric(ks.test(pr_2100_4_INTERACT_sites, sample(pr_2100_4_INTERACT_domain$Value, 496))[c(1:2)])
  pr_4_ks_2020_2100 <- as.numeric(ks.test(sample(pr_2020_4_INTERACT_domain$Value, 496), sample(pr_2100_4_INTERACT_domain$Value, 496))[c(1:2)])
  #pr model_4_ INTERACT medians
  pr_4_I_25_2020 <- as.numeric(quantile(pr_2020_4_INTERACT_sites,probs = c(0.25),na.rm = T))
  pr_4_I_50_2020 <- as.numeric(median(pr_2020_4_INTERACT_sites,na.rm = T))
  pr_4_I_75_2020 <- as.numeric(quantile(pr_2020_4_INTERACT_sites,probs = c(0.75),na.rm = T))
  #pr model_4_ INTERACT (without Russia) medians
  pr_4_W_25_2020 <- as.numeric(quantile(pr_2020_4_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  pr_4_W_50_2020 <- as.numeric(median(pr_2020_4_INTERACT_sites_war,na.rm = T))
  pr_4_W_75_2020 <- as.numeric(quantile(pr_2020_4_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #pr model_4_ Domain medians & 50% CI 2020
  pr_4_D_25_2020 <- as.numeric(quantile(sample(pr_2020_4_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  pr_4_D_50_2020 <- as.numeric(median(sample(pr_2020_4_INTERACT_domain$Value, 496),na.rm = T))
  pr_4_D_75_2020 <- as.numeric(quantile(sample(pr_2020_4_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #pr model_4_ Domain medians & 50% CI 2100
  pr_4_D_25_2100 <- as.numeric(quantile(sample(pr_2100_4_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  pr_4_D_50_2100 <- as.numeric(median(sample(pr_2100_4_INTERACT_domain$Value, 496),na.rm = T))
  pr_4_D_75_2100 <- as.numeric(quantile(sample(pr_2100_4_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  pr_4[i,] <- round(c(pr_4_ks_2020[1],pr_4_ks_2020[2],pr_4_ks_2020_war[1],pr_4_ks_2020_war[2],
                      pr_4_ks_2100[1],pr_4_ks_2100[2],pr_4_ks_2020_2100[1],pr_4_ks_2020_2100[2],
                      pr_4_I_25_2020,pr_4_I_50_2020,pr_4_I_75_2020,pr_4_W_25_2020,pr_4_W_50_2020,pr_4_W_75_2020,
                      pr_4_D_25_2020,pr_4_D_50_2020,pr_4_D_75_2020,pr_4_D_25_2100,pr_4_D_50_2100,pr_4_D_75_2100),4)
  ##############################################################################################################################
  #pr model_5_ KS test stats 
  pr_5_ks_2020 <- as.numeric(ks.test(pr_2020_5_INTERACT_sites, sample(pr_2020_5_INTERACT_domain$Value, 496))[c(1:2)])
  pr_5_ks_2020_war <- as.numeric(ks.test(pr_2020_5_INTERACT_sites_war, sample(pr_2020_5_INTERACT_domain$Value, 496))[c(1:2)])
  pr_5_ks_2100 <- as.numeric(ks.test(pr_2100_5_INTERACT_sites, sample(pr_2100_5_INTERACT_domain$Value, 496))[c(1:2)])
  pr_5_ks_2020_2100 <- as.numeric(ks.test(sample(pr_2020_5_INTERACT_domain$Value, 496), sample(pr_2100_5_INTERACT_domain$Value, 496))[c(1:2)])
  #pr model_5_ INTERACT medians
  pr_5_I_25_2020 <- as.numeric(quantile(pr_2020_5_INTERACT_sites,probs = c(0.25),na.rm = T))
  pr_5_I_50_2020 <- as.numeric(median(pr_2020_5_INTERACT_sites,na.rm = T))
  pr_5_I_75_2020 <- as.numeric(quantile(pr_2020_5_INTERACT_sites,probs = c(0.75),na.rm = T))
  #pr model_5_ INTERACT (without Russia) medians
  pr_5_W_25_2020 <- as.numeric(quantile(pr_2020_5_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  pr_5_W_50_2020 <- as.numeric(median(pr_2020_5_INTERACT_sites_war,na.rm = T))
  pr_5_W_75_2020 <- as.numeric(quantile(pr_2020_5_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #pr model_5_ Domain medians & 50% CI 2020
  pr_5_D_25_2020 <- as.numeric(quantile(sample(pr_2020_5_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  pr_5_D_50_2020 <- as.numeric(median(sample(pr_2020_5_INTERACT_domain$Value, 496),na.rm = T))
  pr_5_D_75_2020 <- as.numeric(quantile(sample(pr_2020_5_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #pr model_5_ Domain medians & 50% CI 2100
  pr_5_D_25_2100 <- as.numeric(quantile(sample(pr_2100_5_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  pr_5_D_50_2100 <- as.numeric(median(sample(pr_2100_5_INTERACT_domain$Value, 496),na.rm = T))
  pr_5_D_75_2100 <- as.numeric(quantile(sample(pr_2100_5_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  pr_5[i,] <- round(c(pr_5_ks_2020[1],pr_5_ks_2020[2],pr_5_ks_2020_war[1],pr_5_ks_2020_war[2],
                      pr_5_ks_2100[1],pr_5_ks_2100[2],pr_5_ks_2020_2100[1],pr_5_ks_2020_2100[2],
                      pr_5_I_25_2020,pr_5_I_50_2020,pr_5_I_75_2020,pr_5_W_25_2020,pr_5_W_50_2020,pr_5_W_75_2020,
                      pr_5_D_25_2020,pr_5_D_50_2020,pr_5_D_75_2020,pr_5_D_25_2100,pr_5_D_50_2100,pr_5_D_75_2100),4)
  ##############################################################################################################################
  #pr model_6_ KS test stats 
  pr_6_ks_2020 <- as.numeric(ks.test(pr_2020_6_INTERACT_sites, sample(pr_2020_6_INTERACT_domain$Value, 496))[c(1:2)])
  pr_6_ks_2020_war <- as.numeric(ks.test(pr_2020_6_INTERACT_sites_war, sample(pr_2020_6_INTERACT_domain$Value, 496))[c(1:2)])
  pr_6_ks_2100 <- as.numeric(ks.test(pr_2100_6_INTERACT_sites, sample(pr_2100_6_INTERACT_domain$Value, 496))[c(1:2)])
  pr_6_ks_2020_2100 <- as.numeric(ks.test(sample(pr_2020_6_INTERACT_domain$Value, 496), sample(pr_2100_6_INTERACT_domain$Value, 496))[c(1:2)])
  #pr model_6_ INTERACT medians
  pr_6_I_25_2020 <- as.numeric(quantile(pr_2020_6_INTERACT_sites,probs = c(0.25),na.rm = T))
  pr_6_I_50_2020 <- as.numeric(median(pr_2020_6_INTERACT_sites,na.rm = T))
  pr_6_I_75_2020 <- as.numeric(quantile(pr_2020_6_INTERACT_sites,probs = c(0.75),na.rm = T))
  #pr model_6_ INTERACT (without Russia) medians
  pr_6_W_25_2020 <- as.numeric(quantile(pr_2020_6_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  pr_6_W_50_2020 <- as.numeric(median(pr_2020_6_INTERACT_sites_war,na.rm = T))
  pr_6_W_75_2020 <- as.numeric(quantile(pr_2020_6_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #pr model_6_ Domain medians & 50% CI 2020
  pr_6_D_25_2020 <- as.numeric(quantile(sample(pr_2020_6_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  pr_6_D_50_2020 <- as.numeric(median(sample(pr_2020_6_INTERACT_domain$Value, 496),na.rm = T))
  pr_6_D_75_2020 <- as.numeric(quantile(sample(pr_2020_6_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #pr model_6_ Domain medians & 50% CI 2100
  pr_6_D_25_2100 <- as.numeric(quantile(sample(pr_2100_6_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  pr_6_D_50_2100 <- as.numeric(median(sample(pr_2100_6_INTERACT_domain$Value, 496),na.rm = T))
  pr_6_D_75_2100 <- as.numeric(quantile(sample(pr_2100_6_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  pr_6[i,] <- round(c(pr_6_ks_2020[1],pr_6_ks_2020[2],pr_6_ks_2020_war[1],pr_6_ks_2020_war[2],
                      pr_6_ks_2100[1],pr_6_ks_2100[2],pr_6_ks_2020_2100[1],pr_6_ks_2020_2100[2],
                      pr_6_I_25_2020,pr_6_I_50_2020,pr_6_I_75_2020,pr_6_W_25_2020,pr_6_W_50_2020,pr_6_W_75_2020,
                      pr_6_D_25_2020,pr_6_D_50_2020,pr_6_D_75_2020,pr_6_D_25_2100,pr_6_D_50_2100,pr_6_D_75_2100),4)
  ##############################################################################################################################
  #pr model_7_ KS test stats 
  pr_7_ks_2020 <- as.numeric(ks.test(pr_2020_7_INTERACT_sites, sample(pr_2020_7_INTERACT_domain$Value, 496))[c(1:2)])
  pr_7_ks_2020_war <- as.numeric(ks.test(pr_2020_7_INTERACT_sites_war, sample(pr_2020_7_INTERACT_domain$Value, 496))[c(1:2)])
  pr_7_ks_2100 <- as.numeric(ks.test(pr_2100_7_INTERACT_sites, sample(pr_2100_7_INTERACT_domain$Value, 496))[c(1:2)])
  pr_7_ks_2020_2100 <- as.numeric(ks.test(sample(pr_2020_7_INTERACT_domain$Value, 496), sample(pr_2100_7_INTERACT_domain$Value, 496))[c(1:2)])
  #pr model_7_ INTERACT medians
  pr_7_I_25_2020 <- as.numeric(quantile(pr_2020_7_INTERACT_sites,probs = c(0.25),na.rm = T))
  pr_7_I_50_2020 <- as.numeric(median(pr_2020_7_INTERACT_sites,na.rm = T))
  pr_7_I_75_2020 <- as.numeric(quantile(pr_2020_7_INTERACT_sites,probs = c(0.75),na.rm = T))
  #pr model_7_ INTERACT (without Russia) medians
  pr_7_W_25_2020 <- as.numeric(quantile(pr_2020_7_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  pr_7_W_50_2020 <- as.numeric(median(pr_2020_7_INTERACT_sites_war,na.rm = T))
  pr_7_W_75_2020 <- as.numeric(quantile(pr_2020_7_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #pr model_7_ Domain medians & 50% CI 2020
  pr_7_D_25_2020 <- as.numeric(quantile(sample(pr_2020_7_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  pr_7_D_50_2020 <- as.numeric(median(sample(pr_2020_7_INTERACT_domain$Value, 496),na.rm = T))
  pr_7_D_75_2020 <- as.numeric(quantile(sample(pr_2020_7_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #pr model_7_ Domain medians & 50% CI 2100
  pr_7_D_25_2100 <- as.numeric(quantile(sample(pr_2100_7_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  pr_7_D_50_2100 <- as.numeric(median(sample(pr_2100_7_INTERACT_domain$Value, 496),na.rm = T))
  pr_7_D_75_2100 <- as.numeric(quantile(sample(pr_2100_7_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  pr_7[i,] <- round(c(pr_7_ks_2020[1],pr_7_ks_2020[2],pr_7_ks_2020_war[1],pr_7_ks_2020_war[2],
                      pr_7_ks_2100[1],pr_7_ks_2100[2],pr_7_ks_2020_2100[1],pr_7_ks_2020_2100[2],
                      pr_7_I_25_2020,pr_7_I_50_2020,pr_7_I_75_2020,pr_7_W_25_2020,pr_7_W_50_2020,pr_7_W_75_2020,
                      pr_7_D_25_2020,pr_7_D_50_2020,pr_7_D_75_2020,pr_7_D_25_2100,pr_7_D_50_2100,pr_7_D_75_2100),4)
  
  ##################################################################################################################################
  #pr model_8_ KS test stats 
  pr_8_ks_2020 <- as.numeric(ks.test(pr_2020_8_INTERACT_sites, sample(pr_2020_8_INTERACT_domain$Value, 496))[c(1:2)])
  pr_8_ks_2020_war <- as.numeric(ks.test(pr_2020_8_INTERACT_sites_war, sample(pr_2020_8_INTERACT_domain$Value, 496))[c(1:2)])
  pr_8_ks_2100 <- as.numeric(ks.test(pr_2100_8_INTERACT_sites, sample(pr_2100_8_INTERACT_domain$Value, 496))[c(1:2)])
  pr_8_ks_2020_2100 <- as.numeric(ks.test(sample(pr_2020_8_INTERACT_domain$Value, 496), sample(pr_2100_8_INTERACT_domain$Value, 496))[c(1:2)])
  #pr model_8_ INTERACT medians
  pr_8_I_25_2020 <- as.numeric(quantile(pr_2020_8_INTERACT_sites,probs = c(0.25),na.rm = T))
  pr_8_I_50_2020 <- as.numeric(median(pr_2020_8_INTERACT_sites,na.rm = T))
  pr_8_I_75_2020 <- as.numeric(quantile(pr_2020_8_INTERACT_sites,probs = c(0.75),na.rm = T))
  #pr model_8_ INTERACT (without Russia) medians
  pr_8_W_25_2020 <- as.numeric(quantile(pr_2020_8_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  pr_8_W_50_2020 <- as.numeric(median(pr_2020_8_INTERACT_sites_war,na.rm = T))
  pr_8_W_75_2020 <- as.numeric(quantile(pr_2020_8_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #pr model_8_ Domain medians & 50% CI 2020
  pr_8_D_25_2020 <- as.numeric(quantile(sample(pr_2020_8_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  pr_8_D_50_2020 <- as.numeric(median(sample(pr_2020_8_INTERACT_domain$Value, 496),na.rm = T))
  pr_8_D_75_2020 <- as.numeric(quantile(sample(pr_2020_8_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #pr model_8_ Domain medians & 50% CI 2100
  pr_8_D_25_2100 <- as.numeric(quantile(sample(pr_2100_8_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  pr_8_D_50_2100 <- as.numeric(median(sample(pr_2100_8_INTERACT_domain$Value, 496),na.rm = T))
  pr_8_D_75_2100 <- as.numeric(quantile(sample(pr_2100_8_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  pr_8[i,] <- round(c(pr_8_ks_2020[1],pr_8_ks_2020[2],pr_8_ks_2020_war[1],pr_8_ks_2020_war[2],
                      pr_8_ks_2100[1],pr_8_ks_2100[2],pr_8_ks_2020_2100[1],pr_8_ks_2020_2100[2],
                      pr_8_I_25_2020,pr_8_I_50_2020,pr_8_I_75_2020,pr_8_W_25_2020,pr_8_W_50_2020,pr_8_W_75_2020,
                      pr_8_D_25_2020,pr_8_D_50_2020,pr_8_D_75_2020,pr_8_D_25_2100,pr_8_D_50_2100,pr_8_D_75_2100),4)
  ##################################################################################################################################
  #snd model_1_ KS test stats 
  snd_1_ks_2020 <- as.numeric(ks.test(snd_2020_1_INTERACT_sites, sample(snd_2020_1_INTERACT_domain$Value, 496))[c(1:2)])
  snd_1_ks_2020_war <- as.numeric(ks.test(snd_2020_1_INTERACT_sites_war, sample(snd_2020_1_INTERACT_domain$Value, 496))[c(1:2)])
  snd_1_ks_2100 <- as.numeric(ks.test(snd_2100_1_INTERACT_sites, sample(snd_2100_1_INTERACT_domain$Value, 496))[c(1:2)])
  snd_1_ks_2020_2100 <- as.numeric(ks.test(sample(snd_2020_1_INTERACT_domain$Value, 496), sample(snd_2100_1_INTERACT_domain$Value, 496))[c(1:2)])
  #snd model_1_ INTERACT medians
  snd_1_I_25_2020 <- as.numeric(quantile(snd_2020_1_INTERACT_sites,probs = c(0.25),na.rm = T))
  snd_1_I_50_2020 <- as.numeric(median(snd_2020_1_INTERACT_sites,na.rm = T))
  snd_1_I_75_2020 <- as.numeric(quantile(snd_2020_1_INTERACT_sites,probs = c(0.75),na.rm = T))
  #snd model_1_ INTERACT (without Russia) medians
  snd_1_W_25_2020 <- as.numeric(quantile(snd_2020_1_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  snd_1_W_50_2020 <- as.numeric(median(snd_2020_1_INTERACT_sites_war,na.rm = T))
  snd_1_W_75_2020 <- as.numeric(quantile(snd_2020_1_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #snd model_1_ Domain medians & 50% CI 2020
  snd_1_D_25_2020 <- as.numeric(quantile(sample(snd_2020_1_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  snd_1_D_50_2020 <- as.numeric(median(sample(snd_2020_1_INTERACT_domain$Value, 496),na.rm = T))
  snd_1_D_75_2020 <- as.numeric(quantile(sample(snd_2020_1_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #snd model_1_ Domain medians & 50% CI 2100
  snd_1_D_25_2100 <- as.numeric(quantile(sample(snd_2100_1_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  snd_1_D_50_2100 <- as.numeric(median(sample(snd_2100_1_INTERACT_domain$Value, 496),na.rm = T))
  snd_1_D_75_2100 <- as.numeric(quantile(sample(snd_2100_1_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  snd_1[i,] <- round(c(snd_1_ks_2020[1],snd_1_ks_2020[2],snd_1_ks_2020_war[1],snd_1_ks_2020_war[2],
                       snd_1_ks_2100[1],snd_1_ks_2100[2],snd_1_ks_2020_2100[1],snd_1_ks_2020_2100[2],
                       snd_1_I_25_2020,snd_1_I_50_2020,snd_1_I_75_2020,snd_1_W_25_2020,snd_1_W_50_2020,snd_1_W_75_2020,
                       snd_1_D_25_2020,snd_1_D_50_2020,snd_1_D_75_2020,snd_1_D_25_2100,snd_1_D_50_2100,snd_1_D_75_2100),4)
  ##################################################################################################################################
  #snd model_2_ KS test stats 
  snd_2_ks_2020 <- as.numeric(ks.test(snd_2020_2_INTERACT_sites, sample(snd_2020_2_INTERACT_domain$Value, 496))[c(1:2)])
  snd_2_ks_2020_war <- as.numeric(ks.test(snd_2020_2_INTERACT_sites_war, sample(snd_2020_2_INTERACT_domain$Value, 496))[c(1:2)])
  snd_2_ks_2100 <- as.numeric(ks.test(snd_2100_2_INTERACT_sites, sample(snd_2100_2_INTERACT_domain$Value, 496))[c(1:2)])
  snd_2_ks_2020_2100 <- as.numeric(ks.test(sample(snd_2020_2_INTERACT_domain$Value, 496), sample(snd_2100_2_INTERACT_domain$Value, 496))[c(1:2)])
  #snd model_2_ INTERACT medians
  snd_2_I_25_2020 <- as.numeric(quantile(snd_2020_2_INTERACT_sites,probs = c(0.25),na.rm = T))
  snd_2_I_50_2020 <- as.numeric(median(snd_2020_2_INTERACT_sites,na.rm = T))
  snd_2_I_75_2020 <- as.numeric(quantile(snd_2020_2_INTERACT_sites,probs = c(0.75),na.rm = T))
  #snd model_2_ INTERACT (without Russia) medians
  snd_2_W_25_2020 <- as.numeric(quantile(snd_2020_2_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  snd_2_W_50_2020 <- as.numeric(median(snd_2020_2_INTERACT_sites_war,na.rm = T))
  snd_2_W_75_2020 <- as.numeric(quantile(snd_2020_2_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #snd model_2_ Domain medians & 50% CI 2020
  snd_2_D_25_2020 <- as.numeric(quantile(sample(snd_2020_2_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  snd_2_D_50_2020 <- as.numeric(median(sample(snd_2020_2_INTERACT_domain$Value, 496),na.rm = T))
  snd_2_D_75_2020 <- as.numeric(quantile(sample(snd_2020_2_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #snd model_2_ Domain medians & 50% CI 2100
  snd_2_D_25_2100 <- as.numeric(quantile(sample(snd_2100_2_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  snd_2_D_50_2100 <- as.numeric(median(sample(snd_2100_2_INTERACT_domain$Value, 496),na.rm = T))
  snd_2_D_75_2100 <- as.numeric(quantile(sample(snd_2100_2_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  snd_2[i,] <- round(c(snd_2_ks_2020[1],snd_2_ks_2020[2],snd_2_ks_2020_war[1],snd_2_ks_2020_war[2],
                       snd_2_ks_2100[1],snd_2_ks_2100[2],snd_2_ks_2020_2100[1],snd_2_ks_2020_2100[2],
                       snd_2_I_25_2020,snd_2_I_50_2020,snd_2_I_75_2020,snd_2_W_25_2020,snd_2_W_50_2020,snd_2_W_75_2020,
                       snd_2_D_25_2020,snd_2_D_50_2020,snd_2_D_75_2020,snd_2_D_25_2100,snd_2_D_50_2100,snd_2_D_75_2100),4)
  ##############################################################################################################################
  #snd model_3_ KS test stats 
  snd_3_ks_2020 <- as.numeric(ks.test(snd_2020_3_INTERACT_sites, sample(snd_2020_3_INTERACT_domain$Value, 496))[c(1:2)])
  snd_3_ks_2020_war <- as.numeric(ks.test(snd_2020_3_INTERACT_sites_war, sample(snd_2020_3_INTERACT_domain$Value, 496))[c(1:2)])
  snd_3_ks_2100 <- as.numeric(ks.test(snd_2100_3_INTERACT_sites, sample(snd_2100_3_INTERACT_domain$Value, 496))[c(1:2)])
  snd_3_ks_2020_2100 <- as.numeric(ks.test(sample(snd_2020_3_INTERACT_domain$Value, 496), sample(snd_2100_3_INTERACT_domain$Value, 496))[c(1:2)])
  #snd model_3_ INTERACT medians
  snd_3_I_25_2020 <- as.numeric(quantile(snd_2020_3_INTERACT_sites,probs = c(0.25),na.rm = T))
  snd_3_I_50_2020 <- as.numeric(median(snd_2020_3_INTERACT_sites,na.rm = T))
  snd_3_I_75_2020 <- as.numeric(quantile(snd_2020_3_INTERACT_sites,probs = c(0.75),na.rm = T))
  #snd model_3_ INTERACT (without Russia) medians
  snd_3_W_25_2020 <- as.numeric(quantile(snd_2020_3_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  snd_3_W_50_2020 <- as.numeric(median(snd_2020_3_INTERACT_sites_war,na.rm = T))
  snd_3_W_75_2020 <- as.numeric(quantile(snd_2020_3_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #snd model_3_ Domain medians & 50% CI 2020
  snd_3_D_25_2020 <- as.numeric(quantile(sample(snd_2020_3_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  snd_3_D_50_2020 <- as.numeric(median(sample(snd_2020_3_INTERACT_domain$Value, 496),na.rm = T))
  snd_3_D_75_2020 <- as.numeric(quantile(sample(snd_2020_3_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #snd model_3_ Domain medians & 50% CI 2100
  snd_3_D_25_2100 <- as.numeric(quantile(sample(snd_2100_3_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  snd_3_D_50_2100 <- as.numeric(median(sample(snd_2100_3_INTERACT_domain$Value, 496),na.rm = T))
  snd_3_D_75_2100 <- as.numeric(quantile(sample(snd_2100_3_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  snd_3[i,] <- round(c(snd_3_ks_2020[1],snd_3_ks_2020[2],snd_3_ks_2020_war[1],snd_3_ks_2020_war[2],
                       snd_3_ks_2100[1],snd_3_ks_2100[2],snd_3_ks_2020_2100[1],snd_3_ks_2020_2100[2],
                       snd_3_I_25_2020,snd_3_I_50_2020,snd_3_I_75_2020,snd_3_W_25_2020,snd_3_W_50_2020,snd_3_W_75_2020,
                       snd_3_D_25_2020,snd_3_D_50_2020,snd_3_D_75_2020,snd_3_D_25_2100,snd_3_D_50_2100,snd_3_D_75_2100),4)
  ##############################################################################################################################
  #snd model_4_ KS test stats 
  #snd_4_ks_2020 <- as.numeric(ks.test(snd_2020_4_INTERACT_sites, sample(snd_2020_4_INTERACT_domain$Value, 496))[c(1:2)])
  #snd_4_ks_2020_war <- as.numeric(ks.test(snd_2020_4_INTERACT_sites_war, sample(snd_2020_4_INTERACT_domain$Value, 496))[c(1:2)])
  #snd_4_ks_2100 <- as.numeric(ks.test(snd_2100_4_INTERACT_sites, sample(snd_2100_4_INTERACT_domain$Value, 496))[c(1:2)])
  #snd_4_ks_2020_2100 <- as.numeric(ks.test(sample(snd_2020_4_INTERACT_domain$Value, 496), sample(snd_2100_4_INTERACT_domain$Value, 496))[c(1:2)])
  #snd model_4_ INTERACT medians
  #snd_4_I_25_2020 <- as.numeric(quantile(snd_2020_4_INTERACT_sites,probs = c(0.25),na.rm = T))
  #snd_4_I_50_2020 <- as.numeric(median(snd_2020_4_INTERACT_sites,na.rm = T))
  #snd_4_I_75_2020 <- as.numeric(quantile(snd_2020_4_INTERACT_sites,probs = c(0.75),na.rm = T))
  #snd model_4_ INTERACT (without Russia) medians
  #snd_4_W_25_2020 <- as.numeric(quantile(snd_2020_4_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  #snd_4_W_50_2020 <- as.numeric(median(snd_2020_4_INTERACT_sites_war,na.rm = T))
  #snd_4_W_75_2020 <- as.numeric(quantile(snd_2020_4_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #snd model_4_ Domain medians & 50% CI 2020
  #snd_4_D_25_2020 <- as.numeric(quantile(sample(snd_2020_4_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  #snd_4_D_50_2020 <- as.numeric(median(sample(snd_2020_4_INTERACT_domain$Value, 496),na.rm = T))
  ##snd_4_D_75_2020 <- as.numeric(quantile(sample(snd_2020_4_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #snd model_4_ Domain medians & 50% CI 2100
  #snd_4_D_25_2100 <- as.numeric(quantile(sample(snd_2100_4_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  #snd_4_D_50_2100 <- as.numeric(median(sample(snd_2100_4_INTERACT_domain$Value, 496),na.rm = T))
  #snd_4_D_75_2100 <- as.numeric(quantile(sample(snd_2100_4_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  #snd_4[i,] <- round(c(snd_4_ks_2020[1],snd_4_ks_2020[2],snd_4_ks_2020_war[1],snd_4_ks_2020_war[2],
  #                     snd_4_ks_2100[1],snd_4_ks_2100[2],snd_4_ks_2020_2100[1],snd_4_ks_2020_2100[2],
  #                     snd_4_I_25_2020,snd_4_I_50_2020,snd_4_I_75_2020,snd_4_W_25_2020,snd_4_W_50_2020,snd_4_W_75_2020,
  #                     snd_4_D_25_2020,snd_4_D_50_2020,snd_4_D_75_2020,snd_4_D_25_2100,snd_4_D_50_2100,snd_4_D_75_2100),4)
  # ##############################################################################################################################
  # #snd model_5_ KS test stats
  # #snd_5_ks_2020 <- as.numeric(ks.test(snd_2020_5_INTERACT_sites, sample(snd_2020_5_INTERACT_domain$Value, 496))[c(1:2)])
  # #snd_5_ks_2020_war <- as.numeric(ks.test(snd_2020_5_INTERACT_sites_war, sample(snd_2020_5_INTERACT_domain$Value, 496))[c(1:2)])
  # #snd_5_ks_2100 <- as.numeric(ks.test(snd_2100_5_INTERACT_sites, sample(snd_2100_5_INTERACT_domain$Value, 496))[c(1:2)])
  # #snd_5_ks_2020_2100 <- as.numeric(ks.test(sample(snd_2020_5_INTERACT_domain$Value, 496), sample(snd_2100_5_INTERACT_domain$Value, 496))[c(1:2)])
  # #snd model_5_ INTERACT medians
  # snd_5_I_25_2020 <- as.numeric(quantile(snd_2020_5_INTERACT_sites,probs = c(0.25),na.rm = T))
  # snd_5_I_50_2020 <- as.numeric(median(snd_2020_5_INTERACT_sites,na.rm = T))
  # snd_5_I_75_2020 <- as.numeric(quantile(snd_2020_5_INTERACT_sites,probs = c(0.75),na.rm = T))
  # #snd model_5_ INTERACT (without Russia) medians
  # snd_5_W_25_2020 <- as.numeric(quantile(snd_2020_5_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  # snd_5_W_50_2020 <- as.numeric(median(snd_2020_5_INTERACT_sites_war,na.rm = T))
  # snd_5_W_75_2020 <- as.numeric(quantile(snd_2020_5_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  # #snd model_5_ Domain medians & 50% CI 2020
  # snd_5_D_25_2020 <- as.numeric(quantile(sample(snd_2020_5_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  # snd_5_D_50_2020 <- as.numeric(median(sample(snd_2020_5_INTERACT_domain$Value, 496),na.rm = T))
  # snd_5_D_75_2020 <- as.numeric(quantile(sample(snd_2020_5_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  # #snd model_5_ Domain medians & 50% CI 2100
  # snd_5_D_25_2100 <- as.numeric(quantile(sample(snd_2100_5_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  # snd_5_D_50_2100 <- as.numeric(median(sample(snd_2100_5_INTERACT_domain$Value, 496),na.rm = T))
  # snd_5_D_75_2100 <- as.numeric(quantile(sample(snd_2100_5_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  # #add all to same df
  # snd_5[i,] <- round(c(snd_5_ks_2020[1],snd_5_ks_2020[2],snd_5_ks_2020_war[1],snd_5_ks_2020_war[2],
  #                      snd_5_ks_2100[1],snd_5_ks_2100[2],snd_5_ks_2020_2100[1],snd_5_ks_2020_2100[2],
  #                      snd_5_I_25_2020,snd_5_I_50_2020,snd_5_I_75_2020,snd_5_W_25_2020,snd_5_W_50_2020,snd_5_W_75_2020,
  #                      snd_5_D_25_2020,snd_5_D_50_2020,snd_5_D_75_2020,snd_5_D_25_2100,snd_5_D_50_2100,snd_5_D_75_2100),4)
  # ##############################################################################################################################
  #snd model_6_ KS test stats 
  snd_6_ks_2020 <- as.numeric(ks.test(snd_2020_6_INTERACT_sites, sample(snd_2020_6_INTERACT_domain$Value, 496))[c(1:2)])
  snd_6_ks_2020_war <- as.numeric(ks.test(snd_2020_6_INTERACT_sites_war, sample(snd_2020_6_INTERACT_domain$Value, 496))[c(1:2)])
  snd_6_ks_2100 <- as.numeric(ks.test(snd_2100_6_INTERACT_sites, sample(snd_2100_6_INTERACT_domain$Value, 496))[c(1:2)])
  snd_6_ks_2020_2100 <- as.numeric(ks.test(sample(snd_2020_6_INTERACT_domain$Value, 496), sample(snd_2100_6_INTERACT_domain$Value, 496))[c(1:2)])
  #snd model_6_ INTERACT medians
  snd_6_I_25_2020 <- as.numeric(quantile(snd_2020_6_INTERACT_sites,probs = c(0.25),na.rm = T))
  snd_6_I_50_2020 <- as.numeric(median(snd_2020_6_INTERACT_sites,na.rm = T))
  snd_6_I_75_2020 <- as.numeric(quantile(snd_2020_6_INTERACT_sites,probs = c(0.75),na.rm = T))
  #snd model_6_ INTERACT (without Russia) medians
  snd_6_W_25_2020 <- as.numeric(quantile(snd_2020_6_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  snd_6_W_50_2020 <- as.numeric(median(snd_2020_6_INTERACT_sites_war,na.rm = T))
  snd_6_W_75_2020 <- as.numeric(quantile(snd_2020_6_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #snd model_6_ Domain medians & 50% CI 2020
  snd_6_D_25_2020 <- as.numeric(quantile(sample(snd_2020_6_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  snd_6_D_50_2020 <- as.numeric(median(sample(snd_2020_6_INTERACT_domain$Value, 496),na.rm = T))
  snd_6_D_75_2020 <- as.numeric(quantile(sample(snd_2020_6_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #snd model_6_ Domain medians & 50% CI 2100
  snd_6_D_25_2100 <- as.numeric(quantile(sample(snd_2100_6_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  snd_6_D_50_2100 <- as.numeric(median(sample(snd_2100_6_INTERACT_domain$Value, 496),na.rm = T))
  snd_6_D_75_2100 <- as.numeric(quantile(sample(snd_2100_6_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  snd_6[i,] <- round(c(snd_6_ks_2020[1],snd_6_ks_2020[2],snd_6_ks_2020_war[1],snd_6_ks_2020_war[2],
                       snd_6_ks_2100[1],snd_6_ks_2100[2],snd_6_ks_2020_2100[1],snd_6_ks_2020_2100[2],
                       snd_6_I_25_2020,snd_6_I_50_2020,snd_6_I_75_2020,snd_6_W_25_2020,snd_6_W_50_2020,snd_6_W_75_2020,
                       snd_6_D_25_2020,snd_6_D_50_2020,snd_6_D_75_2020,snd_6_D_25_2100,snd_6_D_50_2100,snd_6_D_75_2100),4)
  ##############################################################################################################################
  #snd model_7_ KS test stats 
  snd_7_ks_2020 <- as.numeric(ks.test(snd_2020_7_INTERACT_sites, sample(snd_2020_7_INTERACT_domain$Value, 496))[c(1:2)])
  snd_7_ks_2020_war <- as.numeric(ks.test(snd_2020_7_INTERACT_sites_war, sample(snd_2020_7_INTERACT_domain$Value, 496))[c(1:2)])
  snd_7_ks_2100 <- as.numeric(ks.test(snd_2100_7_INTERACT_sites, sample(snd_2100_7_INTERACT_domain$Value, 496))[c(1:2)])
  snd_7_ks_2020_2100 <- as.numeric(ks.test(sample(snd_2020_7_INTERACT_domain$Value, 496), sample(snd_2100_7_INTERACT_domain$Value, 496))[c(1:2)])
  #snd model_7_ INTERACT medians
  snd_7_I_25_2020 <- as.numeric(quantile(snd_2020_7_INTERACT_sites,probs = c(0.25),na.rm = T))
  snd_7_I_50_2020 <- as.numeric(median(snd_2020_7_INTERACT_sites,na.rm = T))
  snd_7_I_75_2020 <- as.numeric(quantile(snd_2020_7_INTERACT_sites,probs = c(0.75),na.rm = T))
  #snd model_7_ INTERACT (without Russia) medians
  snd_7_W_25_2020 <- as.numeric(quantile(snd_2020_7_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  snd_7_W_50_2020 <- as.numeric(median(snd_2020_7_INTERACT_sites_war,na.rm = T))
  snd_7_W_75_2020 <- as.numeric(quantile(snd_2020_7_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #snd model_7_ Domain medians & 50% CI 2020
  snd_7_D_25_2020 <- as.numeric(quantile(sample(snd_2020_7_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  snd_7_D_50_2020 <- as.numeric(median(sample(snd_2020_7_INTERACT_domain$Value, 496),na.rm = T))
  snd_7_D_75_2020 <- as.numeric(quantile(sample(snd_2020_7_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #snd model_7_ Domain medians & 50% CI 2100
  snd_7_D_25_2100 <- as.numeric(quantile(sample(snd_2100_7_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  snd_7_D_50_2100 <- as.numeric(median(sample(snd_2100_7_INTERACT_domain$Value, 496),na.rm = T))
  snd_7_D_75_2100 <- as.numeric(quantile(sample(snd_2100_7_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  snd_7[i,] <- round(c(snd_7_ks_2020[1],snd_7_ks_2020[2],snd_7_ks_2020_war[1],snd_7_ks_2020_war[2],
                       snd_7_ks_2100[1],snd_7_ks_2100[2],snd_7_ks_2020_2100[1],snd_7_ks_2020_2100[2],
                       snd_7_I_25_2020,snd_7_I_50_2020,snd_7_I_75_2020,snd_7_W_25_2020,snd_7_W_50_2020,snd_7_W_75_2020,
                       snd_7_D_25_2020,snd_7_D_50_2020,snd_7_D_75_2020,snd_7_D_25_2100,snd_7_D_50_2100,snd_7_D_75_2100),4)
  
  ##################################################################################################################################
  #snd model_8_ KS test stats 
  snd_8_ks_2020 <- as.numeric(ks.test(snd_2020_8_INTERACT_sites, sample(snd_2020_8_INTERACT_domain$Value, 496))[c(1:2)])
  snd_8_ks_2020_war <- as.numeric(ks.test(snd_2020_8_INTERACT_sites_war, sample(snd_2020_8_INTERACT_domain$Value, 496))[c(1:2)])
  snd_8_ks_2100 <- as.numeric(ks.test(snd_2100_8_INTERACT_sites, sample(snd_2100_8_INTERACT_domain$Value, 496))[c(1:2)])
  snd_8_ks_2020_2100 <- as.numeric(ks.test(sample(snd_2020_8_INTERACT_domain$Value, 496), sample(snd_2100_8_INTERACT_domain$Value, 496))[c(1:2)])
  #snd model_8_ INTERACT medians
  snd_8_I_25_2020 <- as.numeric(quantile(snd_2020_8_INTERACT_sites,probs = c(0.25),na.rm = T))
  snd_8_I_50_2020 <- as.numeric(median(snd_2020_8_INTERACT_sites,na.rm = T))
  snd_8_I_75_2020 <- as.numeric(quantile(snd_2020_8_INTERACT_sites,probs = c(0.75),na.rm = T))
  #snd model_8_ INTERACT (without Russia) medians
  snd_8_W_25_2020 <- as.numeric(quantile(snd_2020_8_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  snd_8_W_50_2020 <- as.numeric(median(snd_2020_8_INTERACT_sites_war,na.rm = T))
  snd_8_W_75_2020 <- as.numeric(quantile(snd_2020_8_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #snd model_8_ Domain medians & 50% CI 2020
  snd_8_D_25_2020 <- as.numeric(quantile(sample(snd_2020_8_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  snd_8_D_50_2020 <- as.numeric(median(sample(snd_2020_8_INTERACT_domain$Value, 496),na.rm = T))
  snd_8_D_75_2020 <- as.numeric(quantile(sample(snd_2020_8_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #snd model_8_ Domain medians & 50% CI 2100
  snd_8_D_25_2100 <- as.numeric(quantile(sample(snd_2100_8_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  snd_8_D_50_2100 <- as.numeric(median(sample(snd_2100_8_INTERACT_domain$Value, 496),na.rm = T))
  snd_8_D_75_2100 <- as.numeric(quantile(sample(snd_2100_8_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  snd_8[i,] <- round(c(snd_8_ks_2020[1],snd_8_ks_2020[2],snd_8_ks_2020_war[1],snd_8_ks_2020_war[2],
                       snd_8_ks_2100[1],snd_8_ks_2100[2],snd_8_ks_2020_2100[1],snd_8_ks_2020_2100[2],
                       snd_8_I_25_2020,snd_8_I_50_2020,snd_8_I_75_2020,snd_8_W_25_2020,snd_8_W_50_2020,snd_8_W_75_2020,
                       snd_8_D_25_2020,snd_8_D_50_2020,snd_8_D_75_2020,snd_8_D_25_2100,snd_8_D_50_2100,snd_8_D_75_2100),4)
  ##################################################################################################################################
  #mrsos model_1_ KS test stats 
  mrsos_1_ks_2020 <- as.numeric(ks.test(mrsos_2020_1_INTERACT_sites, sample(mrsos_2020_1_INTERACT_domain$Value, 496))[c(1:2)])
  mrsos_1_ks_2020_war <- as.numeric(ks.test(mrsos_2020_1_INTERACT_sites_war, sample(mrsos_2020_1_INTERACT_domain$Value, 496))[c(1:2)])
  mrsos_1_ks_2100 <- as.numeric(ks.test(mrsos_2100_1_INTERACT_sites, sample(mrsos_2100_1_INTERACT_domain$Value, 496))[c(1:2)])
  mrsos_1_ks_2020_2100 <- as.numeric(ks.test(sample(mrsos_2020_1_INTERACT_domain$Value, 496), sample(mrsos_2100_1_INTERACT_domain$Value, 496))[c(1:2)])
  #mrsos model_1_ INTERACT medians
  mrsos_1_I_25_2020 <- as.numeric(quantile(mrsos_2020_1_INTERACT_sites,probs = c(0.25),na.rm = T))
  mrsos_1_I_50_2020 <- as.numeric(median(mrsos_2020_1_INTERACT_sites,na.rm = T))
  mrsos_1_I_75_2020 <- as.numeric(quantile(mrsos_2020_1_INTERACT_sites,probs = c(0.75),na.rm = T))
  #mrsos model_1_ INTERACT (without Russia) medians
  mrsos_1_W_25_2020 <- as.numeric(quantile(mrsos_2020_1_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  mrsos_1_W_50_2020 <- as.numeric(median(mrsos_2020_1_INTERACT_sites_war,na.rm = T))
  mrsos_1_W_75_2020 <- as.numeric(quantile(mrsos_2020_1_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #mrsos model_1_ Domain medians & 50% CI 2020
  mrsos_1_D_25_2020 <- as.numeric(quantile(sample(mrsos_2020_1_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  mrsos_1_D_50_2020 <- as.numeric(median(sample(mrsos_2020_1_INTERACT_domain$Value, 496),na.rm = T))
  mrsos_1_D_75_2020 <- as.numeric(quantile(sample(mrsos_2020_1_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #mrsos model_1_ Domain medians & 50% CI 2100
  mrsos_1_D_25_2100 <- as.numeric(quantile(sample(mrsos_2100_1_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  mrsos_1_D_50_2100 <- as.numeric(median(sample(mrsos_2100_1_INTERACT_domain$Value, 496),na.rm = T))
  mrsos_1_D_75_2100 <- as.numeric(quantile(sample(mrsos_2100_1_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  mrsos_1[i,] <- round(c(mrsos_1_ks_2020[1],mrsos_1_ks_2020[2],mrsos_1_ks_2020_war[1],mrsos_1_ks_2020_war[2],
                         mrsos_1_ks_2100[1],mrsos_1_ks_2100[2],mrsos_1_ks_2020_2100[1],mrsos_1_ks_2020_2100[2],
                         mrsos_1_I_25_2020,mrsos_1_I_50_2020,mrsos_1_I_75_2020,mrsos_1_W_25_2020,mrsos_1_W_50_2020,mrsos_1_W_75_2020,
                         mrsos_1_D_25_2020,mrsos_1_D_50_2020,mrsos_1_D_75_2020,mrsos_1_D_25_2100,mrsos_1_D_50_2100,mrsos_1_D_75_2100),4)
  ##################################################################################################################################
  #mrsos model_2_ KS test stats 
  mrsos_2_ks_2020 <- as.numeric(ks.test(mrsos_2020_2_INTERACT_sites, sample(mrsos_2020_2_INTERACT_domain$Value, 496))[c(1:2)])
  mrsos_2_ks_2020_war <- as.numeric(ks.test(mrsos_2020_2_INTERACT_sites_war, sample(mrsos_2020_2_INTERACT_domain$Value, 496))[c(1:2)])
  mrsos_2_ks_2100 <- as.numeric(ks.test(mrsos_2100_2_INTERACT_sites, sample(mrsos_2100_2_INTERACT_domain$Value, 496))[c(1:2)])
  mrsos_2_ks_2020_2100 <- as.numeric(ks.test(sample(mrsos_2020_2_INTERACT_domain$Value, 496), sample(mrsos_2100_2_INTERACT_domain$Value, 496))[c(1:2)])
  #mrsos model_2_ INTERACT medians
  mrsos_2_I_25_2020 <- as.numeric(quantile(mrsos_2020_2_INTERACT_sites,probs = c(0.25),na.rm = T))
  mrsos_2_I_50_2020 <- as.numeric(median(mrsos_2020_2_INTERACT_sites,na.rm = T))
  mrsos_2_I_75_2020 <- as.numeric(quantile(mrsos_2020_2_INTERACT_sites,probs = c(0.75),na.rm = T))
  #mrsos model_2_ INTERACT (without Russia) medians
  mrsos_2_W_25_2020 <- as.numeric(quantile(mrsos_2020_2_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  mrsos_2_W_50_2020 <- as.numeric(median(mrsos_2020_2_INTERACT_sites_war,na.rm = T))
  mrsos_2_W_75_2020 <- as.numeric(quantile(mrsos_2020_2_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #mrsos model_2_ Domain medians & 50% CI 2020
  mrsos_2_D_25_2020 <- as.numeric(quantile(sample(mrsos_2020_2_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  mrsos_2_D_50_2020 <- as.numeric(median(sample(mrsos_2020_2_INTERACT_domain$Value, 496),na.rm = T))
  mrsos_2_D_75_2020 <- as.numeric(quantile(sample(mrsos_2020_2_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #mrsos model_2_ Domain medians & 50% CI 2100
  mrsos_2_D_25_2100 <- as.numeric(quantile(sample(mrsos_2100_2_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  mrsos_2_D_50_2100 <- as.numeric(median(sample(mrsos_2100_2_INTERACT_domain$Value, 496),na.rm = T))
  mrsos_2_D_75_2100 <- as.numeric(quantile(sample(mrsos_2100_2_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  mrsos_2[i,] <- round(c(mrsos_2_ks_2020[1],mrsos_2_ks_2020[2],mrsos_2_ks_2020_war[1],mrsos_2_ks_2020_war[2],
                         mrsos_2_ks_2100[1],mrsos_2_ks_2100[2],mrsos_2_ks_2020_2100[1],mrsos_2_ks_2020_2100[2],
                         mrsos_2_I_25_2020,mrsos_2_I_50_2020,mrsos_2_I_75_2020,mrsos_2_W_25_2020,mrsos_2_W_50_2020,mrsos_2_W_75_2020,
                         mrsos_2_D_25_2020,mrsos_2_D_50_2020,mrsos_2_D_75_2020,mrsos_2_D_25_2100,mrsos_2_D_50_2100,mrsos_2_D_75_2100),4)
  ##############################################################################################################################
  #mrsos model_3_ KS test stats 
  mrsos_3_ks_2020 <- as.numeric(ks.test(mrsos_2020_3_INTERACT_sites, sample(mrsos_2020_3_INTERACT_domain$Value, 496))[c(1:2)])
  mrsos_3_ks_2020_war <- as.numeric(ks.test(mrsos_2020_3_INTERACT_sites_war, sample(mrsos_2020_3_INTERACT_domain$Value, 496))[c(1:2)])
  mrsos_3_ks_2100 <- as.numeric(ks.test(mrsos_2100_3_INTERACT_sites, sample(mrsos_2100_3_INTERACT_domain$Value, 496))[c(1:2)])
  mrsos_3_ks_2020_2100 <- as.numeric(ks.test(sample(mrsos_2020_3_INTERACT_domain$Value, 496), sample(mrsos_2100_3_INTERACT_domain$Value, 496))[c(1:2)])
  #mrsos model_3_ INTERACT medians
  mrsos_3_I_25_2020 <- as.numeric(quantile(mrsos_2020_3_INTERACT_sites,probs = c(0.25),na.rm = T))
  mrsos_3_I_50_2020 <- as.numeric(median(mrsos_2020_3_INTERACT_sites,na.rm = T))
  mrsos_3_I_75_2020 <- as.numeric(quantile(mrsos_2020_3_INTERACT_sites,probs = c(0.75),na.rm = T))
  #mrsos model_3_ INTERACT (without Russia) medians
  mrsos_3_W_25_2020 <- as.numeric(quantile(mrsos_2020_3_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  mrsos_3_W_50_2020 <- as.numeric(median(mrsos_2020_3_INTERACT_sites_war,na.rm = T))
  mrsos_3_W_75_2020 <- as.numeric(quantile(mrsos_2020_3_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #mrsos model_3_ Domain medians & 50% CI 2020
  mrsos_3_D_25_2020 <- as.numeric(quantile(sample(mrsos_2020_3_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  mrsos_3_D_50_2020 <- as.numeric(median(sample(mrsos_2020_3_INTERACT_domain$Value, 496),na.rm = T))
  mrsos_3_D_75_2020 <- as.numeric(quantile(sample(mrsos_2020_3_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #mrsos model_3_ Domain medians & 50% CI 2100
  mrsos_3_D_25_2100 <- as.numeric(quantile(sample(mrsos_2100_3_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  mrsos_3_D_50_2100 <- as.numeric(median(sample(mrsos_2100_3_INTERACT_domain$Value, 496),na.rm = T))
  mrsos_3_D_75_2100 <- as.numeric(quantile(sample(mrsos_2100_3_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  mrsos_3[i,] <- round(c(mrsos_3_ks_2020[1],mrsos_3_ks_2020[2],mrsos_3_ks_2020_war[1],mrsos_3_ks_2020_war[2],
                         mrsos_3_ks_2100[1],mrsos_3_ks_2100[2],mrsos_3_ks_2020_2100[1],mrsos_3_ks_2020_2100[2],
                         mrsos_3_I_25_2020,mrsos_3_I_50_2020,mrsos_3_I_75_2020,mrsos_3_W_25_2020,mrsos_3_W_50_2020,mrsos_3_W_75_2020,
                         mrsos_3_D_25_2020,mrsos_3_D_50_2020,mrsos_3_D_75_2020,mrsos_3_D_25_2100,mrsos_3_D_50_2100,mrsos_3_D_75_2100),4)
  ##############################################################################################################################
  #mrsos model_4_ KS test stats 
  mrsos_4_ks_2020 <- as.numeric(ks.test(mrsos_2020_4_INTERACT_sites, sample(mrsos_2020_4_INTERACT_domain$Value, 496))[c(1:2)])
  mrsos_4_ks_2020_war <- as.numeric(ks.test(mrsos_2020_4_INTERACT_sites_war, sample(mrsos_2020_4_INTERACT_domain$Value, 496))[c(1:2)])
  mrsos_4_ks_2100 <- as.numeric(ks.test(mrsos_2100_4_INTERACT_sites, sample(mrsos_2100_4_INTERACT_domain$Value, 496))[c(1:2)])
  mrsos_4_ks_2020_2100 <- as.numeric(ks.test(sample(mrsos_2020_4_INTERACT_domain$Value, 496), sample(mrsos_2100_4_INTERACT_domain$Value, 496))[c(1:2)])
  #mrsos model_4_ INTERACT medians
  mrsos_4_I_25_2020 <- as.numeric(quantile(mrsos_2020_4_INTERACT_sites,probs = c(0.25),na.rm = T))
  mrsos_4_I_50_2020 <- as.numeric(median(mrsos_2020_4_INTERACT_sites,na.rm = T))
  mrsos_4_I_75_2020 <- as.numeric(quantile(mrsos_2020_4_INTERACT_sites,probs = c(0.75),na.rm = T))
  #mrsos model_4_ INTERACT (without Russia) medians
  mrsos_4_W_25_2020 <- as.numeric(quantile(mrsos_2020_4_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  mrsos_4_W_50_2020 <- as.numeric(median(mrsos_2020_4_INTERACT_sites_war,na.rm = T))
  mrsos_4_W_75_2020 <- as.numeric(quantile(mrsos_2020_4_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #mrsos model_4_ Domain medians & 50% CI 2020
  mrsos_4_D_25_2020 <- as.numeric(quantile(sample(mrsos_2020_4_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  mrsos_4_D_50_2020 <- as.numeric(median(sample(mrsos_2020_4_INTERACT_domain$Value, 496),na.rm = T))
  mrsos_4_D_75_2020 <- as.numeric(quantile(sample(mrsos_2020_4_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #mrsos model_4_ Domain medians & 50% CI 2100
  mrsos_4_D_25_2100 <- as.numeric(quantile(sample(mrsos_2100_4_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  mrsos_4_D_50_2100 <- as.numeric(median(sample(mrsos_2100_4_INTERACT_domain$Value, 496),na.rm = T))
  mrsos_4_D_75_2100 <- as.numeric(quantile(sample(mrsos_2100_4_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  mrsos_4[i,] <- round(c(mrsos_4_ks_2020[1],mrsos_4_ks_2020[2],mrsos_4_ks_2020_war[1],mrsos_4_ks_2020_war[2],
                         mrsos_4_ks_2100[1],mrsos_4_ks_2100[2],mrsos_4_ks_2020_2100[1],mrsos_4_ks_2020_2100[2],
                         mrsos_4_I_25_2020,mrsos_4_I_50_2020,mrsos_4_I_75_2020,mrsos_4_W_25_2020,mrsos_4_W_50_2020,mrsos_4_W_75_2020,
                         mrsos_4_D_25_2020,mrsos_4_D_50_2020,mrsos_4_D_75_2020,mrsos_4_D_25_2100,mrsos_4_D_50_2100,mrsos_4_D_75_2100),4)
  ##############################################################################################################################
  #mrsos model_5_ KS test stats 
  mrsos_5_ks_2020 <- as.numeric(ks.test(mrsos_2020_5_INTERACT_sites, sample(mrsos_2020_5_INTERACT_domain$Value, 496))[c(1:2)])
  mrsos_5_ks_2020_war <- as.numeric(ks.test(mrsos_2020_5_INTERACT_sites_war, sample(mrsos_2020_5_INTERACT_domain$Value, 496))[c(1:2)])
  mrsos_5_ks_2100 <- as.numeric(ks.test(mrsos_2100_5_INTERACT_sites, sample(mrsos_2100_5_INTERACT_domain$Value, 496))[c(1:2)])
  mrsos_5_ks_2020_2100 <- as.numeric(ks.test(sample(mrsos_2020_5_INTERACT_domain$Value, 496), sample(mrsos_2100_5_INTERACT_domain$Value, 496))[c(1:2)])
  #mrsos model_5_ INTERACT medians
  mrsos_5_I_25_2020 <- as.numeric(quantile(mrsos_2020_5_INTERACT_sites,probs = c(0.25),na.rm = T))
  mrsos_5_I_50_2020 <- as.numeric(median(mrsos_2020_5_INTERACT_sites,na.rm = T))
  mrsos_5_I_75_2020 <- as.numeric(quantile(mrsos_2020_5_INTERACT_sites,probs = c(0.75),na.rm = T))
  #mrsos model_5_ INTERACT (without Russia) medians
  mrsos_5_W_25_2020 <- as.numeric(quantile(mrsos_2020_5_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  mrsos_5_W_50_2020 <- as.numeric(median(mrsos_2020_5_INTERACT_sites_war,na.rm = T))
  mrsos_5_W_75_2020 <- as.numeric(quantile(mrsos_2020_5_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #mrsos model_5_ Domain medians & 50% CI 2020
  mrsos_5_D_25_2020 <- as.numeric(quantile(sample(mrsos_2020_5_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  mrsos_5_D_50_2020 <- as.numeric(median(sample(mrsos_2020_5_INTERACT_domain$Value, 496),na.rm = T))
  mrsos_5_D_75_2020 <- as.numeric(quantile(sample(mrsos_2020_5_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #mrsos model_5_ Domain medians & 50% CI 2100
  mrsos_5_D_25_2100 <- as.numeric(quantile(sample(mrsos_2100_5_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  mrsos_5_D_50_2100 <- as.numeric(median(sample(mrsos_2100_5_INTERACT_domain$Value, 496),na.rm = T))
  mrsos_5_D_75_2100 <- as.numeric(quantile(sample(mrsos_2100_5_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  mrsos_5[i,] <- round(c(mrsos_5_ks_2020[1],mrsos_5_ks_2020[2],mrsos_5_ks_2020_war[1],mrsos_5_ks_2020_war[2],
                         mrsos_5_ks_2100[1],mrsos_5_ks_2100[2],mrsos_5_ks_2020_2100[1],mrsos_5_ks_2020_2100[2],
                         mrsos_5_I_25_2020,mrsos_5_I_50_2020,mrsos_5_I_75_2020,mrsos_5_W_25_2020,mrsos_5_W_50_2020,mrsos_5_W_75_2020,
                         mrsos_5_D_25_2020,mrsos_5_D_50_2020,mrsos_5_D_75_2020,mrsos_5_D_25_2100,mrsos_5_D_50_2100,mrsos_5_D_75_2100),4)
  ##############################################################################################################################
  #mrsos model_6_ KS test stats 
  mrsos_6_ks_2020 <- as.numeric(ks.test(mrsos_2020_6_INTERACT_sites, sample(mrsos_2020_6_INTERACT_domain$Value, 496))[c(1:2)])
  mrsos_6_ks_2020_war <- as.numeric(ks.test(mrsos_2020_6_INTERACT_sites_war, sample(mrsos_2020_6_INTERACT_domain$Value, 496))[c(1:2)])
  mrsos_6_ks_2100 <- as.numeric(ks.test(mrsos_2100_6_INTERACT_sites, sample(mrsos_2100_6_INTERACT_domain$Value, 496))[c(1:2)])
  mrsos_6_ks_2020_2100 <- as.numeric(ks.test(sample(mrsos_2020_6_INTERACT_domain$Value, 496), sample(mrsos_2100_6_INTERACT_domain$Value, 496))[c(1:2)])
  #mrsos model_6_ INTERACT medians
  mrsos_6_I_25_2020 <- as.numeric(quantile(mrsos_2020_6_INTERACT_sites,probs = c(0.25),na.rm = T))
  mrsos_6_I_50_2020 <- as.numeric(median(mrsos_2020_6_INTERACT_sites,na.rm = T))
  mrsos_6_I_75_2020 <- as.numeric(quantile(mrsos_2020_6_INTERACT_sites,probs = c(0.75),na.rm = T))
  #mrsos model_6_ INTERACT (without Russia) medians
  mrsos_6_W_25_2020 <- as.numeric(quantile(mrsos_2020_6_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  mrsos_6_W_50_2020 <- as.numeric(median(mrsos_2020_6_INTERACT_sites_war,na.rm = T))
  mrsos_6_W_75_2020 <- as.numeric(quantile(mrsos_2020_6_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #mrsos model_6_ Domain medians & 50% CI 2020
  mrsos_6_D_25_2020 <- as.numeric(quantile(sample(mrsos_2020_6_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  mrsos_6_D_50_2020 <- as.numeric(median(sample(mrsos_2020_6_INTERACT_domain$Value, 496),na.rm = T))
  mrsos_6_D_75_2020 <- as.numeric(quantile(sample(mrsos_2020_6_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #mrsos model_6_ Domain medians & 50% CI 2100
  mrsos_6_D_25_2100 <- as.numeric(quantile(sample(mrsos_2100_6_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  mrsos_6_D_50_2100 <- as.numeric(median(sample(mrsos_2100_6_INTERACT_domain$Value, 496),na.rm = T))
  mrsos_6_D_75_2100 <- as.numeric(quantile(sample(mrsos_2100_6_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  mrsos_6[i,] <- round(c(mrsos_6_ks_2020[1],mrsos_6_ks_2020[2],mrsos_6_ks_2020_war[1],mrsos_6_ks_2020_war[2],
                         mrsos_6_ks_2100[1],mrsos_6_ks_2100[2],mrsos_6_ks_2020_2100[1],mrsos_6_ks_2020_2100[2],
                         mrsos_6_I_25_2020,mrsos_6_I_50_2020,mrsos_6_I_75_2020,mrsos_6_W_25_2020,mrsos_6_W_50_2020,mrsos_6_W_75_2020,
                         mrsos_6_D_25_2020,mrsos_6_D_50_2020,mrsos_6_D_75_2020,mrsos_6_D_25_2100,mrsos_6_D_50_2100,mrsos_6_D_75_2100),4)
  ##############################################################################################################################
  #mrsos model_7_ KS test stats 
  mrsos_7_ks_2020 <- as.numeric(ks.test(mrsos_2020_7_INTERACT_sites, sample(mrsos_2020_7_INTERACT_domain$Value, 496))[c(1:2)])
  mrsos_7_ks_2020_war <- as.numeric(ks.test(mrsos_2020_7_INTERACT_sites_war, sample(mrsos_2020_7_INTERACT_domain$Value, 496))[c(1:2)])
  mrsos_7_ks_2100 <- as.numeric(ks.test(mrsos_2100_7_INTERACT_sites, sample(mrsos_2100_7_INTERACT_domain$Value, 496))[c(1:2)])
  mrsos_7_ks_2020_2100 <- as.numeric(ks.test(sample(mrsos_2020_7_INTERACT_domain$Value, 496), sample(mrsos_2100_7_INTERACT_domain$Value, 496))[c(1:2)])
  #mrsos model_7_ INTERACT medians
  mrsos_7_I_25_2020 <- as.numeric(quantile(mrsos_2020_7_INTERACT_sites,probs = c(0.25),na.rm = T))
  mrsos_7_I_50_2020 <- as.numeric(median(mrsos_2020_7_INTERACT_sites,na.rm = T))
  mrsos_7_I_75_2020 <- as.numeric(quantile(mrsos_2020_7_INTERACT_sites,probs = c(0.75),na.rm = T))
  #mrsos model_7_ INTERACT (without Russia) medians
  mrsos_7_W_25_2020 <- as.numeric(quantile(mrsos_2020_7_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  mrsos_7_W_50_2020 <- as.numeric(median(mrsos_2020_7_INTERACT_sites_war,na.rm = T))
  mrsos_7_W_75_2020 <- as.numeric(quantile(mrsos_2020_7_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #mrsos model_7_ Domain medians & 50% CI 2020
  mrsos_7_D_25_2020 <- as.numeric(quantile(sample(mrsos_2020_7_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  mrsos_7_D_50_2020 <- as.numeric(median(sample(mrsos_2020_7_INTERACT_domain$Value, 496),na.rm = T))
  mrsos_7_D_75_2020 <- as.numeric(quantile(sample(mrsos_2020_7_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #mrsos model_7_ Domain medians & 50% CI 2100
  mrsos_7_D_25_2100 <- as.numeric(quantile(sample(mrsos_2100_7_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  mrsos_7_D_50_2100 <- as.numeric(median(sample(mrsos_2100_7_INTERACT_domain$Value, 496),na.rm = T))
  mrsos_7_D_75_2100 <- as.numeric(quantile(sample(mrsos_2100_7_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  mrsos_7[i,] <- round(c(mrsos_7_ks_2020[1],mrsos_7_ks_2020[2],mrsos_7_ks_2020_war[1],mrsos_7_ks_2020_war[2],
                         mrsos_7_ks_2100[1],mrsos_7_ks_2100[2],mrsos_7_ks_2020_2100[1],mrsos_7_ks_2020_2100[2],
                         mrsos_7_I_25_2020,mrsos_7_I_50_2020,mrsos_7_I_75_2020,mrsos_7_W_25_2020,mrsos_7_W_50_2020,mrsos_7_W_75_2020,
                         mrsos_7_D_25_2020,mrsos_7_D_50_2020,mrsos_7_D_75_2020,mrsos_7_D_25_2100,mrsos_7_D_50_2100,mrsos_7_D_75_2100),4)
  
  ##################################################################################################################################
  #mrsos model_8_ KS test stats 
  mrsos_8_ks_2020 <- as.numeric(ks.test(mrsos_2020_8_INTERACT_sites, sample(mrsos_2020_8_INTERACT_domain$Value, 496))[c(1:2)])
  mrsos_8_ks_2020_war <- as.numeric(ks.test(mrsos_2020_8_INTERACT_sites_war, sample(mrsos_2020_8_INTERACT_domain$Value, 496))[c(1:2)])
  mrsos_8_ks_2100 <- as.numeric(ks.test(mrsos_2100_8_INTERACT_sites, sample(mrsos_2100_8_INTERACT_domain$Value, 496))[c(1:2)])
  mrsos_8_ks_2020_2100 <- as.numeric(ks.test(sample(mrsos_2020_8_INTERACT_domain$Value, 496), sample(mrsos_2100_8_INTERACT_domain$Value, 496))[c(1:2)])
  #mrsos model_8_ INTERACT medians
  mrsos_8_I_25_2020 <- as.numeric(quantile(mrsos_2020_8_INTERACT_sites,probs = c(0.25),na.rm = T))
  mrsos_8_I_50_2020 <- as.numeric(median(mrsos_2020_8_INTERACT_sites,na.rm = T))
  mrsos_8_I_75_2020 <- as.numeric(quantile(mrsos_2020_8_INTERACT_sites,probs = c(0.75),na.rm = T))
  #mrsos model_8_ INTERACT (without Russia) medians
  mrsos_8_W_25_2020 <- as.numeric(quantile(mrsos_2020_8_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  mrsos_8_W_50_2020 <- as.numeric(median(mrsos_2020_8_INTERACT_sites_war,na.rm = T))
  mrsos_8_W_75_2020 <- as.numeric(quantile(mrsos_2020_8_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #mrsos model_8_ Domain medians & 50% CI 2020
  mrsos_8_D_25_2020 <- as.numeric(quantile(sample(mrsos_2020_8_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  mrsos_8_D_50_2020 <- as.numeric(median(sample(mrsos_2020_8_INTERACT_domain$Value, 496),na.rm = T))
  mrsos_8_D_75_2020 <- as.numeric(quantile(sample(mrsos_2020_8_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #mrsos model_8_ Domain medians & 50% CI 2100
  mrsos_8_D_25_2100 <- as.numeric(quantile(sample(mrsos_2100_8_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  mrsos_8_D_50_2100 <- as.numeric(median(sample(mrsos_2100_8_INTERACT_domain$Value, 496),na.rm = T))
  mrsos_8_D_75_2100 <- as.numeric(quantile(sample(mrsos_2100_8_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  mrsos_8[i,] <- round(c(mrsos_8_ks_2020[1],mrsos_8_ks_2020[2],mrsos_8_ks_2020_war[1],mrsos_8_ks_2020_war[2],
                         mrsos_8_ks_2100[1],mrsos_8_ks_2100[2],mrsos_8_ks_2020_2100[1],mrsos_8_ks_2020_2100[2],
                         mrsos_8_I_25_2020,mrsos_8_I_50_2020,mrsos_8_I_75_2020,mrsos_8_W_25_2020,mrsos_8_W_50_2020,mrsos_8_W_75_2020,
                         mrsos_8_D_25_2020,mrsos_8_D_50_2020,mrsos_8_D_75_2020,mrsos_8_D_25_2100,mrsos_8_D_50_2100,mrsos_8_D_75_2100),4)
  ##################################################################################################################################
  #cveg model_1_ KS test stats 
  cveg_1_ks_2020 <- as.numeric(ks.test(cveg_2020_1_INTERACT_sites, sample(cveg_2020_1_INTERACT_domain$Value, 496))[c(1:2)])
  cveg_1_ks_2020_war <- as.numeric(ks.test(cveg_2020_1_INTERACT_sites_war, sample(cveg_2020_1_INTERACT_domain$Value, 496))[c(1:2)])
  cveg_1_ks_2100 <- as.numeric(ks.test(cveg_2100_1_INTERACT_sites, sample(cveg_2100_1_INTERACT_domain$Value, 496))[c(1:2)])
  cveg_1_ks_2020_2100 <- as.numeric(ks.test(sample(cveg_2020_1_INTERACT_domain$Value, 496), sample(cveg_2100_1_INTERACT_domain$Value, 496))[c(1:2)])
  #cveg model_1_ INTERACT medians
  cveg_1_I_25_2020 <- as.numeric(quantile(cveg_2020_1_INTERACT_sites,probs = c(0.25),na.rm = T))
  cveg_1_I_50_2020 <- as.numeric(median(cveg_2020_1_INTERACT_sites,na.rm = T))
  cveg_1_I_75_2020 <- as.numeric(quantile(cveg_2020_1_INTERACT_sites,probs = c(0.75),na.rm = T))
  #cveg model_1_ INTERACT (without Russia) medians
  cveg_1_W_25_2020 <- as.numeric(quantile(cveg_2020_1_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  cveg_1_W_50_2020 <- as.numeric(median(cveg_2020_1_INTERACT_sites_war,na.rm = T))
  cveg_1_W_75_2020 <- as.numeric(quantile(cveg_2020_1_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #cveg model_1_ Domain medians & 50% CI 2020
  cveg_1_D_25_2020 <- as.numeric(quantile(sample(cveg_2020_1_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  cveg_1_D_50_2020 <- as.numeric(median(sample(cveg_2020_1_INTERACT_domain$Value, 496),na.rm = T))
  cveg_1_D_75_2020 <- as.numeric(quantile(sample(cveg_2020_1_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #cveg model_1_ Domain medians & 50% CI 2100
  cveg_1_D_25_2100 <- as.numeric(quantile(sample(cveg_2100_1_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  cveg_1_D_50_2100 <- as.numeric(median(sample(cveg_2100_1_INTERACT_domain$Value, 496),na.rm = T))
  cveg_1_D_75_2100 <- as.numeric(quantile(sample(cveg_2100_1_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  cveg_1[i,] <- round(c(cveg_1_ks_2020[1],cveg_1_ks_2020[2],cveg_1_ks_2020_war[1],cveg_1_ks_2020_war[2],
                        cveg_1_ks_2100[1],cveg_1_ks_2100[2],cveg_1_ks_2020_2100[1],cveg_1_ks_2020_2100[2],
                        cveg_1_I_25_2020,cveg_1_I_50_2020,cveg_1_I_75_2020,cveg_1_W_25_2020,cveg_1_W_50_2020,cveg_1_W_75_2020,
                        cveg_1_D_25_2020,cveg_1_D_50_2020,cveg_1_D_75_2020,cveg_1_D_25_2100,cveg_1_D_50_2100,cveg_1_D_75_2100),4)
  ##################################################################################################################################
  #cveg model_2_ KS test stats 
  cveg_2_ks_2020 <- as.numeric(ks.test(cveg_2020_2_INTERACT_sites, sample(cveg_2020_2_INTERACT_domain$Value, 496))[c(1:2)])
  cveg_2_ks_2020_war <- as.numeric(ks.test(cveg_2020_2_INTERACT_sites_war, sample(cveg_2020_2_INTERACT_domain$Value, 496))[c(1:2)])
  cveg_2_ks_2100 <- as.numeric(ks.test(cveg_2100_2_INTERACT_sites, sample(cveg_2100_2_INTERACT_domain$Value, 496))[c(1:2)])
  cveg_2_ks_2020_2100 <- as.numeric(ks.test(sample(cveg_2020_2_INTERACT_domain$Value, 496), sample(cveg_2100_2_INTERACT_domain$Value, 496))[c(1:2)])
  #cveg model_2_ INTERACT medians
  cveg_2_I_25_2020 <- as.numeric(quantile(cveg_2020_2_INTERACT_sites,probs = c(0.25),na.rm = T))
  cveg_2_I_50_2020 <- as.numeric(median(cveg_2020_2_INTERACT_sites,na.rm = T))
  cveg_2_I_75_2020 <- as.numeric(quantile(cveg_2020_2_INTERACT_sites,probs = c(0.75),na.rm = T))
  #cveg model_2_ INTERACT (without Russia) medians
  cveg_2_W_25_2020 <- as.numeric(quantile(cveg_2020_2_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  cveg_2_W_50_2020 <- as.numeric(median(cveg_2020_2_INTERACT_sites_war,na.rm = T))
  cveg_2_W_75_2020 <- as.numeric(quantile(cveg_2020_2_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #cveg model_2_ Domain medians & 50% CI 2020
  cveg_2_D_25_2020 <- as.numeric(quantile(sample(cveg_2020_2_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  cveg_2_D_50_2020 <- as.numeric(median(sample(cveg_2020_2_INTERACT_domain$Value, 496),na.rm = T))
  cveg_2_D_75_2020 <- as.numeric(quantile(sample(cveg_2020_2_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #cveg model_2_ Domain medians & 50% CI 2100
  cveg_2_D_25_2100 <- as.numeric(quantile(sample(cveg_2100_2_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  cveg_2_D_50_2100 <- as.numeric(median(sample(cveg_2100_2_INTERACT_domain$Value, 496),na.rm = T))
  cveg_2_D_75_2100 <- as.numeric(quantile(sample(cveg_2100_2_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  cveg_2[i,] <- round(c(cveg_2_ks_2020[1],cveg_2_ks_2020[2],cveg_2_ks_2020_war[1],cveg_2_ks_2020_war[2],
                        cveg_2_ks_2100[1],cveg_2_ks_2100[2],cveg_2_ks_2020_2100[1],cveg_2_ks_2020_2100[2],
                        cveg_2_I_25_2020,cveg_2_I_50_2020,cveg_2_I_75_2020,cveg_2_W_25_2020,cveg_2_W_50_2020,cveg_2_W_75_2020,
                        cveg_2_D_25_2020,cveg_2_D_50_2020,cveg_2_D_75_2020,cveg_2_D_25_2100,cveg_2_D_50_2100,cveg_2_D_75_2100),4)
  ##############################################################################################################################
  #cveg model_3_ KS test stats 
  cveg_3_ks_2020 <- as.numeric(ks.test(cveg_2020_3_INTERACT_sites, sample(cveg_2020_3_INTERACT_domain$Value, 496))[c(1:2)])
  cveg_3_ks_2020_war <- as.numeric(ks.test(cveg_2020_3_INTERACT_sites_war, sample(cveg_2020_3_INTERACT_domain$Value, 496))[c(1:2)])
  cveg_3_ks_2100 <- as.numeric(ks.test(cveg_2100_3_INTERACT_sites, sample(cveg_2100_3_INTERACT_domain$Value, 496))[c(1:2)])
  cveg_3_ks_2020_2100 <- as.numeric(ks.test(sample(cveg_2020_3_INTERACT_domain$Value, 496), sample(cveg_2100_3_INTERACT_domain$Value, 496))[c(1:2)])
  #cveg model_3_ INTERACT medians
  cveg_3_I_25_2020 <- as.numeric(quantile(cveg_2020_3_INTERACT_sites,probs = c(0.25),na.rm = T))
  cveg_3_I_50_2020 <- as.numeric(median(cveg_2020_3_INTERACT_sites,na.rm = T))
  cveg_3_I_75_2020 <- as.numeric(quantile(cveg_2020_3_INTERACT_sites,probs = c(0.75),na.rm = T))
  #cveg model_3_ INTERACT (without Russia) medians
  cveg_3_W_25_2020 <- as.numeric(quantile(cveg_2020_3_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  cveg_3_W_50_2020 <- as.numeric(median(cveg_2020_3_INTERACT_sites_war,na.rm = T))
  cveg_3_W_75_2020 <- as.numeric(quantile(cveg_2020_3_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #cveg model_3_ Domain medians & 50% CI 2020
  cveg_3_D_25_2020 <- as.numeric(quantile(sample(cveg_2020_3_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  cveg_3_D_50_2020 <- as.numeric(median(sample(cveg_2020_3_INTERACT_domain$Value, 496),na.rm = T))
  cveg_3_D_75_2020 <- as.numeric(quantile(sample(cveg_2020_3_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #cveg model_3_ Domain medians & 50% CI 2100
  cveg_3_D_25_2100 <- as.numeric(quantile(sample(cveg_2100_3_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  cveg_3_D_50_2100 <- as.numeric(median(sample(cveg_2100_3_INTERACT_domain$Value, 496),na.rm = T))
  cveg_3_D_75_2100 <- as.numeric(quantile(sample(cveg_2100_3_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  cveg_3[i,] <- round(c(cveg_3_ks_2020[1],cveg_3_ks_2020[2],cveg_3_ks_2020_war[1],cveg_3_ks_2020_war[2],
                        cveg_3_ks_2100[1],cveg_3_ks_2100[2],cveg_3_ks_2020_2100[1],cveg_3_ks_2020_2100[2],
                        cveg_3_I_25_2020,cveg_3_I_50_2020,cveg_3_I_75_2020,cveg_3_W_25_2020,cveg_3_W_50_2020,cveg_3_W_75_2020,
                        cveg_3_D_25_2020,cveg_3_D_50_2020,cveg_3_D_75_2020,cveg_3_D_25_2100,cveg_3_D_50_2100,cveg_3_D_75_2100),4)
  ##############################################################################################################################
  #cveg model_4_ KS test stats 
  cveg_4_ks_2020 <- as.numeric(ks.test(cveg_2020_4_INTERACT_sites, sample(cveg_2020_4_INTERACT_domain$Value, 496))[c(1:2)])
  cveg_4_ks_2020_war <- as.numeric(ks.test(cveg_2020_4_INTERACT_sites_war, sample(cveg_2020_4_INTERACT_domain$Value, 496))[c(1:2)])
  cveg_4_ks_2100 <- as.numeric(ks.test(cveg_2100_4_INTERACT_sites, sample(cveg_2100_4_INTERACT_domain$Value, 496))[c(1:2)])
  cveg_4_ks_2020_2100 <- as.numeric(ks.test(sample(cveg_2020_4_INTERACT_domain$Value, 496), sample(cveg_2100_4_INTERACT_domain$Value, 496))[c(1:2)])
  #cveg model_4_ INTERACT medians
  cveg_4_I_25_2020 <- as.numeric(quantile(cveg_2020_4_INTERACT_sites,probs = c(0.25),na.rm = T))
  cveg_4_I_50_2020 <- as.numeric(median(cveg_2020_4_INTERACT_sites,na.rm = T))
  cveg_4_I_75_2020 <- as.numeric(quantile(cveg_2020_4_INTERACT_sites,probs = c(0.75),na.rm = T))
  #cveg model_4_ INTERACT (without Russia) medians
  cveg_4_W_25_2020 <- as.numeric(quantile(cveg_2020_4_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  cveg_4_W_50_2020 <- as.numeric(median(cveg_2020_4_INTERACT_sites_war,na.rm = T))
  cveg_4_W_75_2020 <- as.numeric(quantile(cveg_2020_4_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #cveg model_4_ Domain medians & 50% CI 2020
  cveg_4_D_25_2020 <- as.numeric(quantile(sample(cveg_2020_4_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  cveg_4_D_50_2020 <- as.numeric(median(sample(cveg_2020_4_INTERACT_domain$Value, 496),na.rm = T))
  cveg_4_D_75_2020 <- as.numeric(quantile(sample(cveg_2020_4_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #cveg model_4_ Domain medians & 50% CI 2100
  cveg_4_D_25_2100 <- as.numeric(quantile(sample(cveg_2100_4_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  cveg_4_D_50_2100 <- as.numeric(median(sample(cveg_2100_4_INTERACT_domain$Value, 496),na.rm = T))
  cveg_4_D_75_2100 <- as.numeric(quantile(sample(cveg_2100_4_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  cveg_4[i,] <- round(c(cveg_4_ks_2020[1],cveg_4_ks_2020[2],cveg_4_ks_2020_war[1],cveg_4_ks_2020_war[2],
                        cveg_4_ks_2100[1],cveg_4_ks_2100[2],cveg_4_ks_2020_2100[1],cveg_4_ks_2020_2100[2],
                        cveg_4_I_25_2020,cveg_4_I_50_2020,cveg_4_I_75_2020,cveg_4_W_25_2020,cveg_4_W_50_2020,cveg_4_W_75_2020,
                        cveg_4_D_25_2020,cveg_4_D_50_2020,cveg_4_D_75_2020,cveg_4_D_25_2100,cveg_4_D_50_2100,cveg_4_D_75_2100),4)
  ##############################################################################################################################
  #cveg model_5_ KS test stats 
  cveg_5_ks_2020 <- as.numeric(ks.test(cveg_2020_5_INTERACT_sites, sample(cveg_2020_5_INTERACT_domain$Value, 496))[c(1:2)])
  cveg_5_ks_2020_war <- as.numeric(ks.test(cveg_2020_5_INTERACT_sites_war, sample(cveg_2020_5_INTERACT_domain$Value, 496))[c(1:2)])
  cveg_5_ks_2100 <- as.numeric(ks.test(cveg_2100_5_INTERACT_sites, sample(cveg_2100_5_INTERACT_domain$Value, 496))[c(1:2)])
  cveg_5_ks_2020_2100 <- as.numeric(ks.test(sample(cveg_2020_5_INTERACT_domain$Value, 496), sample(cveg_2100_5_INTERACT_domain$Value, 496))[c(1:2)])
  #cveg model_5_ INTERACT medians
  cveg_5_I_25_2020 <- as.numeric(quantile(cveg_2020_5_INTERACT_sites,probs = c(0.25),na.rm = T))
  cveg_5_I_50_2020 <- as.numeric(median(cveg_2020_5_INTERACT_sites,na.rm = T))
  cveg_5_I_75_2020 <- as.numeric(quantile(cveg_2020_5_INTERACT_sites,probs = c(0.75),na.rm = T))
  #cveg model_5_ INTERACT (without Russia) medians
  cveg_5_W_25_2020 <- as.numeric(quantile(cveg_2020_5_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  cveg_5_W_50_2020 <- as.numeric(median(cveg_2020_5_INTERACT_sites_war,na.rm = T))
  cveg_5_W_75_2020 <- as.numeric(quantile(cveg_2020_5_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #cveg model_5_ Domain medians & 50% CI 2020
  cveg_5_D_25_2020 <- as.numeric(quantile(sample(cveg_2020_5_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  cveg_5_D_50_2020 <- as.numeric(median(sample(cveg_2020_5_INTERACT_domain$Value, 496),na.rm = T))
  cveg_5_D_75_2020 <- as.numeric(quantile(sample(cveg_2020_5_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #cveg model_5_ Domain medians & 50% CI 2100
  cveg_5_D_25_2100 <- as.numeric(quantile(sample(cveg_2100_5_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  cveg_5_D_50_2100 <- as.numeric(median(sample(cveg_2100_5_INTERACT_domain$Value, 496),na.rm = T))
  cveg_5_D_75_2100 <- as.numeric(quantile(sample(cveg_2100_5_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  cveg_5[i,] <- round(c(cveg_5_ks_2020[1],cveg_5_ks_2020[2],cveg_5_ks_2020_war[1],cveg_5_ks_2020_war[2],
                        cveg_5_ks_2100[1],cveg_5_ks_2100[2],cveg_5_ks_2020_2100[1],cveg_5_ks_2020_2100[2],
                        cveg_5_I_25_2020,cveg_5_I_50_2020,cveg_5_I_75_2020,cveg_5_W_25_2020,cveg_5_W_50_2020,cveg_5_W_75_2020,
                        cveg_5_D_25_2020,cveg_5_D_50_2020,cveg_5_D_75_2020,cveg_5_D_25_2100,cveg_5_D_50_2100,cveg_5_D_75_2100),4)
  ##############################################################################################################################
  #cveg model_6_ KS test stats 
  cveg_6_ks_2020 <- as.numeric(ks.test(cveg_2020_6_INTERACT_sites, sample(cveg_2020_6_INTERACT_domain$Value, 496))[c(1:2)])
  cveg_6_ks_2020_war <- as.numeric(ks.test(cveg_2020_6_INTERACT_sites_war, sample(cveg_2020_6_INTERACT_domain$Value, 496))[c(1:2)])
  cveg_6_ks_2100 <- as.numeric(ks.test(cveg_2100_6_INTERACT_sites, sample(cveg_2100_6_INTERACT_domain$Value, 496))[c(1:2)])
  cveg_6_ks_2020_2100 <- as.numeric(ks.test(sample(cveg_2020_6_INTERACT_domain$Value, 496), sample(cveg_2100_6_INTERACT_domain$Value, 496))[c(1:2)])
  #cveg model_6_ INTERACT medians
  cveg_6_I_25_2020 <- as.numeric(quantile(cveg_2020_6_INTERACT_sites,probs = c(0.25),na.rm = T))
  cveg_6_I_50_2020 <- as.numeric(median(cveg_2020_6_INTERACT_sites,na.rm = T))
  cveg_6_I_75_2020 <- as.numeric(quantile(cveg_2020_6_INTERACT_sites,probs = c(0.75),na.rm = T))
  #cveg model_6_ INTERACT (without Russia) medians
  cveg_6_W_25_2020 <- as.numeric(quantile(cveg_2020_6_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  cveg_6_W_50_2020 <- as.numeric(median(cveg_2020_6_INTERACT_sites_war,na.rm = T))
  cveg_6_W_75_2020 <- as.numeric(quantile(cveg_2020_6_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #cveg model_6_ Domain medians & 50% CI 2020
  cveg_6_D_25_2020 <- as.numeric(quantile(sample(cveg_2020_6_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  cveg_6_D_50_2020 <- as.numeric(median(sample(cveg_2020_6_INTERACT_domain$Value, 496),na.rm = T))
  cveg_6_D_75_2020 <- as.numeric(quantile(sample(cveg_2020_6_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #cveg model_6_ Domain medians & 50% CI 2100
  cveg_6_D_25_2100 <- as.numeric(quantile(sample(cveg_2100_6_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  cveg_6_D_50_2100 <- as.numeric(median(sample(cveg_2100_6_INTERACT_domain$Value, 496),na.rm = T))
  cveg_6_D_75_2100 <- as.numeric(quantile(sample(cveg_2100_6_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  cveg_6[i,] <- round(c(cveg_6_ks_2020[1],cveg_6_ks_2020[2],cveg_6_ks_2020_war[1],cveg_6_ks_2020_war[2],
                        cveg_6_ks_2100[1],cveg_6_ks_2100[2],cveg_6_ks_2020_2100[1],cveg_6_ks_2020_2100[2],
                        cveg_6_I_25_2020,cveg_6_I_50_2020,cveg_6_I_75_2020,cveg_6_W_25_2020,cveg_6_W_50_2020,cveg_6_W_75_2020,
                        cveg_6_D_25_2020,cveg_6_D_50_2020,cveg_6_D_75_2020,cveg_6_D_25_2100,cveg_6_D_50_2100,cveg_6_D_75_2100),4)
  ##############################################################################################################################
  #cveg model_7_ KS test stats 
  cveg_7_ks_2020 <- as.numeric(ks.test(cveg_2020_7_INTERACT_sites, sample(cveg_2020_7_INTERACT_domain$Value, 496))[c(1:2)])
  cveg_7_ks_2020_war <- as.numeric(ks.test(cveg_2020_7_INTERACT_sites_war, sample(cveg_2020_7_INTERACT_domain$Value, 496))[c(1:2)])
  cveg_7_ks_2100 <- as.numeric(ks.test(cveg_2100_7_INTERACT_sites, sample(cveg_2100_7_INTERACT_domain$Value, 496))[c(1:2)])
  cveg_7_ks_2020_2100 <- as.numeric(ks.test(sample(cveg_2020_7_INTERACT_domain$Value, 496), sample(cveg_2100_7_INTERACT_domain$Value, 496))[c(1:2)])
  #cveg model_7_ INTERACT medians
  cveg_7_I_25_2020 <- as.numeric(quantile(cveg_2020_7_INTERACT_sites,probs = c(0.25),na.rm = T))
  cveg_7_I_50_2020 <- as.numeric(median(cveg_2020_7_INTERACT_sites,na.rm = T))
  cveg_7_I_75_2020 <- as.numeric(quantile(cveg_2020_7_INTERACT_sites,probs = c(0.75),na.rm = T))
  #cveg model_7_ INTERACT (without Russia) medians
  cveg_7_W_25_2020 <- as.numeric(quantile(cveg_2020_7_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  cveg_7_W_50_2020 <- as.numeric(median(cveg_2020_7_INTERACT_sites_war,na.rm = T))
  cveg_7_W_75_2020 <- as.numeric(quantile(cveg_2020_7_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #cveg model_7_ Domain medians & 50% CI 2020
  cveg_7_D_25_2020 <- as.numeric(quantile(sample(cveg_2020_7_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  cveg_7_D_50_2020 <- as.numeric(median(sample(cveg_2020_7_INTERACT_domain$Value, 496),na.rm = T))
  cveg_7_D_75_2020 <- as.numeric(quantile(sample(cveg_2020_7_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #cveg model_7_ Domain medians & 50% CI 2100
  cveg_7_D_25_2100 <- as.numeric(quantile(sample(cveg_2100_7_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  cveg_7_D_50_2100 <- as.numeric(median(sample(cveg_2100_7_INTERACT_domain$Value, 496),na.rm = T))
  cveg_7_D_75_2100 <- as.numeric(quantile(sample(cveg_2100_7_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  cveg_7[i,] <- round(c(cveg_7_ks_2020[1],cveg_7_ks_2020[2],cveg_7_ks_2020_war[1],cveg_7_ks_2020_war[2],
                        cveg_7_ks_2100[1],cveg_7_ks_2100[2],cveg_7_ks_2020_2100[1],cveg_7_ks_2020_2100[2],
                        cveg_7_I_25_2020,cveg_7_I_50_2020,cveg_7_I_75_2020,cveg_7_W_25_2020,cveg_7_W_50_2020,cveg_7_W_75_2020,
                        cveg_7_D_25_2020,cveg_7_D_50_2020,cveg_7_D_75_2020,cveg_7_D_25_2100,cveg_7_D_50_2100,cveg_7_D_75_2100),4)
  
  ##################################################################################################################################
  #cveg model_8_ KS test stats 
  cveg_8_ks_2020 <- as.numeric(ks.test(cveg_2020_8_INTERACT_sites, sample(cveg_2020_8_INTERACT_domain$Value, 496))[c(1:2)])
  cveg_8_ks_2020_war <- as.numeric(ks.test(cveg_2020_8_INTERACT_sites_war, sample(cveg_2020_8_INTERACT_domain$Value, 496))[c(1:2)])
  cveg_8_ks_2100 <- as.numeric(ks.test(cveg_2100_8_INTERACT_sites, sample(cveg_2100_8_INTERACT_domain$Value, 496))[c(1:2)])
  cveg_8_ks_2020_2100 <- as.numeric(ks.test(sample(cveg_2020_8_INTERACT_domain$Value, 496), sample(cveg_2100_8_INTERACT_domain$Value, 496))[c(1:2)])
  #cveg model_8_ INTERACT medians
  cveg_8_I_25_2020 <- as.numeric(quantile(cveg_2020_8_INTERACT_sites,probs = c(0.25),na.rm = T))
  cveg_8_I_50_2020 <- as.numeric(median(cveg_2020_8_INTERACT_sites,na.rm = T))
  cveg_8_I_75_2020 <- as.numeric(quantile(cveg_2020_8_INTERACT_sites,probs = c(0.75),na.rm = T))
  #cveg model_8_ INTERACT (without Russia) medians
  cveg_8_W_25_2020 <- as.numeric(quantile(cveg_2020_8_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  cveg_8_W_50_2020 <- as.numeric(median(cveg_2020_8_INTERACT_sites_war,na.rm = T))
  cveg_8_W_75_2020 <- as.numeric(quantile(cveg_2020_8_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #cveg model_8_ Domain medians & 50% CI 2020
  cveg_8_D_25_2020 <- as.numeric(quantile(sample(cveg_2020_8_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  cveg_8_D_50_2020 <- as.numeric(median(sample(cveg_2020_8_INTERACT_domain$Value, 496),na.rm = T))
  cveg_8_D_75_2020 <- as.numeric(quantile(sample(cveg_2020_8_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #cveg model_8_ Domain medians & 50% CI 2100
  cveg_8_D_25_2100 <- as.numeric(quantile(sample(cveg_2100_8_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  cveg_8_D_50_2100 <- as.numeric(median(sample(cveg_2100_8_INTERACT_domain$Value, 496),na.rm = T))
  cveg_8_D_75_2100 <- as.numeric(quantile(sample(cveg_2100_8_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  cveg_8[i,] <- round(c(cveg_8_ks_2020[1],cveg_8_ks_2020[2],cveg_8_ks_2020_war[1],cveg_8_ks_2020_war[2],
                        cveg_8_ks_2100[1],cveg_8_ks_2100[2],cveg_8_ks_2020_2100[1],cveg_8_ks_2020_2100[2],
                        cveg_8_I_25_2020,cveg_8_I_50_2020,cveg_8_I_75_2020,cveg_8_W_25_2020,cveg_8_W_50_2020,cveg_8_W_75_2020,
                        cveg_8_D_25_2020,cveg_8_D_50_2020,cveg_8_D_75_2020,cveg_8_D_25_2100,cveg_8_D_50_2100,cveg_8_D_75_2100),4)
  ##################################################################################################################################
  #csoil model_1_ KS test stats 
  csoil_1_ks_2020 <- as.numeric(ks.test(csoil_2020_1_INTERACT_sites, sample(csoil_2020_1_INTERACT_domain$Value, 496))[c(1:2)])
  csoil_1_ks_2020_war <- as.numeric(ks.test(csoil_2020_1_INTERACT_sites_war, sample(csoil_2020_1_INTERACT_domain$Value, 496))[c(1:2)])
  csoil_1_ks_2100 <- as.numeric(ks.test(csoil_2100_1_INTERACT_sites, sample(csoil_2100_1_INTERACT_domain$Value, 496))[c(1:2)])
  csoil_1_ks_2020_2100 <- as.numeric(ks.test(sample(csoil_2020_1_INTERACT_domain$Value, 496), sample(csoil_2100_1_INTERACT_domain$Value, 496))[c(1:2)])
  #csoil model_1_ INTERACT medians
  csoil_1_I_25_2020 <- as.numeric(quantile(csoil_2020_1_INTERACT_sites,probs = c(0.25),na.rm = T))
  csoil_1_I_50_2020 <- as.numeric(median(csoil_2020_1_INTERACT_sites,na.rm = T))
  csoil_1_I_75_2020 <- as.numeric(quantile(csoil_2020_1_INTERACT_sites,probs = c(0.75),na.rm = T))
  #csoil model_1_ INTERACT (without Russia) medians
  csoil_1_W_25_2020 <- as.numeric(quantile(csoil_2020_1_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  csoil_1_W_50_2020 <- as.numeric(median(csoil_2020_1_INTERACT_sites_war,na.rm = T))
  csoil_1_W_75_2020 <- as.numeric(quantile(csoil_2020_1_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #csoil model_1_ Domain medians & 50% CI 2020
  csoil_1_D_25_2020 <- as.numeric(quantile(sample(csoil_2020_1_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  csoil_1_D_50_2020 <- as.numeric(median(sample(csoil_2020_1_INTERACT_domain$Value, 496),na.rm = T))
  csoil_1_D_75_2020 <- as.numeric(quantile(sample(csoil_2020_1_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #csoil model_1_ Domain medians & 50% CI 2100
  csoil_1_D_25_2100 <- as.numeric(quantile(sample(csoil_2100_1_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  csoil_1_D_50_2100 <- as.numeric(median(sample(csoil_2100_1_INTERACT_domain$Value, 496),na.rm = T))
  csoil_1_D_75_2100 <- as.numeric(quantile(sample(csoil_2100_1_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  csoil_1[i,] <- round(c(csoil_1_ks_2020[1],csoil_1_ks_2020[2],csoil_1_ks_2020_war[1],csoil_1_ks_2020_war[2],
                         csoil_1_ks_2100[1],csoil_1_ks_2100[2],csoil_1_ks_2020_2100[1],csoil_1_ks_2020_2100[2],
                         csoil_1_I_25_2020,csoil_1_I_50_2020,csoil_1_I_75_2020,csoil_1_W_25_2020,csoil_1_W_50_2020,csoil_1_W_75_2020,
                         csoil_1_D_25_2020,csoil_1_D_50_2020,csoil_1_D_75_2020,csoil_1_D_25_2100,csoil_1_D_50_2100,csoil_1_D_75_2100),4)
  ##################################################################################################################################
  #csoil model_2_ KS test stats 
  csoil_2_ks_2020 <- as.numeric(ks.test(csoil_2020_2_INTERACT_sites, sample(csoil_2020_2_INTERACT_domain$Value, 496))[c(1:2)])
  csoil_2_ks_2020_war <- as.numeric(ks.test(csoil_2020_2_INTERACT_sites_war, sample(csoil_2020_2_INTERACT_domain$Value, 496))[c(1:2)])
  csoil_2_ks_2100 <- as.numeric(ks.test(csoil_2100_2_INTERACT_sites, sample(csoil_2100_2_INTERACT_domain$Value, 496))[c(1:2)])
  csoil_2_ks_2020_2100 <- as.numeric(ks.test(sample(csoil_2020_2_INTERACT_domain$Value, 496), sample(csoil_2100_2_INTERACT_domain$Value, 496))[c(1:2)])
  #csoil model_2_ INTERACT medians
  csoil_2_I_25_2020 <- as.numeric(quantile(csoil_2020_2_INTERACT_sites,probs = c(0.25),na.rm = T))
  csoil_2_I_50_2020 <- as.numeric(median(csoil_2020_2_INTERACT_sites,na.rm = T))
  csoil_2_I_75_2020 <- as.numeric(quantile(csoil_2020_2_INTERACT_sites,probs = c(0.75),na.rm = T))
  #csoil model_2_ INTERACT (without Russia) medians
  csoil_2_W_25_2020 <- as.numeric(quantile(csoil_2020_2_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  csoil_2_W_50_2020 <- as.numeric(median(csoil_2020_2_INTERACT_sites_war,na.rm = T))
  csoil_2_W_75_2020 <- as.numeric(quantile(csoil_2020_2_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #csoil model_2_ Domain medians & 50% CI 2020
  csoil_2_D_25_2020 <- as.numeric(quantile(sample(csoil_2020_2_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  csoil_2_D_50_2020 <- as.numeric(median(sample(csoil_2020_2_INTERACT_domain$Value, 496),na.rm = T))
  csoil_2_D_75_2020 <- as.numeric(quantile(sample(csoil_2020_2_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #csoil model_2_ Domain medians & 50% CI 2100
  csoil_2_D_25_2100 <- as.numeric(quantile(sample(csoil_2100_2_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  csoil_2_D_50_2100 <- as.numeric(median(sample(csoil_2100_2_INTERACT_domain$Value, 496),na.rm = T))
  csoil_2_D_75_2100 <- as.numeric(quantile(sample(csoil_2100_2_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  csoil_2[i,] <- round(c(csoil_2_ks_2020[1],csoil_2_ks_2020[2],csoil_2_ks_2020_war[1],csoil_2_ks_2020_war[2],
                         csoil_2_ks_2100[1],csoil_2_ks_2100[2],csoil_2_ks_2020_2100[1],csoil_2_ks_2020_2100[2],
                         csoil_2_I_25_2020,csoil_2_I_50_2020,csoil_2_I_75_2020,csoil_2_W_25_2020,csoil_2_W_50_2020,csoil_2_W_75_2020,
                         csoil_2_D_25_2020,csoil_2_D_50_2020,csoil_2_D_75_2020,csoil_2_D_25_2100,csoil_2_D_50_2100,csoil_2_D_75_2100),4)
  ##############################################################################################################################
  #csoil model_3_ KS test stats 
  csoil_3_ks_2020 <- as.numeric(ks.test(csoil_2020_3_INTERACT_sites, sample(csoil_2020_3_INTERACT_domain$Value, 496))[c(1:2)])
  csoil_3_ks_2020_war <- as.numeric(ks.test(csoil_2020_3_INTERACT_sites_war, sample(csoil_2020_3_INTERACT_domain$Value, 496))[c(1:2)])
  csoil_3_ks_2100 <- as.numeric(ks.test(csoil_2100_3_INTERACT_sites, sample(csoil_2100_3_INTERACT_domain$Value, 496))[c(1:2)])
  csoil_3_ks_2020_2100 <- as.numeric(ks.test(sample(csoil_2020_3_INTERACT_domain$Value, 496), sample(csoil_2100_3_INTERACT_domain$Value, 496))[c(1:2)])
  #csoil model_3_ INTERACT medians
  csoil_3_I_25_2020 <- as.numeric(quantile(csoil_2020_3_INTERACT_sites,probs = c(0.25),na.rm = T))
  csoil_3_I_50_2020 <- as.numeric(median(csoil_2020_3_INTERACT_sites,na.rm = T))
  csoil_3_I_75_2020 <- as.numeric(quantile(csoil_2020_3_INTERACT_sites,probs = c(0.75),na.rm = T))
  #csoil model_3_ INTERACT (without Russia) medians
  csoil_3_W_25_2020 <- as.numeric(quantile(csoil_2020_3_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  csoil_3_W_50_2020 <- as.numeric(median(csoil_2020_3_INTERACT_sites_war,na.rm = T))
  csoil_3_W_75_2020 <- as.numeric(quantile(csoil_2020_3_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #csoil model_3_ Domain medians & 50% CI 2020
  csoil_3_D_25_2020 <- as.numeric(quantile(sample(csoil_2020_3_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  csoil_3_D_50_2020 <- as.numeric(median(sample(csoil_2020_3_INTERACT_domain$Value, 496),na.rm = T))
  csoil_3_D_75_2020 <- as.numeric(quantile(sample(csoil_2020_3_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #csoil model_3_ Domain medians & 50% CI 2100
  csoil_3_D_25_2100 <- as.numeric(quantile(sample(csoil_2100_3_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  csoil_3_D_50_2100 <- as.numeric(median(sample(csoil_2100_3_INTERACT_domain$Value, 496),na.rm = T))
  csoil_3_D_75_2100 <- as.numeric(quantile(sample(csoil_2100_3_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  csoil_3[i,] <- round(c(csoil_3_ks_2020[1],csoil_3_ks_2020[2],csoil_3_ks_2020_war[1],csoil_3_ks_2020_war[2],
                         csoil_3_ks_2100[1],csoil_3_ks_2100[2],csoil_3_ks_2020_2100[1],csoil_3_ks_2020_2100[2],
                         csoil_3_I_25_2020,csoil_3_I_50_2020,csoil_3_I_75_2020,csoil_3_W_25_2020,csoil_3_W_50_2020,csoil_3_W_75_2020,
                         csoil_3_D_25_2020,csoil_3_D_50_2020,csoil_3_D_75_2020,csoil_3_D_25_2100,csoil_3_D_50_2100,csoil_3_D_75_2100),4)
  ##############################################################################################################################
  #csoil model_4_ KS test stats 
  csoil_4_ks_2020 <- as.numeric(ks.test(csoil_2020_4_INTERACT_sites, sample(csoil_2020_4_INTERACT_domain$Value, 496))[c(1:2)])
  csoil_4_ks_2020_war <- as.numeric(ks.test(csoil_2020_4_INTERACT_sites_war, sample(csoil_2020_4_INTERACT_domain$Value, 496))[c(1:2)])
  csoil_4_ks_2100 <- as.numeric(ks.test(csoil_2100_4_INTERACT_sites, sample(csoil_2100_4_INTERACT_domain$Value, 496))[c(1:2)])
  csoil_4_ks_2020_2100 <- as.numeric(ks.test(sample(csoil_2020_4_INTERACT_domain$Value, 496), sample(csoil_2100_4_INTERACT_domain$Value, 496))[c(1:2)])
  #csoil model_4_ INTERACT medians
  csoil_4_I_25_2020 <- as.numeric(quantile(csoil_2020_4_INTERACT_sites,probs = c(0.25),na.rm = T))
  csoil_4_I_50_2020 <- as.numeric(median(csoil_2020_4_INTERACT_sites,na.rm = T))
  csoil_4_I_75_2020 <- as.numeric(quantile(csoil_2020_4_INTERACT_sites,probs = c(0.75),na.rm = T))
  #csoil model_4_ INTERACT (without Russia) medians
  csoil_4_W_25_2020 <- as.numeric(quantile(csoil_2020_4_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  csoil_4_W_50_2020 <- as.numeric(median(csoil_2020_4_INTERACT_sites_war,na.rm = T))
  csoil_4_W_75_2020 <- as.numeric(quantile(csoil_2020_4_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #csoil model_4_ Domain medians & 50% CI 2020
  csoil_4_D_25_2020 <- as.numeric(quantile(sample(csoil_2020_4_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  csoil_4_D_50_2020 <- as.numeric(median(sample(csoil_2020_4_INTERACT_domain$Value, 496),na.rm = T))
  csoil_4_D_75_2020 <- as.numeric(quantile(sample(csoil_2020_4_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #csoil model_4_ Domain medians & 50% CI 2100
  csoil_4_D_25_2100 <- as.numeric(quantile(sample(csoil_2100_4_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  csoil_4_D_50_2100 <- as.numeric(median(sample(csoil_2100_4_INTERACT_domain$Value, 496),na.rm = T))
  csoil_4_D_75_2100 <- as.numeric(quantile(sample(csoil_2100_4_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  csoil_4[i,] <- round(c(csoil_4_ks_2020[1],csoil_4_ks_2020[2],csoil_4_ks_2020_war[1],csoil_4_ks_2020_war[2],
                         csoil_4_ks_2100[1],csoil_4_ks_2100[2],csoil_4_ks_2020_2100[1],csoil_4_ks_2020_2100[2],
                         csoil_4_I_25_2020,csoil_4_I_50_2020,csoil_4_I_75_2020,csoil_4_W_25_2020,csoil_4_W_50_2020,csoil_4_W_75_2020,
                         csoil_4_D_25_2020,csoil_4_D_50_2020,csoil_4_D_75_2020,csoil_4_D_25_2100,csoil_4_D_50_2100,csoil_4_D_75_2100),4)
  ##############################################################################################################################
  #csoil model_5_ KS test stats 
  csoil_5_ks_2020 <- as.numeric(ks.test(csoil_2020_5_INTERACT_sites, sample(csoil_2020_5_INTERACT_domain$Value, 496))[c(1:2)])
  csoil_5_ks_2020_war <- as.numeric(ks.test(csoil_2020_5_INTERACT_sites_war, sample(csoil_2020_5_INTERACT_domain$Value, 496))[c(1:2)])
  csoil_5_ks_2100 <- as.numeric(ks.test(csoil_2100_5_INTERACT_sites, sample(csoil_2100_5_INTERACT_domain$Value, 496))[c(1:2)])
  csoil_5_ks_2020_2100 <- as.numeric(ks.test(sample(csoil_2020_5_INTERACT_domain$Value, 496), sample(csoil_2100_5_INTERACT_domain$Value, 496))[c(1:2)])
  #csoil model_5_ INTERACT medians
  csoil_5_I_25_2020 <- as.numeric(quantile(csoil_2020_5_INTERACT_sites,probs = c(0.25),na.rm = T))
  csoil_5_I_50_2020 <- as.numeric(median(csoil_2020_5_INTERACT_sites,na.rm = T))
  csoil_5_I_75_2020 <- as.numeric(quantile(csoil_2020_5_INTERACT_sites,probs = c(0.75),na.rm = T))
  #csoil model_5_ INTERACT (without Russia) medians
  csoil_5_W_25_2020 <- as.numeric(quantile(csoil_2020_5_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  csoil_5_W_50_2020 <- as.numeric(median(csoil_2020_5_INTERACT_sites_war,na.rm = T))
  csoil_5_W_75_2020 <- as.numeric(quantile(csoil_2020_5_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #csoil model_5_ Domain medians & 50% CI 2020
  csoil_5_D_25_2020 <- as.numeric(quantile(sample(csoil_2020_5_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  csoil_5_D_50_2020 <- as.numeric(median(sample(csoil_2020_5_INTERACT_domain$Value, 496),na.rm = T))
  csoil_5_D_75_2020 <- as.numeric(quantile(sample(csoil_2020_5_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #csoil model_5_ Domain medians & 50% CI 2100
  csoil_5_D_25_2100 <- as.numeric(quantile(sample(csoil_2100_5_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  csoil_5_D_50_2100 <- as.numeric(median(sample(csoil_2100_5_INTERACT_domain$Value, 496),na.rm = T))
  csoil_5_D_75_2100 <- as.numeric(quantile(sample(csoil_2100_5_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  csoil_5[i,] <- round(c(csoil_5_ks_2020[1],csoil_5_ks_2020[2],csoil_5_ks_2020_war[1],csoil_5_ks_2020_war[2],
                         csoil_5_ks_2100[1],csoil_5_ks_2100[2],csoil_5_ks_2020_2100[1],csoil_5_ks_2020_2100[2],
                         csoil_5_I_25_2020,csoil_5_I_50_2020,csoil_5_I_75_2020,csoil_5_W_25_2020,csoil_5_W_50_2020,csoil_5_W_75_2020,
                         csoil_5_D_25_2020,csoil_5_D_50_2020,csoil_5_D_75_2020,csoil_5_D_25_2100,csoil_5_D_50_2100,csoil_5_D_75_2100),4)
  ##############################################################################################################################
  #csoil model_6_ KS test stats 
  csoil_6_ks_2020 <- as.numeric(ks.test(csoil_2020_6_INTERACT_sites, sample(csoil_2020_6_INTERACT_domain$Value, 496))[c(1:2)])
  csoil_6_ks_2020_war <- as.numeric(ks.test(csoil_2020_6_INTERACT_sites_war, sample(csoil_2020_6_INTERACT_domain$Value, 496))[c(1:2)])
  csoil_6_ks_2100 <- as.numeric(ks.test(csoil_2100_6_INTERACT_sites, sample(csoil_2100_6_INTERACT_domain$Value, 496))[c(1:2)])
  csoil_6_ks_2020_2100 <- as.numeric(ks.test(sample(csoil_2020_6_INTERACT_domain$Value, 496), sample(csoil_2100_6_INTERACT_domain$Value, 496))[c(1:2)])
  #csoil model_6_ INTERACT medians
  csoil_6_I_25_2020 <- as.numeric(quantile(csoil_2020_6_INTERACT_sites,probs = c(0.25),na.rm = T))
  csoil_6_I_50_2020 <- as.numeric(median(csoil_2020_6_INTERACT_sites,na.rm = T))
  csoil_6_I_75_2020 <- as.numeric(quantile(csoil_2020_6_INTERACT_sites,probs = c(0.75),na.rm = T))
  #csoil model_6_ INTERACT (without Russia) medians
  csoil_6_W_25_2020 <- as.numeric(quantile(csoil_2020_6_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  csoil_6_W_50_2020 <- as.numeric(median(csoil_2020_6_INTERACT_sites_war,na.rm = T))
  csoil_6_W_75_2020 <- as.numeric(quantile(csoil_2020_6_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #csoil model_6_ Domain medians & 50% CI 2020
  csoil_6_D_25_2020 <- as.numeric(quantile(sample(csoil_2020_6_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  csoil_6_D_50_2020 <- as.numeric(median(sample(csoil_2020_6_INTERACT_domain$Value, 496),na.rm = T))
  csoil_6_D_75_2020 <- as.numeric(quantile(sample(csoil_2020_6_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #csoil model_6_ Domain medians & 50% CI 2100
  csoil_6_D_25_2100 <- as.numeric(quantile(sample(csoil_2100_6_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  csoil_6_D_50_2100 <- as.numeric(median(sample(csoil_2100_6_INTERACT_domain$Value, 496),na.rm = T))
  csoil_6_D_75_2100 <- as.numeric(quantile(sample(csoil_2100_6_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  csoil_6[i,] <- round(c(csoil_6_ks_2020[1],csoil_6_ks_2020[2],csoil_6_ks_2020_war[1],csoil_6_ks_2020_war[2],
                         csoil_6_ks_2100[1],csoil_6_ks_2100[2],csoil_6_ks_2020_2100[1],csoil_6_ks_2020_2100[2],
                         csoil_6_I_25_2020,csoil_6_I_50_2020,csoil_6_I_75_2020,csoil_6_W_25_2020,csoil_6_W_50_2020,csoil_6_W_75_2020,
                         csoil_6_D_25_2020,csoil_6_D_50_2020,csoil_6_D_75_2020,csoil_6_D_25_2100,csoil_6_D_50_2100,csoil_6_D_75_2100),4)
  ##############################################################################################################################
  #csoil model_7_ KS test stats 
  csoil_7_ks_2020 <- as.numeric(ks.test(csoil_2020_7_INTERACT_sites, sample(csoil_2020_7_INTERACT_domain$Value, 496))[c(1:2)])
  csoil_7_ks_2020_war <- as.numeric(ks.test(csoil_2020_7_INTERACT_sites_war, sample(csoil_2020_7_INTERACT_domain$Value, 496))[c(1:2)])
  csoil_7_ks_2100 <- as.numeric(ks.test(csoil_2100_7_INTERACT_sites, sample(csoil_2100_7_INTERACT_domain$Value, 496))[c(1:2)])
  csoil_7_ks_2020_2100 <- as.numeric(ks.test(sample(csoil_2020_7_INTERACT_domain$Value, 496), sample(csoil_2100_7_INTERACT_domain$Value, 496))[c(1:2)])
  #csoil model_7_ INTERACT medians
  csoil_7_I_25_2020 <- as.numeric(quantile(csoil_2020_7_INTERACT_sites,probs = c(0.25),na.rm = T))
  csoil_7_I_50_2020 <- as.numeric(median(csoil_2020_7_INTERACT_sites,na.rm = T))
  csoil_7_I_75_2020 <- as.numeric(quantile(csoil_2020_7_INTERACT_sites,probs = c(0.75),na.rm = T))
  #csoil model_7_ INTERACT (without Russia) medians
  csoil_7_W_25_2020 <- as.numeric(quantile(csoil_2020_7_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  csoil_7_W_50_2020 <- as.numeric(median(csoil_2020_7_INTERACT_sites_war,na.rm = T))
  csoil_7_W_75_2020 <- as.numeric(quantile(csoil_2020_7_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #csoil model_7_ Domain medians & 50% CI 2020
  csoil_7_D_25_2020 <- as.numeric(quantile(sample(csoil_2020_7_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  csoil_7_D_50_2020 <- as.numeric(median(sample(csoil_2020_7_INTERACT_domain$Value, 496),na.rm = T))
  csoil_7_D_75_2020 <- as.numeric(quantile(sample(csoil_2020_7_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #csoil model_7_ Domain medians & 50% CI 2100
  csoil_7_D_25_2100 <- as.numeric(quantile(sample(csoil_2100_7_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  csoil_7_D_50_2100 <- as.numeric(median(sample(csoil_2100_7_INTERACT_domain$Value, 496),na.rm = T))
  csoil_7_D_75_2100 <- as.numeric(quantile(sample(csoil_2100_7_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  csoil_7[i,] <- round(c(csoil_7_ks_2020[1],csoil_7_ks_2020[2],csoil_7_ks_2020_war[1],csoil_7_ks_2020_war[2],
                         csoil_7_ks_2100[1],csoil_7_ks_2100[2],csoil_7_ks_2020_2100[1],csoil_7_ks_2020_2100[2],
                         csoil_7_I_25_2020,csoil_7_I_50_2020,csoil_7_I_75_2020,csoil_7_W_25_2020,csoil_7_W_50_2020,csoil_7_W_75_2020,
                         csoil_7_D_25_2020,csoil_7_D_50_2020,csoil_7_D_75_2020,csoil_7_D_25_2100,csoil_7_D_50_2100,csoil_7_D_75_2100),4)
  
  ##################################################################################################################################
  #csoil model_8_ KS test stats 
  csoil_8_ks_2020 <- as.numeric(ks.test(csoil_2020_8_INTERACT_sites, sample(csoil_2020_8_INTERACT_domain$Value, 496))[c(1:2)])
  csoil_8_ks_2020_war <- as.numeric(ks.test(csoil_2020_8_INTERACT_sites_war, sample(csoil_2020_8_INTERACT_domain$Value, 496))[c(1:2)])
  csoil_8_ks_2100 <- as.numeric(ks.test(csoil_2100_8_INTERACT_sites, sample(csoil_2100_8_INTERACT_domain$Value, 496))[c(1:2)])
  csoil_8_ks_2020_2100 <- as.numeric(ks.test(sample(csoil_2020_8_INTERACT_domain$Value, 496), sample(csoil_2100_8_INTERACT_domain$Value, 496))[c(1:2)])
  #csoil model_8_ INTERACT medians
  csoil_8_I_25_2020 <- as.numeric(quantile(csoil_2020_8_INTERACT_sites,probs = c(0.25),na.rm = T))
  csoil_8_I_50_2020 <- as.numeric(median(csoil_2020_8_INTERACT_sites,na.rm = T))
  csoil_8_I_75_2020 <- as.numeric(quantile(csoil_2020_8_INTERACT_sites,probs = c(0.75),na.rm = T))
  #csoil model_8_ INTERACT (without Russia) medians
  csoil_8_W_25_2020 <- as.numeric(quantile(csoil_2020_8_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  csoil_8_W_50_2020 <- as.numeric(median(csoil_2020_8_INTERACT_sites_war,na.rm = T))
  csoil_8_W_75_2020 <- as.numeric(quantile(csoil_2020_8_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #csoil model_8_ Domain medians & 50% CI 2020
  csoil_8_D_25_2020 <- as.numeric(quantile(sample(csoil_2020_8_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  csoil_8_D_50_2020 <- as.numeric(median(sample(csoil_2020_8_INTERACT_domain$Value, 496),na.rm = T))
  csoil_8_D_75_2020 <- as.numeric(quantile(sample(csoil_2020_8_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #csoil model_8_ Domain medians & 50% CI 2100
  csoil_8_D_25_2100 <- as.numeric(quantile(sample(csoil_2100_8_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  csoil_8_D_50_2100 <- as.numeric(median(sample(csoil_2100_8_INTERACT_domain$Value, 496),na.rm = T))
  csoil_8_D_75_2100 <- as.numeric(quantile(sample(csoil_2100_8_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  csoil_8[i,] <- round(c(csoil_8_ks_2020[1],csoil_8_ks_2020[2],csoil_8_ks_2020_war[1],csoil_8_ks_2020_war[2],
                         csoil_8_ks_2100[1],csoil_8_ks_2100[2],csoil_8_ks_2020_2100[1],csoil_8_ks_2020_2100[2],
                         csoil_8_I_25_2020,csoil_8_I_50_2020,csoil_8_I_75_2020,csoil_8_W_25_2020,csoil_8_W_50_2020,csoil_8_W_75_2020,
                         csoil_8_D_25_2020,csoil_8_D_50_2020,csoil_8_D_75_2020,csoil_8_D_25_2100,csoil_8_D_50_2100,csoil_8_D_75_2100),4)
  ##################################################################################################################################
  #npp model_1_ KS test stats 
  npp_1_ks_2020 <- as.numeric(ks.test(npp_2020_1_INTERACT_sites, sample(npp_2020_1_INTERACT_domain$Value, 496))[c(1:2)])
  npp_1_ks_2020_war <- as.numeric(ks.test(npp_2020_1_INTERACT_sites_war, sample(npp_2020_1_INTERACT_domain$Value, 496))[c(1:2)])
  npp_1_ks_2100 <- as.numeric(ks.test(npp_2100_1_INTERACT_sites, sample(npp_2100_1_INTERACT_domain$Value, 496))[c(1:2)])
  npp_1_ks_2020_2100 <- as.numeric(ks.test(sample(npp_2020_1_INTERACT_domain$Value, 496), sample(npp_2100_1_INTERACT_domain$Value, 496))[c(1:2)])
  #npp model_1_ INTERACT medians
  npp_1_I_25_2020 <- as.numeric(quantile(npp_2020_1_INTERACT_sites,probs = c(0.25),na.rm = T))
  npp_1_I_50_2020 <- as.numeric(median(npp_2020_1_INTERACT_sites,na.rm = T))
  npp_1_I_75_2020 <- as.numeric(quantile(npp_2020_1_INTERACT_sites,probs = c(0.75),na.rm = T))
  #npp model_1_ INTERACT (without Russia) medians
  npp_1_W_25_2020 <- as.numeric(quantile(npp_2020_1_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  npp_1_W_50_2020 <- as.numeric(median(npp_2020_1_INTERACT_sites_war,na.rm = T))
  npp_1_W_75_2020 <- as.numeric(quantile(npp_2020_1_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #npp model_1_ Domain medians & 50% CI 2020
  npp_1_D_25_2020 <- as.numeric(quantile(sample(npp_2020_1_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  npp_1_D_50_2020 <- as.numeric(median(sample(npp_2020_1_INTERACT_domain$Value, 496),na.rm = T))
  npp_1_D_75_2020 <- as.numeric(quantile(sample(npp_2020_1_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #npp model_1_ Domain medians & 50% CI 2100
  npp_1_D_25_2100 <- as.numeric(quantile(sample(npp_2100_1_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  npp_1_D_50_2100 <- as.numeric(median(sample(npp_2100_1_INTERACT_domain$Value, 496),na.rm = T))
  npp_1_D_75_2100 <- as.numeric(quantile(sample(npp_2100_1_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  npp_1[i,] <- round(c(npp_1_ks_2020[1],npp_1_ks_2020[2],npp_1_ks_2020_war[1],npp_1_ks_2020_war[2],
                       npp_1_ks_2100[1],npp_1_ks_2100[2],npp_1_ks_2020_2100[1],npp_1_ks_2020_2100[2],
                       npp_1_I_25_2020,npp_1_I_50_2020,npp_1_I_75_2020,npp_1_W_25_2020,npp_1_W_50_2020,npp_1_W_75_2020,
                       npp_1_D_25_2020,npp_1_D_50_2020,npp_1_D_75_2020,npp_1_D_25_2100,npp_1_D_50_2100,npp_1_D_75_2100),4)
  ##################################################################################################################################
  #npp model_2_ KS test stats 
  npp_2_ks_2020 <- as.numeric(ks.test(npp_2020_2_INTERACT_sites, sample(npp_2020_2_INTERACT_domain$Value, 496))[c(1:2)])
  npp_2_ks_2020_war <- as.numeric(ks.test(npp_2020_2_INTERACT_sites_war, sample(npp_2020_2_INTERACT_domain$Value, 496))[c(1:2)])
  npp_2_ks_2100 <- as.numeric(ks.test(npp_2100_2_INTERACT_sites, sample(npp_2100_2_INTERACT_domain$Value, 496))[c(1:2)])
  npp_2_ks_2020_2100 <- as.numeric(ks.test(sample(npp_2020_2_INTERACT_domain$Value, 496), sample(npp_2100_2_INTERACT_domain$Value, 496))[c(1:2)])
  #npp model_2_ INTERACT medians
  npp_2_I_25_2020 <- as.numeric(quantile(npp_2020_2_INTERACT_sites,probs = c(0.25),na.rm = T))
  npp_2_I_50_2020 <- as.numeric(median(npp_2020_2_INTERACT_sites,na.rm = T))
  npp_2_I_75_2020 <- as.numeric(quantile(npp_2020_2_INTERACT_sites,probs = c(0.75),na.rm = T))
  #npp model_2_ INTERACT (without Russia) medians
  npp_2_W_25_2020 <- as.numeric(quantile(npp_2020_2_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  npp_2_W_50_2020 <- as.numeric(median(npp_2020_2_INTERACT_sites_war,na.rm = T))
  npp_2_W_75_2020 <- as.numeric(quantile(npp_2020_2_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #npp model_2_ Domain medians & 50% CI 2020
  npp_2_D_25_2020 <- as.numeric(quantile(sample(npp_2020_2_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  npp_2_D_50_2020 <- as.numeric(median(sample(npp_2020_2_INTERACT_domain$Value, 496),na.rm = T))
  npp_2_D_75_2020 <- as.numeric(quantile(sample(npp_2020_2_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #npp model_2_ Domain medians & 50% CI 2100
  npp_2_D_25_2100 <- as.numeric(quantile(sample(npp_2100_2_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  npp_2_D_50_2100 <- as.numeric(median(sample(npp_2100_2_INTERACT_domain$Value, 496),na.rm = T))
  npp_2_D_75_2100 <- as.numeric(quantile(sample(npp_2100_2_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  npp_2[i,] <- round(c(npp_2_ks_2020[1],npp_2_ks_2020[2],npp_2_ks_2020_war[1],npp_2_ks_2020_war[2],
                       npp_2_ks_2100[1],npp_2_ks_2100[2],npp_2_ks_2020_2100[1],npp_2_ks_2020_2100[2],
                       npp_2_I_25_2020,npp_2_I_50_2020,npp_2_I_75_2020,npp_2_W_25_2020,npp_2_W_50_2020,npp_2_W_75_2020,
                       npp_2_D_25_2020,npp_2_D_50_2020,npp_2_D_75_2020,npp_2_D_25_2100,npp_2_D_50_2100,npp_2_D_75_2100),4)
  ##############################################################################################################################
  #npp model_3_ KS test stats 
  npp_3_ks_2020 <- as.numeric(ks.test(npp_2020_3_INTERACT_sites, sample(npp_2020_3_INTERACT_domain$Value, 496))[c(1:2)])
  npp_3_ks_2020_war <- as.numeric(ks.test(npp_2020_3_INTERACT_sites_war, sample(npp_2020_3_INTERACT_domain$Value, 496))[c(1:2)])
  npp_3_ks_2100 <- as.numeric(ks.test(npp_2100_3_INTERACT_sites, sample(npp_2100_3_INTERACT_domain$Value, 496))[c(1:2)])
  npp_3_ks_2020_2100 <- as.numeric(ks.test(sample(npp_2020_3_INTERACT_domain$Value, 496), sample(npp_2100_3_INTERACT_domain$Value, 496))[c(1:2)])
  #npp model_3_ INTERACT medians
  npp_3_I_25_2020 <- as.numeric(quantile(npp_2020_3_INTERACT_sites,probs = c(0.25),na.rm = T))
  npp_3_I_50_2020 <- as.numeric(median(npp_2020_3_INTERACT_sites,na.rm = T))
  npp_3_I_75_2020 <- as.numeric(quantile(npp_2020_3_INTERACT_sites,probs = c(0.75),na.rm = T))
  #npp model_3_ INTERACT (without Russia) medians
  npp_3_W_25_2020 <- as.numeric(quantile(npp_2020_3_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  npp_3_W_50_2020 <- as.numeric(median(npp_2020_3_INTERACT_sites_war,na.rm = T))
  npp_3_W_75_2020 <- as.numeric(quantile(npp_2020_3_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #npp model_3_ Domain medians & 50% CI 2020
  npp_3_D_25_2020 <- as.numeric(quantile(sample(npp_2020_3_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  npp_3_D_50_2020 <- as.numeric(median(sample(npp_2020_3_INTERACT_domain$Value, 496),na.rm = T))
  npp_3_D_75_2020 <- as.numeric(quantile(sample(npp_2020_3_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #npp model_3_ Domain medians & 50% CI 2100
  npp_3_D_25_2100 <- as.numeric(quantile(sample(npp_2100_3_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  npp_3_D_50_2100 <- as.numeric(median(sample(npp_2100_3_INTERACT_domain$Value, 496),na.rm = T))
  npp_3_D_75_2100 <- as.numeric(quantile(sample(npp_2100_3_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  npp_3[i,] <- round(c(npp_3_ks_2020[1],npp_3_ks_2020[2],npp_3_ks_2020_war[1],npp_3_ks_2020_war[2],
                       npp_3_ks_2100[1],npp_3_ks_2100[2],npp_3_ks_2020_2100[1],npp_3_ks_2020_2100[2],
                       npp_3_I_25_2020,npp_3_I_50_2020,npp_3_I_75_2020,npp_3_W_25_2020,npp_3_W_50_2020,npp_3_W_75_2020,
                       npp_3_D_25_2020,npp_3_D_50_2020,npp_3_D_75_2020,npp_3_D_25_2100,npp_3_D_50_2100,npp_3_D_75_2100),4)
  ##############################################################################################################################
  #npp model_4_ KS test stats 
  npp_4_ks_2020 <- as.numeric(ks.test(npp_2020_4_INTERACT_sites, sample(npp_2020_4_INTERACT_domain$Value, 496))[c(1:2)])
  npp_4_ks_2020_war <- as.numeric(ks.test(npp_2020_4_INTERACT_sites_war, sample(npp_2020_4_INTERACT_domain$Value, 496))[c(1:2)])
  npp_4_ks_2100 <- as.numeric(ks.test(npp_2100_4_INTERACT_sites, sample(npp_2100_4_INTERACT_domain$Value, 496))[c(1:2)])
  npp_4_ks_2020_2100 <- as.numeric(ks.test(sample(npp_2020_4_INTERACT_domain$Value, 496), sample(npp_2100_4_INTERACT_domain$Value, 496))[c(1:2)])
  #npp model_4_ INTERACT medians
  npp_4_I_25_2020 <- as.numeric(quantile(npp_2020_4_INTERACT_sites,probs = c(0.25),na.rm = T))
  npp_4_I_50_2020 <- as.numeric(median(npp_2020_4_INTERACT_sites,na.rm = T))
  npp_4_I_75_2020 <- as.numeric(quantile(npp_2020_4_INTERACT_sites,probs = c(0.75),na.rm = T))
  #npp model_4_ INTERACT (without Russia) medians
  npp_4_W_25_2020 <- as.numeric(quantile(npp_2020_4_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  npp_4_W_50_2020 <- as.numeric(median(npp_2020_4_INTERACT_sites_war,na.rm = T))
  npp_4_W_75_2020 <- as.numeric(quantile(npp_2020_4_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #npp model_4_ Domain medians & 50% CI 2020
  npp_4_D_25_2020 <- as.numeric(quantile(sample(npp_2020_4_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  npp_4_D_50_2020 <- as.numeric(median(sample(npp_2020_4_INTERACT_domain$Value, 496),na.rm = T))
  npp_4_D_75_2020 <- as.numeric(quantile(sample(npp_2020_4_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #npp model_4_ Domain medians & 50% CI 2100
  npp_4_D_25_2100 <- as.numeric(quantile(sample(npp_2100_4_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  npp_4_D_50_2100 <- as.numeric(median(sample(npp_2100_4_INTERACT_domain$Value, 496),na.rm = T))
  npp_4_D_75_2100 <- as.numeric(quantile(sample(npp_2100_4_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  npp_4[i,] <- round(c(npp_4_ks_2020[1],npp_4_ks_2020[2],npp_4_ks_2020_war[1],npp_4_ks_2020_war[2],
                       npp_4_ks_2100[1],npp_4_ks_2100[2],npp_4_ks_2020_2100[1],npp_4_ks_2020_2100[2],
                       npp_4_I_25_2020,npp_4_I_50_2020,npp_4_I_75_2020,npp_4_W_25_2020,npp_4_W_50_2020,npp_4_W_75_2020,
                       npp_4_D_25_2020,npp_4_D_50_2020,npp_4_D_75_2020,npp_4_D_25_2100,npp_4_D_50_2100,npp_4_D_75_2100),4)
  ##############################################################################################################################
  #npp model_5_ KS test stats 
  npp_5_ks_2020 <- as.numeric(ks.test(npp_2020_5_INTERACT_sites, sample(npp_2020_5_INTERACT_domain$Value, 496))[c(1:2)])
  npp_5_ks_2020_war <- as.numeric(ks.test(npp_2020_5_INTERACT_sites_war, sample(npp_2020_5_INTERACT_domain$Value, 496))[c(1:2)])
  npp_5_ks_2100 <- as.numeric(ks.test(npp_2100_5_INTERACT_sites, sample(npp_2100_5_INTERACT_domain$Value, 496))[c(1:2)])
  npp_5_ks_2020_2100 <- as.numeric(ks.test(sample(npp_2020_5_INTERACT_domain$Value, 496), sample(npp_2100_5_INTERACT_domain$Value, 496))[c(1:2)])
  #npp model_5_ INTERACT medians
  npp_5_I_25_2020 <- as.numeric(quantile(npp_2020_5_INTERACT_sites,probs = c(0.25),na.rm = T))
  npp_5_I_50_2020 <- as.numeric(median(npp_2020_5_INTERACT_sites,na.rm = T))
  npp_5_I_75_2020 <- as.numeric(quantile(npp_2020_5_INTERACT_sites,probs = c(0.75),na.rm = T))
  #npp model_5_ INTERACT (without Russia) medians
  npp_5_W_25_2020 <- as.numeric(quantile(npp_2020_5_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  npp_5_W_50_2020 <- as.numeric(median(npp_2020_5_INTERACT_sites_war,na.rm = T))
  npp_5_W_75_2020 <- as.numeric(quantile(npp_2020_5_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #npp model_5_ Domain medians & 50% CI 2020
  npp_5_D_25_2020 <- as.numeric(quantile(sample(npp_2020_5_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  npp_5_D_50_2020 <- as.numeric(median(sample(npp_2020_5_INTERACT_domain$Value, 496),na.rm = T))
  npp_5_D_75_2020 <- as.numeric(quantile(sample(npp_2020_5_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #npp model_5_ Domain medians & 50% CI 2100
  npp_5_D_25_2100 <- as.numeric(quantile(sample(npp_2100_5_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  npp_5_D_50_2100 <- as.numeric(median(sample(npp_2100_5_INTERACT_domain$Value, 496),na.rm = T))
  npp_5_D_75_2100 <- as.numeric(quantile(sample(npp_2100_5_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  npp_5[i,] <- round(c(npp_5_ks_2020[1],npp_5_ks_2020[2],npp_5_ks_2020_war[1],npp_5_ks_2020_war[2],
                       npp_5_ks_2100[1],npp_5_ks_2100[2],npp_5_ks_2020_2100[1],npp_5_ks_2020_2100[2],
                       npp_5_I_25_2020,npp_5_I_50_2020,npp_5_I_75_2020,npp_5_W_25_2020,npp_5_W_50_2020,npp_5_W_75_2020,
                       npp_5_D_25_2020,npp_5_D_50_2020,npp_5_D_75_2020,npp_5_D_25_2100,npp_5_D_50_2100,npp_5_D_75_2100),4)
  ##############################################################################################################################
  #npp model_6_ KS test stats 
  npp_6_ks_2020 <- as.numeric(ks.test(npp_2020_6_INTERACT_sites, sample(npp_2020_6_INTERACT_domain$Value, 496))[c(1:2)])
  npp_6_ks_2020_war <- as.numeric(ks.test(npp_2020_6_INTERACT_sites_war, sample(npp_2020_6_INTERACT_domain$Value, 496))[c(1:2)])
  npp_6_ks_2100 <- as.numeric(ks.test(npp_2100_6_INTERACT_sites, sample(npp_2100_6_INTERACT_domain$Value, 496))[c(1:2)])
  npp_6_ks_2020_2100 <- as.numeric(ks.test(sample(npp_2020_6_INTERACT_domain$Value, 496), sample(npp_2100_6_INTERACT_domain$Value, 496))[c(1:2)])
  #npp model_6_ INTERACT medians
  npp_6_I_25_2020 <- as.numeric(quantile(npp_2020_6_INTERACT_sites,probs = c(0.25),na.rm = T))
  npp_6_I_50_2020 <- as.numeric(median(npp_2020_6_INTERACT_sites,na.rm = T))
  npp_6_I_75_2020 <- as.numeric(quantile(npp_2020_6_INTERACT_sites,probs = c(0.75),na.rm = T))
  #npp model_6_ INTERACT (without Russia) medians
  npp_6_W_25_2020 <- as.numeric(quantile(npp_2020_6_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  npp_6_W_50_2020 <- as.numeric(median(npp_2020_6_INTERACT_sites_war,na.rm = T))
  npp_6_W_75_2020 <- as.numeric(quantile(npp_2020_6_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #npp model_6_ Domain medians & 50% CI 2020
  npp_6_D_25_2020 <- as.numeric(quantile(sample(npp_2020_6_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  npp_6_D_50_2020 <- as.numeric(median(sample(npp_2020_6_INTERACT_domain$Value, 496),na.rm = T))
  npp_6_D_75_2020 <- as.numeric(quantile(sample(npp_2020_6_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #npp model_6_ Domain medians & 50% CI 2100
  npp_6_D_25_2100 <- as.numeric(quantile(sample(npp_2100_6_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  npp_6_D_50_2100 <- as.numeric(median(sample(npp_2100_6_INTERACT_domain$Value, 496),na.rm = T))
  npp_6_D_75_2100 <- as.numeric(quantile(sample(npp_2100_6_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  npp_6[i,] <- round(c(npp_6_ks_2020[1],npp_6_ks_2020[2],npp_6_ks_2020_war[1],npp_6_ks_2020_war[2],
                       npp_6_ks_2100[1],npp_6_ks_2100[2],npp_6_ks_2020_2100[1],npp_6_ks_2020_2100[2],
                       npp_6_I_25_2020,npp_6_I_50_2020,npp_6_I_75_2020,npp_6_W_25_2020,npp_6_W_50_2020,npp_6_W_75_2020,
                       npp_6_D_25_2020,npp_6_D_50_2020,npp_6_D_75_2020,npp_6_D_25_2100,npp_6_D_50_2100,npp_6_D_75_2100),4)
  ##############################################################################################################################
  #npp model_7_ KS test stats 
  npp_7_ks_2020 <- as.numeric(ks.test(npp_2020_7_INTERACT_sites, sample(npp_2020_7_INTERACT_domain$Value, 488))[c(1:2)])
  npp_7_ks_2020_war <- as.numeric(ks.test(npp_2020_7_INTERACT_sites_war, sample(npp_2020_7_INTERACT_domain$Value, 488))[c(1:2)])
  npp_7_ks_2100 <- as.numeric(ks.test(npp_2100_7_INTERACT_sites, sample(npp_2100_7_INTERACT_domain$Value, 488))[c(1:2)])
  npp_7_ks_2020_2100 <- as.numeric(ks.test(sample(npp_2020_7_INTERACT_domain$Value, 488), sample(npp_2100_7_INTERACT_domain$Value, 488))[c(1:2)])
  #npp model_7_ INTERACT medians
  npp_7_I_25_2020 <- as.numeric(quantile(npp_2020_7_INTERACT_sites,probs = c(0.25),na.rm = T))
  npp_7_I_50_2020 <- as.numeric(median(npp_2020_7_INTERACT_sites,na.rm = T))
  npp_7_I_75_2020 <- as.numeric(quantile(npp_2020_7_INTERACT_sites,probs = c(0.75),na.rm = T))
  #npp model_7_ INTERACT (without Russia) medians
  npp_7_W_25_2020 <- as.numeric(quantile(npp_2020_7_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  npp_7_W_50_2020 <- as.numeric(median(npp_2020_7_INTERACT_sites_war,na.rm = T))
  npp_7_W_75_2020 <- as.numeric(quantile(npp_2020_7_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #npp model_7_ Domain medians & 50% CI 2020
  npp_7_D_25_2020 <- as.numeric(quantile(sample(npp_2020_7_INTERACT_domain$Value, 488),probs = c(0.25),na.rm = T))
  npp_7_D_50_2020 <- as.numeric(median(sample(npp_2020_7_INTERACT_domain$Value, 488),na.rm = T))
  npp_7_D_75_2020 <- as.numeric(quantile(sample(npp_2020_7_INTERACT_domain$Value, 488),probs = c(0.75),na.rm = T))
  #npp model_7_ Domain medians & 50% CI 2100
  npp_7_D_25_2100 <- as.numeric(quantile(sample(npp_2100_7_INTERACT_domain$Value, 488),probs = c(0.25),na.rm = T))
  npp_7_D_50_2100 <- as.numeric(median(sample(npp_2100_7_INTERACT_domain$Value, 488),na.rm = T))
  npp_7_D_75_2100 <- as.numeric(quantile(sample(npp_2100_7_INTERACT_domain$Value, 488),probs = c(0.75),na.rm = T))
  #add all to same df
  npp_7[i,] <- round(c(npp_7_ks_2020[1],npp_7_ks_2020[2],npp_7_ks_2020_war[1],npp_7_ks_2020_war[2],
                       npp_7_ks_2100[1],npp_7_ks_2100[2],npp_7_ks_2020_2100[1],npp_7_ks_2020_2100[2],
                       npp_7_I_25_2020,npp_7_I_50_2020,npp_7_I_75_2020,npp_7_W_25_2020,npp_7_W_50_2020,npp_7_W_75_2020,
                       npp_7_D_25_2020,npp_7_D_50_2020,npp_7_D_75_2020,npp_7_D_25_2100,npp_7_D_50_2100,npp_7_D_75_2100),4)
  
  ##################################################################################################################################
  #npp model_8_ KS test stats 
  npp_8_ks_2020 <- as.numeric(ks.test(npp_2020_8_INTERACT_sites, sample(npp_2020_8_INTERACT_domain$Value, 496))[c(1:2)])
  npp_8_ks_2020_war <- as.numeric(ks.test(npp_2020_8_INTERACT_sites_war, sample(npp_2020_8_INTERACT_domain$Value, 496))[c(1:2)])
  npp_8_ks_2100 <- as.numeric(ks.test(npp_2100_8_INTERACT_sites, sample(npp_2100_8_INTERACT_domain$Value, 496))[c(1:2)])
  npp_8_ks_2020_2100 <- as.numeric(ks.test(sample(npp_2020_8_INTERACT_domain$Value, 496), sample(npp_2100_8_INTERACT_domain$Value, 496))[c(1:2)])
  #npp model_8_ INTERACT medians
  npp_8_I_25_2020 <- as.numeric(quantile(npp_2020_8_INTERACT_sites,probs = c(0.25),na.rm = T))
  npp_8_I_50_2020 <- as.numeric(median(npp_2020_8_INTERACT_sites,na.rm = T))
  npp_8_I_75_2020 <- as.numeric(quantile(npp_2020_8_INTERACT_sites,probs = c(0.75),na.rm = T))
  #npp model_8_ INTERACT (without Russia) medians
  npp_8_W_25_2020 <- as.numeric(quantile(npp_2020_8_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  npp_8_W_50_2020 <- as.numeric(median(npp_2020_8_INTERACT_sites_war,na.rm = T))
  npp_8_W_75_2020 <- as.numeric(quantile(npp_2020_8_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #npp model_8_ Domain medians & 50% CI 2020
  npp_8_D_25_2020 <- as.numeric(quantile(sample(npp_2020_8_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  npp_8_D_50_2020 <- as.numeric(median(sample(npp_2020_8_INTERACT_domain$Value, 496),na.rm = T))
  npp_8_D_75_2020 <- as.numeric(quantile(sample(npp_2020_8_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #npp model_8_ Domain medians & 50% CI 2100
  npp_8_D_25_2100 <- as.numeric(quantile(sample(npp_2100_8_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  npp_8_D_50_2100 <- as.numeric(median(sample(npp_2100_8_INTERACT_domain$Value, 496),na.rm = T))
  npp_8_D_75_2100 <- as.numeric(quantile(sample(npp_2100_8_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  npp_8[i,] <- round(c(npp_8_ks_2020[1],npp_8_ks_2020[2],npp_8_ks_2020_war[1],npp_8_ks_2020_war[2],
                       npp_8_ks_2100[1],npp_8_ks_2100[2],npp_8_ks_2020_2100[1],npp_8_ks_2020_2100[2],
                       npp_8_I_25_2020,npp_8_I_50_2020,npp_8_I_75_2020,npp_8_W_25_2020,npp_8_W_50_2020,npp_8_W_75_2020,
                       npp_8_D_25_2020,npp_8_D_50_2020,npp_8_D_75_2020,npp_8_D_25_2100,npp_8_D_50_2100,npp_8_D_75_2100),4)
  ##################################################################################################################################
  #rh model_1_ KS test stats 
  rh_1_ks_2020 <- as.numeric(ks.test(rh_2020_1_INTERACT_sites, sample(rh_2020_1_INTERACT_domain$Value, 496))[c(1:2)])
  rh_1_ks_2020_war <- as.numeric(ks.test(rh_2020_1_INTERACT_sites_war, sample(rh_2020_1_INTERACT_domain$Value, 496))[c(1:2)])
  rh_1_ks_2100 <- as.numeric(ks.test(rh_2100_1_INTERACT_sites, sample(rh_2100_1_INTERACT_domain$Value, 496))[c(1:2)])
  rh_1_ks_2020_2100 <- as.numeric(ks.test(sample(rh_2020_1_INTERACT_domain$Value, 496), sample(rh_2100_1_INTERACT_domain$Value, 496))[c(1:2)])
  #rh model_1_ INTERACT medians
  rh_1_I_25_2020 <- as.numeric(quantile(rh_2020_1_INTERACT_sites,probs = c(0.25),na.rm = T))
  rh_1_I_50_2020 <- as.numeric(median(rh_2020_1_INTERACT_sites,na.rm = T))
  rh_1_I_75_2020 <- as.numeric(quantile(rh_2020_1_INTERACT_sites,probs = c(0.75),na.rm = T))
  #rh model_1_ INTERACT (without Russia) medians
  rh_1_W_25_2020 <- as.numeric(quantile(rh_2020_1_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  rh_1_W_50_2020 <- as.numeric(median(rh_2020_1_INTERACT_sites_war,na.rm = T))
  rh_1_W_75_2020 <- as.numeric(quantile(rh_2020_1_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #rh model_1_ Domain medians & 50% CI 2020
  rh_1_D_25_2020 <- as.numeric(quantile(sample(rh_2020_1_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  rh_1_D_50_2020 <- as.numeric(median(sample(rh_2020_1_INTERACT_domain$Value, 496),na.rm = T))
  rh_1_D_75_2020 <- as.numeric(quantile(sample(rh_2020_1_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #rh model_1_ Domain medians & 50% CI 2100
  rh_1_D_25_2100 <- as.numeric(quantile(sample(rh_2100_1_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  rh_1_D_50_2100 <- as.numeric(median(sample(rh_2100_1_INTERACT_domain$Value, 496),na.rm = T))
  rh_1_D_75_2100 <- as.numeric(quantile(sample(rh_2100_1_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  rh_1[i,] <- round(c(rh_1_ks_2020[1],rh_1_ks_2020[2],rh_1_ks_2020_war[1],rh_1_ks_2020_war[2],
                      rh_1_ks_2100[1],rh_1_ks_2100[2],rh_1_ks_2020_2100[1],rh_1_ks_2020_2100[2],
                      rh_1_I_25_2020,rh_1_I_50_2020,rh_1_I_75_2020,rh_1_W_25_2020,rh_1_W_50_2020,rh_1_W_75_2020,
                      rh_1_D_25_2020,rh_1_D_50_2020,rh_1_D_75_2020,rh_1_D_25_2100,rh_1_D_50_2100,rh_1_D_75_2100),4)
  ##################################################################################################################################
  #rh model_2_ KS test stats 
  rh_2_ks_2020 <- as.numeric(ks.test(rh_2020_2_INTERACT_sites, sample(rh_2020_2_INTERACT_domain$Value, 496))[c(1:2)])
  rh_2_ks_2020_war <- as.numeric(ks.test(rh_2020_2_INTERACT_sites_war, sample(rh_2020_2_INTERACT_domain$Value, 496))[c(1:2)])
  rh_2_ks_2100 <- as.numeric(ks.test(rh_2100_2_INTERACT_sites, sample(rh_2100_2_INTERACT_domain$Value, 496))[c(1:2)])
  rh_2_ks_2020_2100 <- as.numeric(ks.test(sample(rh_2020_2_INTERACT_domain$Value, 496), sample(rh_2100_2_INTERACT_domain$Value, 496))[c(1:2)])
  #rh model_2_ INTERACT medians
  rh_2_I_25_2020 <- as.numeric(quantile(rh_2020_2_INTERACT_sites,probs = c(0.25),na.rm = T))
  rh_2_I_50_2020 <- as.numeric(median(rh_2020_2_INTERACT_sites,na.rm = T))
  rh_2_I_75_2020 <- as.numeric(quantile(rh_2020_2_INTERACT_sites,probs = c(0.75),na.rm = T))
  #rh model_2_ INTERACT (without Russia) medians
  rh_2_W_25_2020 <- as.numeric(quantile(rh_2020_2_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  rh_2_W_50_2020 <- as.numeric(median(rh_2020_2_INTERACT_sites_war,na.rm = T))
  rh_2_W_75_2020 <- as.numeric(quantile(rh_2020_2_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #rh model_2_ Domain medians & 50% CI 2020
  rh_2_D_25_2020 <- as.numeric(quantile(sample(rh_2020_2_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  rh_2_D_50_2020 <- as.numeric(median(sample(rh_2020_2_INTERACT_domain$Value, 496),na.rm = T))
  rh_2_D_75_2020 <- as.numeric(quantile(sample(rh_2020_2_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #rh model_2_ Domain medians & 50% CI 2100
  rh_2_D_25_2100 <- as.numeric(quantile(sample(rh_2100_2_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  rh_2_D_50_2100 <- as.numeric(median(sample(rh_2100_2_INTERACT_domain$Value, 496),na.rm = T))
  rh_2_D_75_2100 <- as.numeric(quantile(sample(rh_2100_2_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  rh_2[i,] <- round(c(rh_2_ks_2020[1],rh_2_ks_2020[2],rh_2_ks_2020_war[1],rh_2_ks_2020_war[2],
                      rh_2_ks_2100[1],rh_2_ks_2100[2],rh_2_ks_2020_2100[1],rh_2_ks_2020_2100[2],
                      rh_2_I_25_2020,rh_2_I_50_2020,rh_2_I_75_2020,rh_2_W_25_2020,rh_2_W_50_2020,rh_2_W_75_2020,
                      rh_2_D_25_2020,rh_2_D_50_2020,rh_2_D_75_2020,rh_2_D_25_2100,rh_2_D_50_2100,rh_2_D_75_2100),4)
  ##############################################################################################################################
  #rh model_3_ KS test stats 
  rh_3_ks_2020 <- as.numeric(ks.test(rh_2020_3_INTERACT_sites, sample(rh_2020_3_INTERACT_domain$Value, 496))[c(1:2)])
  rh_3_ks_2020_war <- as.numeric(ks.test(rh_2020_3_INTERACT_sites_war, sample(rh_2020_3_INTERACT_domain$Value, 496))[c(1:2)])
  rh_3_ks_2100 <- as.numeric(ks.test(rh_2100_3_INTERACT_sites, sample(rh_2100_3_INTERACT_domain$Value, 496))[c(1:2)])
  rh_3_ks_2020_2100 <- as.numeric(ks.test(sample(rh_2020_3_INTERACT_domain$Value, 496), sample(rh_2100_3_INTERACT_domain$Value, 496))[c(1:2)])
  #rh model_3_ INTERACT medians
  rh_3_I_25_2020 <- as.numeric(quantile(rh_2020_3_INTERACT_sites,probs = c(0.25),na.rm = T))
  rh_3_I_50_2020 <- as.numeric(median(rh_2020_3_INTERACT_sites,na.rm = T))
  rh_3_I_75_2020 <- as.numeric(quantile(rh_2020_3_INTERACT_sites,probs = c(0.75),na.rm = T))
  #rh model_3_ INTERACT (without Russia) medians
  rh_3_W_25_2020 <- as.numeric(quantile(rh_2020_3_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  rh_3_W_50_2020 <- as.numeric(median(rh_2020_3_INTERACT_sites_war,na.rm = T))
  rh_3_W_75_2020 <- as.numeric(quantile(rh_2020_3_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #rh model_3_ Domain medians & 50% CI 2020
  rh_3_D_25_2020 <- as.numeric(quantile(sample(rh_2020_3_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  rh_3_D_50_2020 <- as.numeric(median(sample(rh_2020_3_INTERACT_domain$Value, 496),na.rm = T))
  rh_3_D_75_2020 <- as.numeric(quantile(sample(rh_2020_3_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #rh model_3_ Domain medians & 50% CI 2100
  rh_3_D_25_2100 <- as.numeric(quantile(sample(rh_2100_3_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  rh_3_D_50_2100 <- as.numeric(median(sample(rh_2100_3_INTERACT_domain$Value, 496),na.rm = T))
  rh_3_D_75_2100 <- as.numeric(quantile(sample(rh_2100_3_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  rh_3[i,] <- round(c(rh_3_ks_2020[1],rh_3_ks_2020[2],rh_3_ks_2020_war[1],rh_3_ks_2020_war[2],
                      rh_3_ks_2100[1],rh_3_ks_2100[2],rh_3_ks_2020_2100[1],rh_3_ks_2020_2100[2],
                      rh_3_I_25_2020,rh_3_I_50_2020,rh_3_I_75_2020,rh_3_W_25_2020,rh_3_W_50_2020,rh_3_W_75_2020,
                      rh_3_D_25_2020,rh_3_D_50_2020,rh_3_D_75_2020,rh_3_D_25_2100,rh_3_D_50_2100,rh_3_D_75_2100),4)
  ##############################################################################################################################
  #rh model_4_ KS test stats 
  rh_4_ks_2020 <- as.numeric(ks.test(rh_2020_4_INTERACT_sites, sample(rh_2020_4_INTERACT_domain$Value, 496))[c(1:2)])
  rh_4_ks_2020_war <- as.numeric(ks.test(rh_2020_4_INTERACT_sites_war, sample(rh_2020_4_INTERACT_domain$Value, 496))[c(1:2)])
  rh_4_ks_2100 <- as.numeric(ks.test(rh_2100_4_INTERACT_sites, sample(rh_2100_4_INTERACT_domain$Value, 496))[c(1:2)])
  rh_4_ks_2020_2100 <- as.numeric(ks.test(sample(rh_2020_4_INTERACT_domain$Value, 496), sample(rh_2100_4_INTERACT_domain$Value, 496))[c(1:2)])
  #rh model_4_ INTERACT medians
  rh_4_I_25_2020 <- as.numeric(quantile(rh_2020_4_INTERACT_sites,probs = c(0.25),na.rm = T))
  rh_4_I_50_2020 <- as.numeric(median(rh_2020_4_INTERACT_sites,na.rm = T))
  rh_4_I_75_2020 <- as.numeric(quantile(rh_2020_4_INTERACT_sites,probs = c(0.75),na.rm = T))
  #rh model_4_ INTERACT (without Russia) medians
  rh_4_W_25_2020 <- as.numeric(quantile(rh_2020_4_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  rh_4_W_50_2020 <- as.numeric(median(rh_2020_4_INTERACT_sites_war,na.rm = T))
  rh_4_W_75_2020 <- as.numeric(quantile(rh_2020_4_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #rh model_4_ Domain medians & 50% CI 2020
  rh_4_D_25_2020 <- as.numeric(quantile(sample(rh_2020_4_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  rh_4_D_50_2020 <- as.numeric(median(sample(rh_2020_4_INTERACT_domain$Value, 496),na.rm = T))
  rh_4_D_75_2020 <- as.numeric(quantile(sample(rh_2020_4_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #rh model_4_ Domain medians & 50% CI 2100
  rh_4_D_25_2100 <- as.numeric(quantile(sample(rh_2100_4_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  rh_4_D_50_2100 <- as.numeric(median(sample(rh_2100_4_INTERACT_domain$Value, 496),na.rm = T))
  rh_4_D_75_2100 <- as.numeric(quantile(sample(rh_2100_4_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  rh_4[i,] <- round(c(rh_4_ks_2020[1],rh_4_ks_2020[2],rh_4_ks_2020_war[1],rh_4_ks_2020_war[2],
                      rh_4_ks_2100[1],rh_4_ks_2100[2],rh_4_ks_2020_2100[1],rh_4_ks_2020_2100[2],
                      rh_4_I_25_2020,rh_4_I_50_2020,rh_4_I_75_2020,rh_4_W_25_2020,rh_4_W_50_2020,rh_4_W_75_2020,
                      rh_4_D_25_2020,rh_4_D_50_2020,rh_4_D_75_2020,rh_4_D_25_2100,rh_4_D_50_2100,rh_4_D_75_2100),4)
  ##############################################################################################################################
  #rh model_5_ KS test stats 
  rh_5_ks_2020 <- as.numeric(ks.test(rh_2020_5_INTERACT_sites, sample(rh_2020_5_INTERACT_domain$Value, 496))[c(1:2)])
  rh_5_ks_2020_war <- as.numeric(ks.test(rh_2020_5_INTERACT_sites_war, sample(rh_2020_5_INTERACT_domain$Value, 496))[c(1:2)])
  rh_5_ks_2100 <- as.numeric(ks.test(rh_2100_5_INTERACT_sites, sample(rh_2100_5_INTERACT_domain$Value, 496))[c(1:2)])
  rh_5_ks_2020_2100 <- as.numeric(ks.test(sample(rh_2020_5_INTERACT_domain$Value, 496), sample(rh_2100_5_INTERACT_domain$Value, 496))[c(1:2)])
  #rh model_5_ INTERACT medians
  rh_5_I_25_2020 <- as.numeric(quantile(rh_2020_5_INTERACT_sites,probs = c(0.25),na.rm = T))
  rh_5_I_50_2020 <- as.numeric(median(rh_2020_5_INTERACT_sites,na.rm = T))
  rh_5_I_75_2020 <- as.numeric(quantile(rh_2020_5_INTERACT_sites,probs = c(0.75),na.rm = T))
  #rh model_5_ INTERACT (without Russia) medians
  rh_5_W_25_2020 <- as.numeric(quantile(rh_2020_5_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  rh_5_W_50_2020 <- as.numeric(median(rh_2020_5_INTERACT_sites_war,na.rm = T))
  rh_5_W_75_2020 <- as.numeric(quantile(rh_2020_5_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #rh model_5_ Domain medians & 50% CI 2020
  rh_5_D_25_2020 <- as.numeric(quantile(sample(rh_2020_5_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  rh_5_D_50_2020 <- as.numeric(median(sample(rh_2020_5_INTERACT_domain$Value, 496),na.rm = T))
  rh_5_D_75_2020 <- as.numeric(quantile(sample(rh_2020_5_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #rh model_5_ Domain medians & 50% CI 2100
  rh_5_D_25_2100 <- as.numeric(quantile(sample(rh_2100_5_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  rh_5_D_50_2100 <- as.numeric(median(sample(rh_2100_5_INTERACT_domain$Value, 496),na.rm = T))
  rh_5_D_75_2100 <- as.numeric(quantile(sample(rh_2100_5_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  rh_5[i,] <- round(c(rh_5_ks_2020[1],rh_5_ks_2020[2],rh_5_ks_2020_war[1],rh_5_ks_2020_war[2],
                      rh_5_ks_2100[1],rh_5_ks_2100[2],rh_5_ks_2020_2100[1],rh_5_ks_2020_2100[2],
                      rh_5_I_25_2020,rh_5_I_50_2020,rh_5_I_75_2020,rh_5_W_25_2020,rh_5_W_50_2020,rh_5_W_75_2020,
                      rh_5_D_25_2020,rh_5_D_50_2020,rh_5_D_75_2020,rh_5_D_25_2100,rh_5_D_50_2100,rh_5_D_75_2100),4)
  ##############################################################################################################################
  #rh model_6_ KS test stats 
  rh_6_ks_2020 <- as.numeric(ks.test(rh_2020_6_INTERACT_sites, sample(rh_2020_6_INTERACT_domain$Value, 496))[c(1:2)])
  rh_6_ks_2020_war <- as.numeric(ks.test(rh_2020_6_INTERACT_sites_war, sample(rh_2020_6_INTERACT_domain$Value, 496))[c(1:2)])
  rh_6_ks_2100 <- as.numeric(ks.test(rh_2100_6_INTERACT_sites, sample(rh_2100_6_INTERACT_domain$Value, 496))[c(1:2)])
  rh_6_ks_2020_2100 <- as.numeric(ks.test(sample(rh_2020_6_INTERACT_domain$Value, 496), sample(rh_2100_6_INTERACT_domain$Value, 496))[c(1:2)])
  #rh model_6_ INTERACT medians
  rh_6_I_25_2020 <- as.numeric(quantile(rh_2020_6_INTERACT_sites,probs = c(0.25),na.rm = T))
  rh_6_I_50_2020 <- as.numeric(median(rh_2020_6_INTERACT_sites,na.rm = T))
  rh_6_I_75_2020 <- as.numeric(quantile(rh_2020_6_INTERACT_sites,probs = c(0.75),na.rm = T))
  #rh model_6_ INTERACT (without Russia) medians
  rh_6_W_25_2020 <- as.numeric(quantile(rh_2020_6_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  rh_6_W_50_2020 <- as.numeric(median(rh_2020_6_INTERACT_sites_war,na.rm = T))
  rh_6_W_75_2020 <- as.numeric(quantile(rh_2020_6_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #rh model_6_ Domain medians & 50% CI 2020
  rh_6_D_25_2020 <- as.numeric(quantile(sample(rh_2020_6_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  rh_6_D_50_2020 <- as.numeric(median(sample(rh_2020_6_INTERACT_domain$Value, 496),na.rm = T))
  rh_6_D_75_2020 <- as.numeric(quantile(sample(rh_2020_6_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #rh model_6_ Domain medians & 50% CI 2100
  rh_6_D_25_2100 <- as.numeric(quantile(sample(rh_2100_6_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  rh_6_D_50_2100 <- as.numeric(median(sample(rh_2100_6_INTERACT_domain$Value, 496),na.rm = T))
  rh_6_D_75_2100 <- as.numeric(quantile(sample(rh_2100_6_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  rh_6[i,] <- round(c(rh_6_ks_2020[1],rh_6_ks_2020[2],rh_6_ks_2020_war[1],rh_6_ks_2020_war[2],
                      rh_6_ks_2100[1],rh_6_ks_2100[2],rh_6_ks_2020_2100[1],rh_6_ks_2020_2100[2],
                      rh_6_I_25_2020,rh_6_I_50_2020,rh_6_I_75_2020,rh_6_W_25_2020,rh_6_W_50_2020,rh_6_W_75_2020,
                      rh_6_D_25_2020,rh_6_D_50_2020,rh_6_D_75_2020,rh_6_D_25_2100,rh_6_D_50_2100,rh_6_D_75_2100),4)
  ##############################################################################################################################
  #rh model_7_ KS test stats 
  rh_7_ks_2020 <- as.numeric(ks.test(rh_2020_7_INTERACT_sites, sample(rh_2020_7_INTERACT_domain$Value, 496))[c(1:2)])
  rh_7_ks_2020_war <- as.numeric(ks.test(rh_2020_7_INTERACT_sites_war, sample(rh_2020_7_INTERACT_domain$Value, 496))[c(1:2)])
  rh_7_ks_2100 <- as.numeric(ks.test(rh_2100_7_INTERACT_sites, sample(rh_2100_7_INTERACT_domain$Value, 496))[c(1:2)])
  rh_7_ks_2020_2100 <- as.numeric(ks.test(sample(rh_2020_7_INTERACT_domain$Value, 496), sample(rh_2100_7_INTERACT_domain$Value, 496))[c(1:2)])
  #rh model_7_ INTERACT medians
  rh_7_I_25_2020 <- as.numeric(quantile(rh_2020_7_INTERACT_sites,probs = c(0.25),na.rm = T))
  rh_7_I_50_2020 <- as.numeric(median(rh_2020_7_INTERACT_sites,na.rm = T))
  rh_7_I_75_2020 <- as.numeric(quantile(rh_2020_7_INTERACT_sites,probs = c(0.75),na.rm = T))
  #rh model_7_ INTERACT (without Russia) medians
  rh_7_W_25_2020 <- as.numeric(quantile(rh_2020_7_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  rh_7_W_50_2020 <- as.numeric(median(rh_2020_7_INTERACT_sites_war,na.rm = T))
  rh_7_W_75_2020 <- as.numeric(quantile(rh_2020_7_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #rh model_7_ Domain medians & 50% CI 2020
  rh_7_D_25_2020 <- as.numeric(quantile(sample(rh_2020_7_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  rh_7_D_50_2020 <- as.numeric(median(sample(rh_2020_7_INTERACT_domain$Value, 496),na.rm = T))
  rh_7_D_75_2020 <- as.numeric(quantile(sample(rh_2020_7_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #rh model_7_ Domain medians & 50% CI 2100
  rh_7_D_25_2100 <- as.numeric(quantile(sample(rh_2100_7_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  rh_7_D_50_2100 <- as.numeric(median(sample(rh_2100_7_INTERACT_domain$Value, 496),na.rm = T))
  rh_7_D_75_2100 <- as.numeric(quantile(sample(rh_2100_7_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  rh_7[i,] <- round(c(rh_7_ks_2020[1],rh_7_ks_2020[2],rh_7_ks_2020_war[1],rh_7_ks_2020_war[2],
                      rh_7_ks_2100[1],rh_7_ks_2100[2],rh_7_ks_2020_2100[1],rh_7_ks_2020_2100[2],
                      rh_7_I_25_2020,rh_7_I_50_2020,rh_7_I_75_2020,rh_7_W_25_2020,rh_7_W_50_2020,rh_7_W_75_2020,
                      rh_7_D_25_2020,rh_7_D_50_2020,rh_7_D_75_2020,rh_7_D_25_2100,rh_7_D_50_2100,rh_7_D_75_2100),4)
  
  ##################################################################################################################################
  #rh model_8_ KS test stats 
  rh_8_ks_2020 <- as.numeric(ks.test(rh_2020_8_INTERACT_sites, sample(rh_2020_8_INTERACT_domain$Value, 496))[c(1:2)])
  rh_8_ks_2020_war <- as.numeric(ks.test(rh_2020_8_INTERACT_sites_war, sample(rh_2020_8_INTERACT_domain$Value, 496))[c(1:2)])
  rh_8_ks_2100 <- as.numeric(ks.test(rh_2100_8_INTERACT_sites, sample(rh_2100_8_INTERACT_domain$Value, 496))[c(1:2)])
  rh_8_ks_2020_2100 <- as.numeric(ks.test(sample(rh_2020_8_INTERACT_domain$Value, 496), sample(rh_2100_8_INTERACT_domain$Value, 496))[c(1:2)])
  #rh model_8_ INTERACT medians
  rh_8_I_25_2020 <- as.numeric(quantile(rh_2020_8_INTERACT_sites,probs = c(0.25),na.rm = T))
  rh_8_I_50_2020 <- as.numeric(median(rh_2020_8_INTERACT_sites,na.rm = T))
  rh_8_I_75_2020 <- as.numeric(quantile(rh_2020_8_INTERACT_sites,probs = c(0.75),na.rm = T))
  #rh model_8_ INTERACT (without Russia) medians
  rh_8_W_25_2020 <- as.numeric(quantile(rh_2020_8_INTERACT_sites_war,probs = c(0.25),na.rm = T))
  rh_8_W_50_2020 <- as.numeric(median(rh_2020_8_INTERACT_sites_war,na.rm = T))
  rh_8_W_75_2020 <- as.numeric(quantile(rh_2020_8_INTERACT_sites_war,probs = c(0.75),na.rm = T))
  #rh model_8_ Domain medians & 50% CI 2020
  rh_8_D_25_2020 <- as.numeric(quantile(sample(rh_2020_8_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  rh_8_D_50_2020 <- as.numeric(median(sample(rh_2020_8_INTERACT_domain$Value, 496),na.rm = T))
  rh_8_D_75_2020 <- as.numeric(quantile(sample(rh_2020_8_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #rh model_8_ Domain medians & 50% CI 2100
  rh_8_D_25_2100 <- as.numeric(quantile(sample(rh_2100_8_INTERACT_domain$Value, 496),probs = c(0.25),na.rm = T))
  rh_8_D_50_2100 <- as.numeric(median(sample(rh_2100_8_INTERACT_domain$Value, 496),na.rm = T))
  rh_8_D_75_2100 <- as.numeric(quantile(sample(rh_2100_8_INTERACT_domain$Value, 496),probs = c(0.75),na.rm = T))
  #add all to same df
  rh_8[i,] <- round(c(rh_8_ks_2020[1],rh_8_ks_2020[2],rh_8_ks_2020_war[1],rh_8_ks_2020_war[2],
                      rh_8_ks_2100[1],rh_8_ks_2100[2],rh_8_ks_2020_2100[1],rh_8_ks_2020_2100[2],
                      rh_8_I_25_2020,rh_8_I_50_2020,rh_8_I_75_2020,rh_8_W_25_2020,rh_8_W_50_2020,rh_8_W_75_2020,
                      rh_8_D_25_2020,rh_8_D_50_2020,rh_8_D_75_2020,rh_8_D_25_2100,rh_8_D_50_2100,rh_8_D_75_2100),4)
  
}


tas_1<-as.data.frame(tas_1);names(tas_1)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");tas_1$Variable<-"Tas";tas_1$Model<-"EC-Earth3"
tas_2<-as.data.frame(tas_2);names(tas_2)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");tas_2$Variable<-"Tas";tas_2$Model<-"NorESM2"
tas_3<-as.data.frame(tas_3);names(tas_3)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");tas_3$Variable<-"Tas";tas_3$Model<-"IPSL"
tas_4<-as.data.frame(tas_4);names(tas_4)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");tas_4$Variable<-"Tas";tas_4$Model<-"MPI"
tas_5<-as.data.frame(tas_5);names(tas_5)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");tas_5$Variable<-"Tas";tas_5$Model<-"ACCESS"
tas_6<-as.data.frame(tas_6);names(tas_6)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");tas_6$Variable<-"Tas";tas_6$Model<-"BCC"
tas_7<-as.data.frame(tas_7);names(tas_7)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");tas_7$Variable<-"Tas";tas_7$Model<-"CanESM5"
tas_8<-as.data.frame(tas_8);names(tas_8)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");tas_8$Variable<-"Tas";tas_8$Model<-"CMCC"

pr_1<-as.data.frame(pr_1);names(pr_1)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");pr_1$Variable<-"Pr";pr_1$Model<-"EC-Earth3"
pr_2<-as.data.frame(pr_2);names(pr_2)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");pr_2$Variable<-"Pr";pr_2$Model<-"NorESM2"
pr_3<-as.data.frame(pr_3);names(pr_3)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");pr_3$Variable<-"Pr";pr_3$Model<-"IPSL"
pr_4<-as.data.frame(pr_4);names(pr_4)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");pr_4$Variable<-"Pr";pr_4$Model<-"MPI"
pr_5<-as.data.frame(pr_5);names(pr_5)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");pr_5$Variable<-"Pr";pr_5$Model<-"ACCESS"
pr_6<-as.data.frame(pr_6);names(pr_6)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");pr_6$Variable<-"Pr";pr_6$Model<-"BCC"
pr_7<-as.data.frame(pr_7);names(pr_7)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");pr_7$Variable<-"Pr";pr_7$Model<-"CanESM5"
pr_8<-as.data.frame(pr_8);names(pr_8)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");pr_8$Variable<-"Pr";pr_8$Model<-"CMCC"

snd_1<-as.data.frame(snd_1);names(snd_1)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");snd_1$Variable<-"SnowD";snd_1$Model<-"EC-Earth3"
snd_2<-as.data.frame(snd_2);names(snd_2)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");snd_2$Variable<-"SnowD";snd_2$Model<-"NorESM2"
snd_3<-as.data.frame(snd_3);names(snd_3)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");snd_3$Variable<-"SnowD";snd_3$Model<-"IPSL"
snd_4<-as.data.frame(snd_4);names(snd_4)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");snd_4$Variable<-"SnowD";snd_4$Model<-"MPI"
snd_5<-as.data.frame(snd_5);names(snd_5)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");snd_5$Variable<-"SnowD";snd_5$Model<-"ACCESS"
snd_6<-as.data.frame(snd_6);names(snd_6)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");snd_6$Variable<-"SnowD";snd_6$Model<-"BCC"
snd_7<-as.data.frame(snd_7);names(snd_7)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");snd_7$Variable<-"SnowD";snd_7$Model<-"CanESM5"
snd_8<-as.data.frame(snd_8);names(snd_8)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");snd_8$Variable<-"SnowD";snd_8$Model<-"CMCC"

mrsos_1<-as.data.frame(mrsos_1);names(mrsos_1)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");mrsos_1$Variable<-"SoilM";mrsos_1$Model<-"EC-Earth3"
mrsos_2<-as.data.frame(mrsos_2);names(mrsos_2)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");mrsos_2$Variable<-"SoilM";mrsos_2$Model<-"NorESM2"
mrsos_3<-as.data.frame(mrsos_3);names(mrsos_3)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");mrsos_3$Variable<-"SoilM";mrsos_3$Model<-"IPSL"
mrsos_4<-as.data.frame(mrsos_4);names(mrsos_4)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");mrsos_4$Variable<-"SoilM";mrsos_4$Model<-"MPI"
mrsos_5<-as.data.frame(mrsos_5);names(mrsos_5)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");mrsos_5$Variable<-"SoilM";mrsos_5$Model<-"ACCESS"
mrsos_6<-as.data.frame(mrsos_6);names(mrsos_6)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");mrsos_6$Variable<-"SoilM";mrsos_6$Model<-"BCC"
mrsos_7<-as.data.frame(mrsos_7);names(mrsos_7)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");mrsos_7$Variable<-"SoilM";mrsos_7$Model<-"CanESM5"
mrsos_8<-as.data.frame(mrsos_8);names(mrsos_8)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");mrsos_8$Variable<-"SoilM";mrsos_8$Model<-"CMCC"

lai_1<-as.data.frame(lai_1);names(lai_1)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");lai_1$Variable<-"LAI";lai_1$Model<-"EC-Earth3"
lai_2<-as.data.frame(lai_2);names(lai_2)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");lai_2$Variable<-"LAI";lai_2$Model<-"NorESM2"
lai_3<-as.data.frame(lai_3);names(lai_3)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");lai_3$Variable<-"LAI";lai_3$Model<-"IPSL"
lai_4<-as.data.frame(lai_4);names(lai_4)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");lai_4$Variable<-"LAI";lai_4$Model<-"MPI"
lai_5<-as.data.frame(lai_5);names(lai_5)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");lai_5$Variable<-"LAI";lai_5$Model<-"ACCESS"
lai_6<-as.data.frame(lai_6);names(lai_6)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");lai_6$Variable<-"LAI";lai_6$Model<-"BCC"
lai_7<-as.data.frame(lai_7);names(lai_7)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");lai_7$Variable<-"LAI";lai_7$Model<-"CanESM5"
lai_8<-as.data.frame(lai_8);names(lai_8)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");lai_8$Variable<-"LAI";lai_8$Model<-"CMCC"

cveg_1<-as.data.frame(cveg_1);names(cveg_1)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");cveg_1$Variable<-"Cveg";cveg_1$Model<-"EC-Earth3"
cveg_2<-as.data.frame(cveg_2);names(cveg_2)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");cveg_2$Variable<-"Cveg";cveg_2$Model<-"NorESM2"
cveg_3<-as.data.frame(cveg_3);names(cveg_3)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");cveg_3$Variable<-"Cveg";cveg_3$Model<-"IPSL"
cveg_4<-as.data.frame(cveg_4);names(cveg_4)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");cveg_4$Variable<-"Cveg";cveg_4$Model<-"MPI"
cveg_5<-as.data.frame(cveg_5);names(cveg_5)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");cveg_5$Variable<-"Cveg";cveg_5$Model<-"ACCESS"
cveg_6<-as.data.frame(cveg_6);names(cveg_6)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");cveg_6$Variable<-"Cveg";cveg_6$Model<-"BCC"
cveg_7<-as.data.frame(cveg_7);names(cveg_7)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");cveg_7$Variable<-"Cveg";cveg_7$Model<-"CanESM5"
cveg_8<-as.data.frame(cveg_8);names(cveg_8)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");cveg_8$Variable<-"Cveg";cveg_8$Model<-"CMCC"

csoil_1<-as.data.frame(csoil_1);names(csoil_1)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");csoil_1$Variable<-"Csoil";csoil_1$Model<-"EC-Earth3"
csoil_2<-as.data.frame(csoil_2);names(csoil_2)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");csoil_2$Variable<-"Csoil";csoil_2$Model<-"NorESM2"
csoil_3<-as.data.frame(csoil_3);names(csoil_3)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");csoil_3$Variable<-"Csoil";csoil_3$Model<-"IPSL"
csoil_4<-as.data.frame(csoil_4);names(csoil_4)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");csoil_4$Variable<-"Csoil";csoil_4$Model<-"MPI"
csoil_5<-as.data.frame(csoil_5);names(csoil_5)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");csoil_5$Variable<-"Csoil";csoil_5$Model<-"ACCESS"
csoil_6<-as.data.frame(csoil_6);names(csoil_6)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");csoil_6$Variable<-"Csoil";csoil_6$Model<-"BCC"
csoil_7<-as.data.frame(csoil_7);names(csoil_7)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");csoil_7$Variable<-"Csoil";csoil_7$Model<-"CanESM5"
csoil_8<-as.data.frame(csoil_8);names(csoil_8)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");csoil_8$Variable<-"Csoil";csoil_8$Model<-"CMCC"

npp_1<-as.data.frame(npp_1);names(npp_1)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");npp_1$Variable<-"NPP";npp_1$Model<-"EC-Earth3"
npp_2<-as.data.frame(npp_2);names(npp_2)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");npp_2$Variable<-"NPP";npp_2$Model<-"NorESM2"
npp_3<-as.data.frame(npp_3);names(npp_3)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");npp_3$Variable<-"NPP";npp_3$Model<-"IPSL"
npp_4<-as.data.frame(npp_4);names(npp_4)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");npp_4$Variable<-"NPP";npp_4$Model<-"MPI"
npp_5<-as.data.frame(npp_5);names(npp_5)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");npp_5$Variable<-"NPP";npp_5$Model<-"ACCESS"
npp_6<-as.data.frame(npp_6);names(npp_6)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");npp_6$Variable<-"NPP";npp_6$Model<-"BCC"
npp_7<-as.data.frame(npp_7);names(npp_7)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");npp_7$Variable<-"NPP";npp_7$Model<-"CanESM5"
npp_8<-as.data.frame(npp_8);names(npp_8)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");npp_8$Variable<-"NPP";npp_8$Model<-"CMCC"

rh_1<-as.data.frame(rh_1);names(rh_1)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");rh_1$Variable<-"rH";rh_1$Model<-"EC-Earth3"
rh_2<-as.data.frame(rh_2);names(rh_2)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");rh_2$Variable<-"rH";rh_2$Model<-"NorESM2"
rh_3<-as.data.frame(rh_3);names(rh_3)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");rh_3$Variable<-"rH";rh_3$Model<-"IPSL"
rh_4<-as.data.frame(rh_4);names(rh_4)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");rh_4$Variable<-"rH";rh_4$Model<-"MPI"
rh_5<-as.data.frame(rh_5);names(rh_5)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");rh_5$Variable<-"rH";rh_5$Model<-"ACCESS"
rh_6<-as.data.frame(rh_6);names(rh_6)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");rh_6$Variable<-"rH";rh_6$Model<-"BCC"
rh_7<-as.data.frame(rh_7);names(rh_7)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");rh_7$Variable<-"rH";rh_7$Model<-"CanESM5"
rh_8<-as.data.frame(rh_8);names(rh_8)<-c("D_2020","PV_2020","D_2020_war","PV_2020_war","D_2100","PV_2100","D_2020_2100","PV_2020_2100","CI25_I_2020","CI50_I_2020","CI75_I_2020","CI25_W_2020","CI50_W_2020","CI75_W_2020","CI25_D_2020","CI50_D_2020","CI75_D_2020","CI25_D_2100","CI50_D_2100","CI75_D_2100");rh_8$Variable<-"rH";rh_8$Model<-"CMCC"

CMIP6_Ds_Ms<-rbind(tas_1,tas_2,tas_3,tas_4,tas_5,tas_6,tas_7,tas_8,
                   pr_1,pr_2,pr_3,pr_4,pr_5,pr_6,pr_7,pr_8,
                   snd_1,snd_2,snd_3,snd_6,snd_7,snd_8,
                   mrsos_1,mrsos_2,mrsos_3,mrsos_4,mrsos_5,mrsos_6,mrsos_7,mrsos_8,
                   lai_1,lai_2,lai_3,lai_4,lai_5,lai_6,lai_7,lai_8,
                   cveg_1,cveg_2,cveg_3,cveg_4,cveg_5,cveg_6,cveg_7,cveg_8,
                   csoil_1,csoil_2,csoil_3,csoil_4,csoil_5,csoil_6,csoil_7,csoil_8,
                   npp_1,npp_2,npp_3,npp_4,npp_5,npp_6,npp_7,npp_8,
                   rh_1,rh_2,rh_3,rh_4,rh_5,rh_6,rh_7,rh_8)

#write.csv(CMIP6_Ds_Ms, file = "/Users/elb/OneDrive - Grønlands Naturinstitut/GL2100/CMIP6/ScenarioMIP/PDFs_v2/timeseries/CMIP6_Ds_Ms_V3.csv", row.names = FALSE) 
