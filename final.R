require(bayestestR)
require(parallel)
require(coda)
require(ggplot2)

setwd("~/Desktop/HPCresults")

beta_name <- c("Kf_GLC_Transport                        ",
               "Kr_GLC_Transport                        ",
               "Kf_GLC_SerPool1                         ",
               "Kf_GLC_PYR                              ",
               "Kf_GLC_PRPP                             ",
               "Kf_SerPool3_Synthesis                   ",
               "Kf_SerPool3_Degradation                 ",
               "Kf_GlyPool3_Synthesis                   ",
               "Kf_GlyPool3_Degradation                 ",
               "Kf_SerPool1_GlyPool1                      ",
               "Kr_SerPool1_GlyPool1                      ",
               "Kf_Ser_Transport                        ",
               "Kr_Ser_Transport                        ",
               "Kf_Gly_Transport                        ",
               "Kr_Gly_Transport                        ",
               "Kf_SerPool2_SerPoolMitochon             ",
               "Kr_SerPool2_SerPoolMitochon             ",
               "Kf_GlyPool2_GlyPoolMitochon             ",
               "Kr_GlyPool2_GlyPoolMitochon             ",
               "Kf_SerPoolMitochon_GlyPoolMitochon      ",
               "Kr_SerPoolMitochon_GlyPoolMitochon      ",
               "Kf_GlyPool1_CO2                          ",
               "Kr_GlyPool1_CO2                          ",
               "Kf_MethyleneTHF_FormylTHF_Cytoplasm     ",
               "Kr_MethyleneTHF_FormylTHF_Cytoplasm     ",
               "Kf_FormylTHF_Formate_Cytoplasm          ",
               "Kr_FormylTHF_Formate_Cytoplasm          ",
               "Kf_GlyPoolMitochon_CO2                  ",
               "Kr_GlyPoolMitochon_CO2                  ",
               "Kf_MethyleneTHF_FormylTHF_Mitochon      ",
               "Kr_MethyleneTHF_FormylTHF_Mitochon      ",
               "Kf_FormylTHF_Formate_Mitochon           ",
               "Kr_FormylTHF_Formate_Mitochon           ",
               "Kf_PRPP1_GAR                            ",
               "Kf_GAR_FGAR                             ",
               "Kf_FGAR_AMP                             ",
               "Kf_IMP_AMP                              ",
               "Kf_AMP_Degradation                      ",
               "Kr_MethyleneTHF_Transport               ",
               "Kf_MethyleneTHF_Transport               ",
               "Kr_FormylTHF_Transport                  ",
               "Kf_FormylTHF_Transport                  ",
               "Kr_THF_Transport                        ",
               "Kf_THF_Transport                        ",
               "Kr_Formate_Transport                    ",
               "Kf_Formate_Transport                    ",
               "Kf_PRPP2_GAR                            ",
               "Kf_PRPP3_GAR                            ",
               "Kf_SerPool2_GlyPool2                    ",
               "Kr_SerPool2_GlyPool2                    ",
               "Kf_GlyPool2_CO2                         ",
               "Kr_GlyPool2_CO2                         ")
beta_formatted <- sapply(1:52,function(i) gsub(" ","",beta_name[i]))

optipar <- read.table("OptimizePar.txt")

initial_beta <- log(c(0.000219,
                      0.0864,
                      0.00044,
                      0.1833,
                      0.0001245,
                      7.68e-06,
                      0.0177,
                      1.92e-05,
                      0.00053,
                      264,
                      2.6625,
                      0.0002085,
                      0.010005,
                      9.045e-05,
                      0.001041,
                      0.0039,
                      0.00975,
                      0.00585,
                      0.014625,
                      831,
                      16.5,
                      53.1,
                      0.639,
                      0.32,
                      0.8,
                      0.32,
                      1000000,
                      434,
                      0.42,
                      2.67,
                      0.267,
                      2.67,
                      44400,
                      88,
                      44000,
                      44000,
                      0.0125,
                      0.000625,
                      0.039,
                      0.0078,
                      4.68,
                      0.468,
                      0.78,
                      0.078,
                      4.68,
                      0.468,
                      44,
                      22,
                      264,
                      2.6625,
                      53.1,
                      0.639))
priorsd_beta <- c(log(0.01)-log(0.005),0-log(0.5),log(0.1)-log(0.05),0-log(0.5),log(0.01)-log(0.005),
                  log(5e-4)-log(2.5e-4),log(1)-log(0.5),log(2e-4)-log(1e-4),log(0.02)-log(0.01),
                  log(500)-log(225),log(10)-log(5),log(0.01)-log(0.005),log(0.1)-log(0.05),log(0.001)-log(0.0005),
                  log(0.02)-log(0.01),log(0.02)-log(0.01),log(0.03)-log(0.015),log(0.01)-log(0.005),
                  log(0.05)-log(0.025),log(1000)-log(450),log(20)-log(10),log(100)-log(45),log(2)-log(1),
                  log(3)-log(1.5),log(10)-log(5),log(3)-log(1.5),log(1e7)-log(5e6),log(1e3)-log(500),
                  log(4)-log(2),log(10)-log(5),0-log(0.5),log(10)-log(5),log(1e5)-log(5e4),log(100)-log(50),
                  log(2e5)-log(0.9e5),log(2e5)-log(0.9e5),log(0.25)-log(0.15),log(0.01)-log(0.005),
                  log(0.5)-log(0.25),log(0.3)-log(0.15),log(50)-log(25),log(5)-log(2.5),log(8)-log(4),
                  log(1)-log(0.5),log(10)-log(5),log(10)-log(5),log(100)-log(50),log(50)-log(25),
                  log(500)-log(225),log(500)-log(250),log(100)-log(50),log(3)-log(1.5))*1.5

far_initial_index <- which(abs(initial_beta-log(optipar$V1[1:52]))/priorsd_beta>1)
initial_distance <- abs(initial_beta-log(optipar$V1[1:52]))/priorsd_beta

beta_lb <- c(0,0.1,0,0.1,0,0,0,0,0,50,
             1,0,0,0,0,0,0,0,0,100,
             1,10,0,0,0.1,0,0.2e6,100,0.1,0.2,
             0.02,0.2,1e4,10,2e4,2e4,0,0,0,0,
             1,0.05,0.1,0.01,0.5,0.1,5,1,50,1,
             1,0) #lower bounds
beta_ub <- c(0.01,1,0.1,1,0.01,5e-4,1,2e-4,2e-2,500,
             10,0.01,0.1,0.001,0.02,0.02,0.03,0.01,0.05,1000,
             20,100,2,3,10,3,1e7,1000,4,10,
             1,10,1e5,100,2e5,2e5,0.25,0.01,0.5,0.3,
             50,5,8,1,10,10,100,50,500,500,
             100,3) #upper bounds
beta_true_nonc <- c(0.000465165222015355, 0.242665045280639, 2.676613e-05, 0.518326766571717, 0.0002, 1.957725e-05, 4.511945e-02,0.000075, 0.002521119, 130.0,
                    0.6, 0.00018 , 0.075, 0.00005, 0.002653635 , 0.001,0.0008 , 0.00974237918290052,   0.0278492039049790 , 780.0 ,
                    0.6, 54.2492566486518, 0.03 , 0.95 ,   2.00461377104628, 0.398825587055585, 999999.999999978,  434.005733100667,  1.05698640506409,2.24239831194499 ,
                    0.551424353965021, 2.31291778438686, 44400.0000092611,  80.0, 44000.0001051987,  44000.0000615634, 0.0268575860192817,0.000861408841917836,   0.107341967375275, 0.00756915225856982,
                    11.6618608584811,  0.499324532242705, 2.01627621832653 ,  0.0682897697038028,4.57835999241956, 0.579450741578679, 42.0, 21.0, 120.0, 1.2,
                    65.0,  0.72)

beta_index <- c( 3,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
                 20, 21, 22, 23, 28, 29, 34, 49, 50, 51,52)

######################################################################################
setwd("~/Desktop/HPCresults/Mac/holly/dim10/Thesis/final/est_full2/")
# full est random error var

mytrue_var <- log(c(1.034238e-04, 2.170496e-01, 1.886177e-05, 2.602462e-01, 1.863844e-03, 5.148025e-06,2.519172e-02, 3.582511e-05, 6.400460e-04, 8.990407e+01,
                    8.611367e+00, 4.696572e-05,2.616996e-02, 2.759432e-05, 7.075162e-03, 2.511795e-03, 9.039947e-03, 3.221557e-03,1.537229e-02, 3.222696e+02,
                    8.592349e+00, 3.588606e+01, 1.946228e-02, 1.141591e+00,6.272787e-01, 2.958951e-01, 4.659520e+05, 2.454078e+02, 2.430246e-01, 2.552116e+00,
                    2.044861e-01, 2.167994e+00, 2.213579e+04 ,4.008122e+01, 4.075222e+04, 4.303686e+04,1.282449e-02, 6.241349e-04, 3.921026e-02, 6.875464e-03,
                    3.184057e+00, 3.513642e-01,8.493559e-01, 4.747550e-02, 2.876920e+00, 2.857837e-01, 2.376757e+01, 1.531322e+01,2.374717e+02, 2.254443e+00,
                    7.039808e+01, 1.270505e-02))
mytrue_err <- c(0.02074755, 0.05413806, 0.06939170, 0.06251849, 0.03294391, 0.13394740, 0.03289497, 0.09676033, 0.04452552,0.02447726,
                0.07792948, 0.03227189 ,0.06168760, 0.01979036, 0.03940984, 0.22289147, 0.10456421, 0.03607298,0.02864494, 0.06800813,
                0.03760034, 0.09821494, 0.07877421, 0.01880766, 0.03435694, 0.07666611, 0.04778826,0.13032726, 0.08356618, 0.01840211,
                0.07138771, 0.03556571, 0.03126719, 0.07504247, 0.03244967, 0.05636749,0.02147176, 0.02265899, 0.02782587, 0.02698518,
                0.05686491, 0.05463770, 0.03416760, 0.02445783, 0.02205466,0.04129188 ,0.05514810, 0.08329250, 0.06319228, 0.03519235,
                0.01885804, 0.20336788, 0.16840891, 0.01908970, 0.02491993, 0.01412577, 0.03858249, 0.14570720, 0.03384012, 0.11420980,
                0.17640609, 0.09572480, 0.05760146, 0.03794238, 0.02795174, 0.23228771, 0.05406738, 0.09107526 ,0.04821502, 0.24135744,
                0.17391265, 0.03032879, 0.05459746, 0.07919798, 0.04053577, 0.03486330, 0.03000994, 0.15013457, 0.09402209, 0.02377423,
                0.03409032, 0.07075026, 0.08589932, 0.09530419, 0.05244093, 0.01114397, 0.10159183, 0.04403797, 0.10084598, 0.05777371,
                0.03569479, 0.05652702, 0.03582980, 0.03513353, 0.09804653)

s <- seq(1,30,1)
# len <- c(16,16,17,16,17,16,15,17,17,15,
#          17,16,17,17,16,17,18,17,17,16,
#          17,18,17,18,16,18,18,16,17,16)
len <- rep(15,30)
data_list <- list(NA)
for(i in 1:length(s)){
  file_names <- lapply(1:len[i], function(x) paste("mysample_",s[i],"curr_",x,".txt", sep = ""))
  your_data_frame <- do.call(rbind,lapply(file_names,read.table,header = T))
  data_list[[i]] <- your_data_frame
}

png("est5_full_cut.png",height = 1050,width = 750)
par(mfrow= c(13,4),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
for(i in 1:52){
  plot(log(data_list[[5]][-(1:5000),i]),type = "l", ylab = i)
  lines(x = c(0,27000),y = c(mytrue_var[i],mytrue_var[i]),col='red')
}
dev.off()

png("err_full.png",height = 2750,width = 1350)
par(mfrow= c(19,5),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
for(i in 53:147){
  plot(data_list[[5]][,i],type = "l", ylab = i-52)
  lines(x = c(0,27000),y = c(mytrue_err[i-52],mytrue_err[i-52]),col='red',lwd=1.5)
}
dev.off()

#credible interval
cov <- list(NA)
for(j in 1:52){
  cov[[j]] <- rep(NA,60)
  for(i in 1:length(s)){
    sample_CI <- hdi(log(data_list[[i]][-(1:5000),j]), ci = 0.95)
    cov[[j]][i] <- sample_CI$CI_low 
    cov[[j]][i+length(s)] <- sample_CI$CI_high
  }
}

# cov <- list(NA)
# for(j in 1:52){
#   cov[[j]] <- rep(NA,60)
#   for(i in 1:length(s)){
#     sample_CI <- quantile(log(data_list[[i]][-(1:4000),j]), c(0.25,0.75))
#     cov[[j]][i] <- sample_CI[1] 
#     cov[[j]][i+length(s)] <- sample_CI[2]
#   }
# }

the52cov <- sapply(1:52, function(j)
  mean(sapply(1:length(s), function(x) cov[[j]][x]<=mytrue_var[j]&&
                cov[[j]][x+30]>=mytrue_var[j] )) )
round(the52cov,3)
mean(the52cov)


theSD<-NA
thebias<-NA
for(dim in 1:52){
  sum <- NA
  sum1 <- NA
  for(i in 1:length(s)){
    # sum <- c(sum, sqrt(mean((log(data_list[[i]][-(1:5000),dim])-mytrue_var[dim])^2)
    # ) )
    sum <- c(sum, sqrt(var(log(data_list[[i]][-(1:5000),dim]))))
    sum1 <- c(sum1, mean(log(data_list[[i]][-(1:5000),dim]))-mytrue_var[dim] )
  }
  theSD <- c(theSD,mean(sum[-1]))
  thebias <- c(thebias,abs(mean(sum1[-1])))
}

round(theSD[-1],3)
round(thebias[-1],3)
round(sqrt(themse[-1]^2-thebias[-1]^2),3)
range(theSD[-1])
range(thebias[-1])

sumdata <- matrix(c(theSD[-1],thebias[-1]),ncol = 2,byrow = F)
pdf("est_full_mse_bias.pdf",height = 2.5,width = 4)
par(mfrow= c(1,1),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
boxplot(sumdata, ylim = range(sumdata), main="",names=c("pSD","ABias"),
     cex.lab = 0.5,cex.axis = 0.5,col = c("red","orange"),border = "brown",
     horizontal = T,boxwex = 0.5)
dev.off()



thebias_rel<-NA
for(dim in 1:52){
  # sum <- NA
  sum1 <- NA
  for(i in 1:length(s)){
    # sum <- c(sum, sqrt(mean((log(data_list[[i]][-(1:5000),dim])-mytrue_var[dim])^2)
    # ) )
    # sum <- c(sum, sqrt(var(log(data_list[[i]][-(1:5000),dim]))))
    sum1 <- c(sum1, (mean(log(data_list[[i]][-(1:5000),dim]))-mytrue_var[dim])/sd(log(data_list[[i]][-(1:5000),dim])) )
  }
  # theSD <- c(theSD,mean(sum[-1]))
  thebias_rel <- c(thebias_rel,mean(sum1[-1]))
}

round(thebias_rel[-1],3)
range(thebias_rel[-1])

sumdata <- matrix(c(thebias_rel[-1]),ncol = 1,byrow = F)
pdf("est_full_mse_bias_relative.pdf",height = 2.5,width = 4)
par(mfrow= c(1,1),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
boxplot(sumdata, ylim = c(-1,1.2), main="",names=c("RBias"),
        cex.lab = 0.5,cex.axis = 0.5,col = c("red"),border = "brown",
        horizontal = T,boxwex = 0.5)
dev.off()



pdf("est_full_bias_all_beta.pdf",height = 26,width = 12)
par(mfrow= c(13,4),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
for(i in 1:52){
  the_data <- sapply(1:30, function(x) mean(log(data_list[[x]][-(1:5000),i])))
  plot(the_data,type = "p", xlab="Simulations", ylab = paste("Post. Mean of beta_",i),ylim = range(the_data,mytrue_var[i]))
  lines(x = c(0,27000),y = c(mytrue_var[i],mytrue_var[i]),col='red')
}
dev.off()

the_data <- matrix(NA,52,30)
for(i in 1:52){
  the_data[i,] <- sapply(1:30, function(x) mean(log(data_list[[x]][-(1:5000),i])))
}
write.table(the_data,
            'postmean_30sims.csv',
            sep=',',row.names = FALSE,col.names = FALSE)


# d <- data.frame(x =seq(1,52,1),
#                cov = the52cov,
#                themse = themse[-1],
#                thebias = thebias[-1])
#                # logp = signif(-log10(runif(10)), 2))
# pdf("est_full_three_plot.pdf",height = 6,width = 8)
# par(mar = c(5,5,2,5))
# with(d, plot(x, themse, type="l", col="red", 
#              xlab = "index",
#              ylab="RMSE or ABias",
#              ylim=c(0,1.4)))
# with(d, lines(x, thebias, col="green",
#              ylim=c(0,1.4)))
# grid()
# par(new = T)
# with(d, plot(x, abs(cov-0.5), pch=16, axes=F, xlab=NA, ylab=NA, cex=0.8))
# axis(side = 4)
# mtext(side = 4, line = 3, '|Coverage Probability-0.5|')
# with(d, lines(x=x,y=abs(cov-0.5),col="black",lty = 1))
# # with(d, lines(x=c(0,52),y=c(0.5,0.5),col="gray",lty = 2))
# legend("topleft",
#        cex = 0.8,
#        legend=c("RMSE","ABias", "|Coverage Probability-0.5|"),
#        lty=c(1,1,0), pch=c(NA,NA,16), col=c("red","green", "black"))
# dev.off()

stability <- sapply(1:52, function(j) geweke.diag(log(data_list[[5]][-(1:5000),j]))$z )
plot(stability)
lines(x=c(0,55),y=c(-1.96,-1.96),col='red')
lines(x=c(0,55),y=c(1.96,1.96),col='red')

stable_check <- rep(NA,30)
for(i in 1:30){
  stability <- sapply(1:52, function(j) geweke.diag(log(data_list[[i]][-(1:5000),j]))$z )
  stable_check[i] <- sum((stability>qnorm(0.025))*(stability<qnorm(0.975)))/52
}

# > stable_check
# [1] 0.6153846 0.7692308 0.6538462 0.7692308 0.5192308 0.6730769 0.6923077 0.8076923 0.7500000
# [10] 0.8269231 0.8076923 0.6923077 0.7692308 0.8076923 0.7307692 0.7115385 0.7307692 0.7115385
# [19] 0.6923077 0.8846154 0.6153846 0.5192308 0.8076923 0.6730769 0.4615385 0.5000000 0.6346154
# [28] 0.7115385 0.7692308 0.6538462

######################################################################################


######################################################################################
setwd("~/Desktop/HPCresults/Mac/holly/dim10/Thesis/final/real_data/cancer/")

file_names <- lapply(1:37, function(x) paste("mysample_curr_",x,".txt", sep = "")) #14
your_data_frame <- do.call(rbind,lapply(file_names,read.table,header = T))

# tune_index <- seq(4300,length(your_data_frame[,1]),5)
png("real_est_full.png",height = 1050,width = 750)
par(mfrow= c(13,4),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
for(i in 1:52){
  plot(log(your_data_frame[,i]),type = "l", ylab = i)
  # lines(x=c(0,20000),y=c(mean(log(your_data_frame[(2000:4000),i])),
  #                        mean(log(your_data_frame[(2000:4000),i]))),
  #       col="red")
}
dev.off()

png("real_est_obs.png",height = 1050,width = 550)
par(mfrow= c(13,2),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
for(i in 1:length(beta_index)){
  plot(log(your_data_frame[-(1:4000),beta_index[i]]),type = "l", ylab = beta_index[i])
  lines(x=c(0,20000),y=c(mean(log(your_data_frame[-(1:4000),beta_index[i]])),
                         mean(log(your_data_frame[-(1:4000),beta_index[i]]))),
        col="red")
}
dev.off()

# the_distance <-(sapply(1:52, function(i) mean(log(your_data_frame[-(1:3000),i])))
#                 - initial_beta)/priorsd_beta
# 
# sapply(1:length(beta_index), function(i) 
#   mean(log(your_data_frame[-(1:3000),beta_index[i]])))

#real data 24h y concentrations
therealy24_c <- matrix(c(log(16.72720275),log(0.008553777),log(0.008553777),log(0.008553777),log(133.3070782),log(0.008553777),log(0.316733723),log(2895.962684),log(8352.975818),log(0.008553777),
                       log(14.613105),log(0.008553777),log(0.008553777),log(0.008553777),log(86.80758989),log(0.008553777),log(0.316733723),log(2330.656011),log(7553.664549),log(0.008553777),
                       log(15.75248504),log(0.008553777),log(0.008553777),log(0.030020138),log(140.9276324),log(0.008553777),log(0.316733723),log(3425.228351),log(12198.03665),log(0.008553777)),nrow=3,byrow = T)

therealy24_cp <- sapply(1:10, function(x) mean(therealy24_c[,x]))
therealy24_cu <- sapply(1:10, function(x) mean(therealy24_c[,x])-sd(therealy24_c[,x]))
therealy24_cd <- sapply(1:10, function(x) mean(therealy24_c[,x])+sd(therealy24_c[,x]))

#lcc/s5/, data distance

###12/08/2022
write.table(your_data_frame[-(1:5000),1:52],
  '~/Desktop/HPCresults/Mac/holly/dim10/Thesis/final/real_data/cancer/results08dec2022/cppinput.csv',
            sep=',',row.names = FALSE,col.names = FALSE)
file_names <- list("ode_est_data10_2022_2501_5000.txt","ode_est_data10_2022_5001_8000.txt") #2
theesty24_ind <- do.call(rbind,lapply(file_names,read.table,header = F))
#theesty24_ind <- read.table("ode_est_data10_2022.txt", header = F) #one row for one sample result.


write.table(your_data_frame[7501:13000,1:52],
            '~/Desktop/HPCresults/Mac/holly/dim10/Thesis/final/real_data/cancer/beta52_mcmc_iterations5500.csv',
            sep=',',row.names = FALSE,col.names = FALSE)
thepostmean <- sapply(1:52, function(i) exp(mean(log(your_data_frame[7501:13000,i]))))
# thepost_p_sd <- sapply(1:52, function(i) exp(mean(log(your_data_frame[7501:13000,i]))+
#                                                sd(log(your_data_frame[7501:13000,i]))))
# thepost_m_sd <- sapply(1:52, function(i) exp(mean(log(your_data_frame[7501:13000,i]))-
#                                                sd(log(your_data_frame[7501:13000,i]))))
# theesty24_c <- read.table("ode_est_data10.txt", header = F)
# theesty24_cu <- read.table("ode_est_data10_u.txt", header = F)
# theesty24_cd <- read.table("ode_est_data10_d.txt", header = F)
#all concentrations at 24h for sd calculation
your_flux24 <- theesty24_ind

theesty24_c <- sapply(1:10, function(i) mean(your_flux24[,i]))
theesty24_cu <- sapply(1:10, function(x) mean(your_flux24[,x])-sd(your_flux24[,x]))
theesty24_cd <- sapply(1:10, function(x) mean(your_flux24[,x])+sd(your_flux24[,x]))

range(therealy24_cp,therealy24_cu,therealy24_cd,theesty24_c,theesty24_cu, theesty24_cd)

# pdf("real_est_comp.pdf",height = 5,width = 5)
# par(mfrow= c(1,1),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
# plot(x=therealy24_cp, y= theesty24_c, ylab = "Estimated log concentration",
#      xlab = "Real data mean (+-sd) log concentration",main = "",col="black",bg="black",
#      xlim=c(-4.92,9.41),ylim=c(-4.9,9.5),cex=0.7,pch=16)
# for(i in 1:10){
#   lines(x=c(therealy24_cd[i],therealy24_cu[i]), y= c(theesty24_c[i],theesty24_c[i]),
#         col="red",bg="red",pch=24,cex=0.4)
# }
# lines(x=c(-10,32),y=c(-10,32),col="gray",lty="dashed")
# dev.off()

(sum(abs(therealy24_c[1,]-theesty24_c))+
  sum(abs(therealy24_c[2,]-theesty24_c))+
  sum(abs(therealy24_c[3,]-theesty24_c)))/3 
 #5.573834



write.table(your_data_frame[7501:13000,1:52],
            '~/Desktop/HPCresults/Mac/holly/dim10/Thesis/final/real_data/cancer/results/flux_band/cppinput2.csv',
            sep=',',row.names = FALSE,col.names = FALSE)

the6flux <- read.table("est_flux.txt",header = F)

v5 <- thepostmean[3] * the6flux[,2]/6.02214129e14
v6 <- thepostmean[3] * the6flux[,3]/6.02214129e14

vv5 <- v5/sd(v5)
vv6 <- v6/sd(v6)

range(vv5,vv6)

pdf("the6flux_est.pdf",height = 5,width = 5)
par(mfrow= c(1,1),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
plot(x=the6flux[,1], y=vv5,type="l", ylab = "Estimated flux standardized by SD",
     xlab = "Time (min)",main = "",col=3,cex=0.4,ylim=c(0,27),xlim=c(0,15),
     lwd=2,lty=1)
lines(x=the6flux[,1],y=vv6,col=6,cex=0.4,lty=2,lwd=2)
legend(5, 10,
       legend = c("Vf_Glc_SerPool1_13C0_13C0D000", "Vf_Glc_SerPool1_13C1_13C3D000"), 
       col = c(3,6), 
       lty = c(1,2),
       lwd = c(2,2),
       cex = 0.8, 
       text.col = "black", 
       box.lty = 1,box.col = "black",
       horiz = F , 
       inset = c(1, 1))
dev.off()


range(v5,v6)

pdf("the6flux_est_no_sd.pdf",height = 5,width = 5)
par(mfrow= c(1,1),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
plot(x=the6flux[,1], y=v5,type="l", ylab = "Estimated flux standardized by SD",
     xlab = "Time (min)",main = "",col=3,cex=0.4,ylim=c(0,0.0035),xlim=c(0,15),
     lwd=2,lty=1)
lines(x=the6flux[,1],y=v6,col=6,cex=0.4,lty=2,lwd=2)
legend(5, 0.003,
       legend = c("Vf_Glc_SerPool1_13C0_13C0D000", "Vf_Glc_SerPool1_13C1_13C3D000"), 
       col = c(3,6), 
       lty = c(1,2),
       lwd = c(2,2),
       cex = 0.8, 
       text.col = "black", 
       box.lty = 1,box.col = "black",
       horiz = F , 
       inset = c(1, 1))
dev.off()

# pdf("the6flux_est.pdf",height = 5,width = 5)
# par(mfrow= c(1,1),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
# plot(x=the6flux[,1], y=vv5,type="l", ylab = "Estimated flux standardized by SD",
#      xlab = "Time (min)",main = "",col=3,cex=0.4,ylim=c(0,50),xlim=c(0,400),
#      lwd=2,lty=1)
# lines(x=the6flux[,1],y=vv6,col=6,cex=0.4,lty=2,lwd=2)
# legend(130, 5,
#        legend = c("Vf_Glc_SerPool1_13C0_13C0D000", "Vf_Glc_SerPool1_13C1_13C3D000"), 
#        col = c(3,6), 
#        lty = c(1,2),
#        lwd = c(2,2),
#        cex = 0.8, 
#        text.col = "black", 
#        box.lty = 1,box.col = "black",
#        horiz = F , 
#        inset = c(1, 1))
# dev.off()

stability_y_conc <- sapply(1:10, function(j) geweke.diag((theesty24_ind[,j]),0.15,0.5)$z )
pnorm(stability_y_conc)
stable_check_y <- sum((stability_y_conc>qnorm(0.025/10))*(stability_y_conc<qnorm(1-0.025/10)))/10


stability_c <- sapply(1:52, function(j) geweke.diag(log(your_data_frame[7501:13000,j]),0.15,0.5)$z )
pval <- (1-pnorm(abs(stability_c)))*2
adj_pval <- p.adjust(pval, method = "bonferroni")
sum(adj_pval>0.05)
stable_check_c <- sum((stability_c>qnorm(0.025/52))*(stability_c<qnorm(1-0.025/52)))/52

stack_vec <- NA
flag0 <- seq(15000,18500,1000)
flag1 <- seq(4000, 14000,1000)

flag0 <- seq(15500,18500,500)
flag1 <- seq(4000, 15000,500)
flag2 <- seq(0.1,0.45,0.05)
flag3 <- seq(0.15,0.5,0.05)
for (i1 in 1:length(flag0)){
  for(i2 in 1:length(flag1)){

        stability <- sapply(1:52, function(j)
          heidel.diag(log(your_data_frame[(flag1[i2]+1):flag0[i1],j]),0.1,0.05)[1,1]
        )
          # geweke.diag(log(your_data_frame[(flag1[i2]+1):flag0[i1],j]),flag2[i3],flag3[i4])$z )
        if(sum(stability)>50){
        stack_vec <- c(stack_vec,i1,i2,i3,i4,sum(stability)/52)
        print(c(i1,i2,sum(stability)))
        }
  }
}

# 13001-15500,0.15,0.35, 50/52

png("real_est_full2_1129.png",height = 1050,width = 750)
par(mfrow= c(13,4),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
for(i in 1:52){
  plot(log(your_data_frame[10001:15000,i]),type = "l", ylab = i)
  # lines(x=c(0,20000),y=c(mean(log(your_data_frame[(2000:4000),i])),
  #                        mean(log(your_data_frame[(2000:4000),i]))),
  #       col="red")
}
dev.off()

stability <- sapply(1:52, function(j) heidel.diag(log(your_data_frame[13001:15500,j]),0.1,0.05)[1,1] )
stable_check <- sum(stability)/52 


#sensitivity
setwd("~/Desktop/HPCresults/Mac/holly/dim10/Thesis/final/real_data/cancer/larger_var2")

file_names <- lapply(1:14, function(x) paste("mysample_curr_",x,".txt", sep = ""))
your_data_frame <- do.call(rbind,lapply(file_names,read.table,header = T))

# tune_index <- seq(4300,length(your_data_frame[,1]),5)
png("real_est_full.png",height = 1050,width = 750)
par(mfrow= c(13,4),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
for(i in 1:52){
  plot(log(your_data_frame[,i]),type = "l", ylab = i)
  # lines(x=c(0,20000),y=c(mean(log(your_data_frame[(2000:4000),i])),
  #                        mean(log(your_data_frame[(2000:4000),i]))),
  #       col="red")
}
dev.off()

png("real_est_obs.png",height = 1050,width = 550)
par(mfrow= c(13,2),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
for(i in 1:length(beta_index)){
  plot(log(your_data_frame[-(1:4000),beta_index[i]]),type = "l", ylab = beta_index[i])
  lines(x=c(0,20000),y=c(mean(log(your_data_frame[-(1:4000),beta_index[i]])),
                         mean(log(your_data_frame[-(1:4000),beta_index[i]]))),
        col="red")
}
dev.off()

#####################################################################################



######################################################################################
setwd("~/Desktop/HPCresults/Mac/holly/dim10/Thesis/final/real_data/lcc/noncancer/") 
#real data est noncancer

file_names <- lapply(1:67, function(x) paste("mysample_curr_",x,".txt", sep = ""))
your_data_frame <- do.call(rbind,lapply(file_names,read.table,header = T))

write.table(your_data_frame[-(1:5000),1:52],
            '~/Desktop/HPCresults/Mac/holly/dim10/Thesis/final/real_data/cancer/results08dec2022/cppinput_nc.csv',
            sep=',',row.names = FALSE,col.names = FALSE)

# # tune_index <- seq(4300,length(your_data_frame[,1]),5)
# png("real_est_full_nc.png",height = 1050,width = 750)
# par(mfrow= c(13,4),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
# for(i in 1:52){
#   plot(log(your_data_frame[,i]),type = "l", ylab = i)
#   # lines(x=c(0,20000),y=c(mean(log(your_data_frame[(2000:4000),i])),
#   #                        mean(log(your_data_frame[(2000:4000),i]))),
#   #       col="red")
# }
# dev.off()
# 
# png("real_est_obs_nc.png",height = 1050,width = 550)
# par(mfrow= c(13,2),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
# for(i in 1:length(beta_index)){
#   plot(log(your_data_frame[-(1:8000),beta_index[i]]),type = "l", ylab = beta_index[i])
#   # lines(x=c(0,20000),y=c(mean(log(your_data_frame[-(1:6000),beta_index[i]])),
#   #                        mean(log(your_data_frame[-(1:6000),beta_index[i]]))),
#   #       col="red")
# }
# dev.off()

# the_distance <-(sapply(1:52, function(i) mean(log(your_data_frame[-(1:8000),i])))
#                 - initial_beta)/priorsd_beta
# 
# sapply(1:length(beta_index), function(i) 
#   mean(log(your_data_frame[-(1:8000),beta_index[i]])))

#real data 24h y concentrations
therealy24 <- matrix(c(log(10.49472718),log(0.008553777),log(0.008553777),log(0.008553777),log(66.09298315),log(0.008553777),log(0.008553777),log(3548.833023),log(7893.239601),log(0.008553777),
                       log(15.46839942),log(0.008553777),log(0.008553777),log(0.008553777),log(82.72434009),log(0.008553777),log(0.008553777),log(2443.951195),log(8555.906373),log(0.008553777),
                       log(7.970413169),log(0.008553777),log(0.008553777),log(0.008553777),log(25.62146494),log(0.008553777),log(0.008553777),log(2992.618267),log(6553.527511),log(0.008553777)),nrow=3,byrow = T)

therealy24_p <- sapply(1:10, function(x) mean(therealy24[,x]))
therealy24_u <- sapply(1:10, function(x) mean(therealy24[,x])-sd(therealy24[,x]))
therealy24_d <- sapply(1:10, function(x) mean(therealy24[,x])+sd(therealy24[,x]))

#lcc/s5/, data distance
thepostmean <- sapply(1:52, function(i) exp(mean(log(your_data_frame[7501:13000,i]))))
# thepost_p_sd <- sapply(1:52, function(i) exp(mean(log(your_data_frame[7501:13000,i]))+
#                                                sd(log(your_data_frame[7501:13000,i]))))
# thepost_m_sd <- sapply(1:52, function(i) exp(mean(log(your_data_frame[7501:13000,i]))-
#                                                sd(log(your_data_frame[7501:13000,i]))))
# theesty24 <- read.table("ode_est_data10.txt", header = F)
# theesty24_u <- read.table("ode_est_data10_u.txt", header = F)
# theesty24_d <- read.table("ode_est_data10_d.txt", header = F)

#all concentrations at 24h for sd calculation
file_names <- list("ode_est_data10_2022_nc_2501_5000.txt","ode_est_data10_2022_nc_5001_8000.txt")
your_flux24 <- do.call(rbind,lapply(file_names,read.table,header = F))

theesty24 <- sapply(1:10, function(i) mean(your_flux24[,i]))
theesty24_u <- sapply(1:10, function(x) mean(your_flux24[,x])-sd(your_flux24[,x]))
theesty24_d <- sapply(1:10, function(x) mean(your_flux24[,x])+sd(your_flux24[,x]))

# pdf("real_est_comp_n.pdf",height = 5,width = 5)
# par(mfrow= c(1,1),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
# plot(x=therealy24[1,], y= theesty24, ylab = "Estimated log concentration",
#      xlab = "Real data log concentration",main = "",col="green",bg="green",
#      xlim=c(-4.77,9.06),ylim=c(-8.66,8.43),cex=0.4,pch=24)
# points(x=therealy24[2,], y= theesty24,col="blue",bg="blue",pch=24,cex=0.4)
# points(x=therealy24[3,], y= theesty24,col="red",bg="red",pch=24,cex=0.4)
# lines(x=c(-10,32),y=c(-10,32),col="gray",lty="dashed")
# dev.off()

# sum(abs(therealy24[1,]-theesty24))+
#   sum(abs(therealy24[2,]-theesty24))+
#   sum(abs(therealy24[3,]-theesty24)) 
# 27.14973/3 

#overlap adjustment
therealy24_cp_s <- therealy24_cp
therealy24_cp_s[which(therealy24_cp == therealy24_cp[2])] <- therealy24_cp[which(therealy24_cp == therealy24_cp[2])]+
  seq(0,length(which(therealy24_cp== therealy24_cp[2]))*0.02-0.02, by=0.02)

therealy24_p_s <- therealy24_p
therealy24_p_s[which(therealy24_p == therealy24_p[2])] <- therealy24_p[which(therealy24_p == therealy24_p[2])]+
  seq(0,length(which(therealy24_p== therealy24_p[2]))*0.02-0.02, by=0.02)

range(therealy24_cp_s)
range(theesty24_c)
range(theesty24_u)
range(theesty24_d)

range(therealy24_p)
range(theesty24)
range(theesty24_u)
range(theesty24_d)

pdf("real_est_comp_both.pdf",height = 4,width = 8)
par(mfrow= c(1,2),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))

plot(x=therealy24_cp_s, y= theesty24_c, ylab = "Estimated log concentration (+/- psd) in cancer",
     xlab = "Real data mean (+/- sd) log concentration",main = "",col="black",bg="black",
     xlim=c(-5,10),ylim=c(-5,10),cex=0.5,pch=16)
for(i in 1:10){
  lines(x=c(therealy24_cd[i],therealy24_cu[i]), y= c(theesty24_c[i],theesty24_c[i]),
        col=2,bg=2,pch=24,cex=0.4)
  lines(x=c(therealy24_cp_s[i],therealy24_cp_s[i]), y=c(theesty24_cd[i],theesty24_cu[i]),
        col="blue",bg="blue",pch=24,cex=0.4)
}
lines(x=c(-30,32),y=c(-30,32),col="gray",lty="dashed")

plot(x=therealy24_p_s, y= theesty24, ylab = "Estimated log concentration (+/- psd) in noncancer",
     xlab = "Real data mean (+/- sd) log concentration",main = "",col="black",bg="black",
     xlim=c(-9,10),ylim=c(-9,10),cex=0.7,pch=16)
for(i in 1:10){
  lines(x=c(therealy24_d[i],therealy24_u[i]), y=c(theesty24[i],theesty24[i]),
        col=2,bg=2,pch=24,cex=0.4)
  lines(x=c(therealy24_p_s[i],therealy24_p_s[i]), y=c(theesty24_d[i],theesty24_u[i]),
        col="blue",bg="blue",pch=24,cex=0.4)
}
lines(x=c(-30,32),y=c(-30,32),col="gray",lty="dashed")

dev.off()
par(mfrow=c(1,1))

# the6flux <- read.table("est_flux.txt",header = F)
# v1 <- thepostmean[1] * the6flux[,2]/6.02214129e14
# v2 <- thepostmean[1] * the6flux[,3]/6.02214129e14
# v3 <- thepostmean[2] * the6flux[,4]/6.02214129e14
# v4 <- thepostmean[2] * the6flux[,5]/6.02214129e14
# v5 <- thepostmean[3] * the6flux[,4]/6.02214129e14
# v6 <- thepostmean[3] * the6flux[,5]/6.02214129e14
# 
# vv1 <- v1/sd(v1)
# vv2 <- v2/sd(v2)
# vv3 <- v3/sd(v3)
# vv4 <- v4/sd(v4)
# vv5 <- v5/sd(v5)
# vv6 <- v6/sd(v6)
# range(vv1,vv2,vv3,vv4,vv5,vv6)
# 
# pdf("the6flux_est.pdf",height = 5,width = 5)
# par(mfrow= c(1,1),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
# plot(vv5,type="l", ylab = "Estimated flux standardized by SD",
#      xlab = "Time",main = "",col=3,cex=0.4,ylim=c(0,20),xlim=c(0,400),lty=1)
# lines(vv6,col=6,cex=0.4,lty=2)
# legend(130, 5,
#        legend = c("Vf_Glc_SerPool1_13C0_13C0D000", "Vf_Glc_SerPool1_13C1_13C3D000"), 
#        col = c(3,6), 
#        lty = c(1,2), 
#        cex = 0.8, 
#        text.col = "black", 
#        box.lty = 1,box.col = "black",
#        horiz = F , 
#        inset = c(1, 1))
# dev.off()


# stability_nc <- sapply(1:52, function(j) geweke.diag(log(your_data_frame[8501:11000,j]),0.2,0.5)$z )
# pnorm(stability_nc)
# stable_check_nc <- sum((stability_nc>qnorm(0.025/52))*(stability_nc<qnorm(1-0.025/52)))/52
#non-bonferroni, 46

stability_y_conc <- sapply(1:10, function(j) geweke.diag((your_flux24[,j]),0.15,0.5)$z )
pnorm(stability_y_conc)
stable_check_y <- sum((stability_y_conc>qnorm(0.025/10))*(stability_y_conc<qnorm(1-0.025/10)))/10

stability_nc <- sapply(1:52, function(j) geweke.diag(log(your_data_frame[7501:13000,j]),0.15,0.5)$z )
pval <- (1-pnorm(abs(stability_nc)))*2
adj_pval <- p.adjust(pval, method = "bonferroni")
sum(adj_pval>0.05)

pnorm(stability_nc)
(1-pnorm(abs(stability_nc)))<0.05/2/52
stable_check_nc <- sum((stability_nc>qnorm(0.025/52))*(stability_nc<qnorm(1-0.025/52)))/52

stack_vec <- NA
flag0 <- seq(12500,13000,500)
flag1 <- seq(4000, 9000,500)
flag2 <- seq(0.1,0.45,0.05)
flag3 <- seq(0.15,0.5,0.05)
for (i1 in 1:length(flag0)){
  for (i2 in 1:length(flag1)){
    for(i3 in 1:length(flag2)){
      for(i4 in 1:length(flag3)){
    stability <- sapply(1:52, function(j)
     geweke.diag(log(your_data_frame[(flag1[i2]+1):flag0[i1],j]),flag2[i3],flag3[i4])$z )
    if(sum((pnorm(stability)>0.025/52)*(pnorm(stability)<1-0.025/52))>51){
      print(c(flag1[i2]+1,flag0[i1],flag2[i3],flag3[i4],sum((stability>qnorm(0.025/52))*(stability<qnorm(1-0.025/52)))))
    }
    }
  }
  }
}



#####################################################################################



######################################################################################
setwd("~/Desktop/HPCresults/Mac/holly/dim10/Thesis/final/htest_h0_3_low/")

file_names <- lapply(1:7, function(x) paste("prehtsample_curr39_",x,".txt", sep = ""))
your_data_frame <- do.call(rbind,lapply(file_names,read.table,header = T))
file_names2 <- lapply(1:7, function(x) paste("prehtsample2_curr39_",x,".txt", sep = ""))
your_data_frame2 <- do.call(rbind,lapply(file_names2,read.table,header = T))

mytrue_var <- log(optipar$V1[1:52])
mytrue_var2 <- mytrue_var
mytrue_var2[3] <- mytrue_var[3] - mytrue_var2[3]

png("39par1.png",height = 550,width = 350)
par(mfrow= c(2,1),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
i=1
plot(your_data_frame[,i],type = "l", ylab = "beta3",main="")
lines(x = c(0,26000),y = c(mytrue_var[3],mytrue_var[3]),col='red')
i=2
plot(your_data_frame[,i],type = "l", ylab = "beta22",main="")
lines(x = c(0,26000),y = c(mytrue_var[22],mytrue_var[22]),col='red')
dev.off()

png("39par2.png",height = 550,width = 350)
par(mfrow= c(2,1),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
i=1
plot(your_data_frame2[,i],type = "l", ylab = "difference of beta3",main="")
lines(x = c(0,26000),y = c(mytrue_var2[3],mytrue_var2[3]),col='red')
i=2
plot(your_data_frame2[,i],type = "l", ylab = "beta22",main="")
lines(x = c(0,26000),y = c(mytrue_var2[22],mytrue_var2[22]),col='red')
dev.off()
par(mfrow= c(1,1))

png("39diff.png",height = 350,width = 350)
par(mfrow= c(1,1),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
plot(your_data_frame2[,1],type = "l", ylab = "difference of beta3", main = "")
lines(x = c(0,26000),y = c(mytrue_var2[3],mytrue_var2[3]),col='red')
dev.off()
par(mfrow= c(1,1))


s <- seq(1,100,1)
data_list <- list(NA)
for(i in 1:length(s)){
  file_names <- lapply(1:6, function(x) paste("prehtsample2_curr",s[i],"_",x,".txt", sep = ""))
  your_data_frame <- do.call(rbind,lapply(file_names,read.table,header = T))
  data_list[[i]] <- your_data_frame
}
#credible interval
cov <- list(NA)
for(i in 1:length(s)){
  sample_CI <- hdi(data_list[[i]][-(1:1000),1], ci = 0.95)
  cov[[i]] <- c(sample_CI$CI_low, sample_CI$CI_high)
}

b3diff <- range(c(sapply(1:length(s), function(x) cov[[x]][1:2]),mytrue_var2[3]))

png("ht_h0_cov.png",height = 350,width = 350)
par(mfrow= c(1,1),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
plot(x=c(1,1), y=c(cov[[1]][1],cov[[1]][2]), type = "l",xlim = c(1,length(s)),
     ylim = c(b3diff[1]-0.001,b3diff[2]+0.001),xlab = "",ylab = "95% credible interval for beta3 difference")
for(i in 2:length(s)){
  lines(x=c(i,i),y=c(cov[[i]][1],cov[[i]][2]))
}
lines(x=c(0,101),y=c(mytrue_var2[3],mytrue_var2[3]), col="red")
dev.off()

mean(sapply(1:length(s), function(x) 
  cov[[x]][1]<=mytrue_var2[3]&&cov[[x]][2]>=mytrue_var2[3] )) 
#0.94

# p-value
dim <- 1
vec <- seq(0.001,1,0.002)
mypval <- function(i){
  cout <- rep(NA,length(vec))
  cout <- sapply(1:length(vec), function(j){
    credInt <- hdi(data_list[[i]][-(1:1000),dim], ci=vec[j])
    c1 <- credInt$CI_low
    c2 <- credInt$CI_high
    if(c1*c2<=0) return(min(abs(c1),abs(c2)))
    else return(NA)
  })
  index <- which(cout == min(cout,na.rm = T))
  return(1-min(vec[index]))
}

cov <- mclapply(1:length(s), mypval, mc.cores = 4)
# plot(log(data_list[[63]][,1]),type = "l")
# cov2[which(cov2==-Inf)] <- 0
range(cov)
png("ht_h0_pvalue.png",height = 350,width = 350)
par(mfrow= c(1,1),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
hist(unlist(cov), breaks = seq(0.0,1,0.05),freq = T,xlab = "credible value",
     ylab = "credible value histogram of 100 simulations",main="")
dev.off()
######################################################################################


######################################################################################
setwd("~/Desktop/HPCresults/Mac/holly/dim10/Thesis/final/htest_h1_3_0.35_low/")
setwd("~/Desktop/HPCresults/Mac/holly/dim10/Thesis/final/htest_h1_3_0.06_low/")

file_names <- lapply(1:7, function(x) paste("prehtsample_curr52_",x,".txt", sep = ""))
your_data_frame <- do.call(rbind,lapply(file_names,read.table,header = T))
file_names2 <- lapply(1:7, function(x) paste("prehtsample2_curr52_",x,".txt", sep = ""))
your_data_frame2 <- do.call(rbind,lapply(file_names2,read.table,header = T))

mytrue_var <- log(optipar$V1[1:52])

#for htest_h1_3_0.35_low
mytrue_var2 <- log(c(0.000465165222015355, 0.242665045280639, 0.00009790512, 0.518326766571717, 0.0002, 0.00003, 0.2,
                     0.000075, 0.075, 130.0,   0.6, 0.00018 , 0.075, 0.00005, 0.12 , 0.001,0.0008 , 0.00974237918290052,   0.0278492039049790 , 780.0 ,
                     0.6, 54.2492566486518, 0.03 , 0.95 ,   2.00461377104628, 0.398825587055585, 999999.999999978,  434.005733100667,  1.05698640506409,
                     2.24239831194499 , 0.551424353965021, 2.31291778438686, 44400.0000092611,  80.0, 44000.0001051987,  44000.0000615634, 0.0268575860192817,
                     0.000861408841917836,   0.107341967375275, 0.00756915225856982, 11.6618608584811,  0.499324532242705, 2.01627621832653 ,  0.0682897697038028,
                     4.57835999241956, 0.579450741578679, 42.0, 21.0, 120.0,1.2,  65.0,  0.72))

#for htest_h1_3_0.06_low
mytrue_var2 <- log(c(0.000465165222015355, 0.242665045280639, 0.0001308431, 0.518326766571717, 0.0002, 0.00003, 0.2,
                     0.000075, 0.075, 130.0,   0.6, 0.00018 , 0.075, 0.00005, 0.12 , 0.001,0.0008 , 0.00974237918290052,   0.0278492039049790 , 780.0 ,
                     0.6, 54.2492566486518, 0.03 , 0.95 ,   2.00461377104628, 0.398825587055585, 999999.999999978,  434.005733100667,  1.05698640506409,
                     2.24239831194499 , 0.551424353965021, 2.31291778438686, 44400.0000092611,  80.0, 44000.0001051987,  44000.0000615634, 0.0268575860192817,
                     0.000861408841917836,   0.107341967375275, 0.00756915225856982, 11.6618608584811,  0.499324532242705, 2.01627621832653 ,  0.0682897697038028,
                     4.57835999241956, 0.579450741578679, 42.0, 21.0, 120.0,1.2,  65.0,  0.72))

mytrue_var2[3] <- mytrue_var[3] - mytrue_var2[3]

png("ht_h1_52par1.png",height = 550,width = 350)
par(mfrow= c(2,1),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
i=1
plot(your_data_frame[,i],type = "l", ylab = "beta3",main="")
lines(x = c(0,26000),y = c(mytrue_var[3],mytrue_var[3]),col='red')
i=2
plot(your_data_frame[,i],type = "l", ylab = "beta22",main="")
lines(x = c(0,26000),y = c(mytrue_var[22],mytrue_var[22]),col='red')
dev.off()

png("ht_h1_52par2.png",height = 550,width = 350)
par(mfrow= c(2,1),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
i=1
plot(your_data_frame2[,i],type = "l", ylab = "difference of beta3",main="")
lines(x = c(0,26000),y = c(mytrue_var2[3],mytrue_var2[3]),col='red')
i=2
plot(your_data_frame2[,i],type = "l", ylab = "beta22",main="")
lines(x = c(0,26000),y = c(mytrue_var2[22],mytrue_var2[22]),col='red')
dev.off()
par(mfrow= c(1,1))

png("52diff.png",height = 350,width = 350)
par(mfrow= c(1,1),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
plot(your_data_frame2[,1],type = "l", ylab = "", main = "Difference of beta3")
lines(x = c(0,26000),y = c(mytrue_var2[3],mytrue_var2[3]),col='red')
dev.off()
par(mfrow= c(1,1))


s <- seq(1,100,1)
data_list <- list(NA)
for(i in 1:length(s)){
  file_names <- lapply(1:6, function(x) {
    paste("prehtsample2_curr",s[i],"_",x,".txt", sep = "")
  })
  your_data_frame <- do.call(rbind,lapply(file_names,read.table,header = T))
  data_list[[i]] <- your_data_frame
}
#credible interval
cov <- list(NA)
for(i in 1:length(s)){
  sample_CI <- hdi(data_list[[i]][-(1:1000),1], ci = 0.95)
  cov[[i]] <- c(sample_CI$CI_low, sample_CI$CI_high)
}

mytrue_var2[3]
diff <- range(c(sapply(1:length(s), function(x) cov[[x]][1:2]),mytrue_var2[3]))

png("ht_h1_cov.png",height = 350,width = 350)
par(mfrow= c(1,1),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
plot(x=c(1,1), y=c(cov[[1]][1],cov[[1]][2]), type = "l",xlim = c(1,length(s)),
     ylim = c(diff[1]-0.001,diff[2]+0.001),xlab = "",ylab = "95% credible interval for beta3 difference")
for(i in 2:length(s)){
  lines(x=c(i,i),y=c(cov[[i]][1],cov[[i]][2]))
}
lines(x=c(0,101),y=c(0,0), col="green")
lines(x=c(0,101),y=c(mytrue_var2[3],mytrue_var2[3]), col="red")
dev.off()

mean(sapply(1:length(s), function(x) 
  cov[[x]][1]<=mytrue_var2[3]&&cov[[x]][2]>=mytrue_var2[3] )) #0.95

mean(sapply(1:length(s), function(x) 
  cov[[x]][1]<=0&&cov[[x]][2]>=0 )) #0.06

#p-value
dim <- 1
vec <- seq(0.001,1,0.002)
mypval <- function(i){
  cout <- rep(NA,length(vec))
  cout <- sapply(1:length(vec), function(j){
    credInt <- hdi(data_list[[i]][-(1:1000),1], ci=vec[j])
    c1 <- credInt$CI_low
    c2 <- credInt$CI_high
    if(c1*c2<=0) return(min(abs(c1),abs(c2)))
    else return(NA)
  })
  if(sum(cout,na.rm=T)==0){return(0)}
  else{
    index <- which(cout == min(cout,na.rm = T))
    return(1-min(vec[index]))}
}

cov <- mclapply(1:length(s), mypval, mc.cores = 4)
# cov[which(cov==-Inf)] <- 0
range(cov)
png("ht_h1_pvalue.png",height = 350,width = 350)
par(mfrow= c(1,1),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
hist(unlist(cov), breaks = seq(0.0,1,0.05),freq = T,xlab = "credible value",
     ylab = "credible value histogram of 100 simulations", main="")
dev.off()
######################################################################################



######################################################################################
setwd("~/Desktop/HPCresults/Mac/holly/dim10/Thesis/final/htest_b10_h0_3/")

mytrue_var <- log(c(0.000461657750477753  ,  0.427503086639990  ,  0.000138933974778344 ,   0.912957765988904 ,   0.000302264234311588  ,
                    0.000264216988034605 ,   0.528433976069209  ,  0.000178629396143585  ,  0.0105686795213842 ,   389.997313257504 ,
                    3.65988124056536 ,   0.000422747180855367 ,   0.0496192767065513   , 7.44289150598270e-05 ,   0.0178629396143585  ,
                    0.0109460124414652 ,   0.0140352064043982 ,   0.00966781987350446  ,  0.0490562874450916  ,  830.998723620212,
                    17.0631986313531 ,   54.1720082557871,    1.25536770067575 ,   1.67337425755250 ,   2.03451173487994,    0.399124270787164 ,
                    1000000.00000000 ,   434.000150742915 ,   1.64330339153493  ,  2.22492836752176 ,   0.600557183706194 ,   2.29538774044506 ,
                    44400.0000000432  ,  89.9889179513641  ,  44000.0000016264 ,   43999.9999979280  ,  0.0266554754467594 ,   0.00151731008995339,
                    0.189003218774087  ,  0.00751136610783774  ,  11.6715973413360  ,  0.495498097191622 ,   2.51038222698186  ,
                    0.0677700081258078  ,  4.57755449879482,    0.579744450285708 ,   42.5463885164964 ,   21.0649543138687 ,   85.0037419701686  ,
                    120.001190410154 ,   65.0980133729533  ,  1.78629396143585))
mytrue_var2 <- mytrue_var
mytrue_var2[10] <- mytrue_var[10] - mytrue_var2[10]

s <- seq(1,20,1)
len <- c(22,24,25,26,21,20,21,20,25,25,
         26,23,23,27,22,19,23,25,24,23)
data_list <- list(NA)
for(i in 1:length(s)){
  file_names <- lapply(1:len[i], function(x) {
    paste("htsample2_curr_",s[i],"_",x,".txt", sep = "")
  })
  your_data_frame <- do.call(rbind,lapply(file_names,read.table,header = T))
  data_list[[i]] <- your_data_frame
}

# png("3par2.png",height = 1050,width = 750)
# par(mfrow= c(13,4),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
# for(i in 1:52){
#   plot(data_list[[3]][,i],type = "l", ylab = i,main="")
#   lines(x = c(0,26000),y = c(mytrue_var2[i],mytrue_var2[i]),col='red')
# }
# dev.off()
# par(mfrow= c(1,1))

#credible interval
cover <- list(NA)
for(i in 1:length(s)){
  sample_CI <- hdi(data_list[[i]][-(1:5000),10], ci = 0.95)
  cover[[i]] <- c(sample_CI$CI_low, sample_CI$CI_high)
}

mytrue_var2[10]
diff <- range(sapply(1:length(s), function(x) cover[[x]][1:2]),
              mytrue_var2[10])

# png("full_h0_10_cov.png",height = 350,width = 350)
# par(mfrow= c(1,1),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
# plot(x=c(1,1), y=c(cover[[1]][1],cover[[1]][2]), type = "l",xlim = c(1,length(s)),
#      ylim = c(diff[1]-0.001,diff[2]+0.001),xlab = "",ylab = "50% credible interval for beta10 difference")
# for(i in 2:length(s)){
#   lines(x=c(i,i),y=c(cover[[i]][1],cover[[i]][2]))
# }
# lines(x=c(0,101),y=c(mytrue_var2[10],mytrue_var2[10]), col="red")#lines(x=c(0,101),y=c(0,0), col="red")
# dev.off()

mean(sapply(1:length(s), function(x) 
  cover[[x]][1]<=0&&cover[[x]][2]>=0 )) #0.95

#p-value
dim <- 10
vec <- seq(0.01,1,0.01)
mypval <- function(i){
  cout <- rep(NA,length(vec))
  cout <- sapply(1:length(vec), function(j){
    credInt <- hdi(data_list[[i]][-(1:5000),dim], ci=vec[j])
    c1 <- credInt$CI_low
    c2 <- credInt$CI_high
    if(c1*c2>=0) return(min(abs(c1),abs(c2)))
    else return(NA)
  })
  if(sum(cout,na.rm=T)==0){return(0)}
  else{
    index <- which(cout == min(cout,na.rm = T))
    return(1-max(vec[index]))}
}

cov <- mclapply(1:length(s), mypval, mc.cores = 4)
# cov2 <- sapply(1:length(s), function(x) max(cov[[x]]))
# cov2[which(cov2==-Inf)] <- 0
range(cov)
# png("full_h0_10_pvalue.png",height = 350,width = 350)
# par(mfrow= c(1,1),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
# hist(unlist(cov), breaks = seq(0.0,1,0.05),freq = T,xlab = "credible value",
#      ylab="frequency for 20 simulations", main="")
# dev.off()


pdf("/Users/zhangxu/Desktop/HPCresults/Mac/holly/dim10/Thesis/final/sim_test_all.pdf",height = 12,width = 8)
par(mfrow= c(3,2),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
# pdf("full_h0_10_cov_p.pdf",height = 4,width = 8)
# par(mfrow= c(1,2),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
plot(x=c(1,1), y=c(cover[[1]][1],cover[[1]][2]), type = "l",xlim = c(1,length(s)),
     ylim = c(diff[1]-0.001,diff[2]+0.001),main="(A)",xlab = "Simulation #",ylab = "95% credible interval for beta10 difference")
for(i in 2:length(s)){
  lines(x=c(i,i),y=c(cover[[i]][1],cover[[i]][2]))
}
lines(x=c(0,101),y=c(mytrue_var2[10],mytrue_var2[10]), col="red")#lines(x=c(0,101),y=c(0,0), col="red")

hist(unlist(cov), breaks = seq(0.0,1,0.05),freq = T,main="(B)",xlab = "credible value", 
     ylab="frequency for 20 simulations")
# dev.off()
# par(mfrow=c(1,1))
######################################################################################



######################################################################################
setwd("~/Desktop/HPCresults/Mac/holly/dim10/Thesis/final/htest_b10_h1_02_3/")
#full dimension test
mytrue_var <- log(c(0.000461657750477753  ,  0.427503086639990  ,  0.000138933974778344 ,   0.912957765988904 ,   0.000302264234311588  ,
                    0.000264216988034605 ,   0.528433976069209  ,  0.000178629396143585  ,  0.0105686795213842 ,   389.997313257504 ,
                    3.65988124056536 ,   0.000422747180855367 ,   0.0496192767065513   , 7.44289150598270e-05 ,   0.0178629396143585  ,
                    0.0109460124414652 ,   0.0140352064043982 ,   0.00966781987350446  ,  0.0490562874450916  ,  830.998723620212,
                    17.0631986313531 ,   54.1720082557871,    1.25536770067575 ,   1.67337425755250 ,   2.03451173487994,    0.399124270787164 ,
                    1000000.00000000 ,   434.000150742915 ,   1.64330339153493  ,  2.22492836752176 ,   0.600557183706194 ,   2.29538774044506 ,
                    44400.0000000432  ,  89.9889179513641  ,  44000.0000016264 ,   43999.9999979280  ,  0.0266554754467594 ,   0.00151731008995339,
                    0.189003218774087  ,  0.00751136610783774  ,  11.6715973413360  ,  0.495498097191622 ,   2.51038222698186  ,
                    0.0677700081258078  ,  4.57755449879482,    0.579744450285708 ,   42.5463885164964 ,   21.0649543138687 ,   85.0037419701686  ,
                    120.001190410154 ,   65.0980133729533  ,  1.78629396143585 ))

mytrue_var2 <- log(c(0.000465165222015355, 0.242665045280639, 0.00009790512, 0.518326766571717, 0.0002, 0.00009, 0.2,
                     0.000075, 0.075, 319.3028,
                     0.6, 0.00018 , 0.075, 0.00005, 0.12 , 0.001,0.0008 , 0.00974237918290052,   0.0278492039049790 , 780.0 ,
                     0.6, 54.2492566486518, 0.03 , 0.95 ,   2.00461377104628, 0.398825587055585, 999999.999999978,  434.005733100667,  1.05698640506409,
                     2.24239831194499 , 0.551424353965021, 2.31291778438686, 44400.0000092611,  80.0, 44000.0001051987,  44000.0000615634, 0.0268575860192817,
                     0.000861408841917836,   0.107341967375275, 0.00756915225856982, 11.6618608584811,  0.499324532242705, 2.01627621832653 ,  0.0682897697038028,
                     4.57835999241956, 0.579450741578679, 42.0, 21.0, 120.0, 1.2,  65.0,  0.72))
mytrue_var2[10] <- mytrue_var[10] - mytrue_var2[10] #0.2

s <- seq(1,20,1)
len <- c(22,20,21,21,18,22,21,22,24,23,
         27,21,25,24,20,21,23,27,24,23)
data_list <- list(NA)
for(i in 1:length(s)){
  file_names <- lapply(1:len[i], function(x) {
    paste("htsample2_curr_",s[i],"_",x,".txt", sep = "")
  })
  your_data_frame <- do.call(rbind,lapply(file_names,read.table,header = T))
  data_list[[i]] <- your_data_frame
}

# png("18par2.png",height = 1050,width = 750)
# par(mfrow= c(13,4),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
# for(i in 1:52){
#   plot(data_list[[18]][,i],type = "l", ylab = i,main="")
#   lines(x = c(0,26000),y = c(mytrue_var2[i],mytrue_var2[i]),col='red')
# }
# dev.off()
# par(mfrow= c(1,1))

#credible interval
cover <- list(NA)
for(i in 1:length(s)){
  sample_CI <- hdi(data_list[[i]][-(1:5000),10], ci = 0.95)
  cover[[i]] <- c(sample_CI$CI_low, sample_CI$CI_high)
}

mytrue_var2[10]
diff <- range(sapply(1:length(s), function(x) cover[[x]][1:2]),mytrue_var2[10],0)

# png("full_h1_10_cov.png",height = 350,width = 350)
# par(mfrow= c(1,1),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
# plot(x=c(1,1), y=c(cover[[1]][1],cover[[1]][2]), type = "l",xlim = c(1,length(s)),
#      ylim = c(diff[1]-0.001,diff[2]+0.001),xlab = "",ylab = "95% credible interval for beta10 difference")
# for(i in 2:length(s)){
#   lines(x=c(i,i),y=c(cover[[i]][1],cover[[i]][2]))
# }
# lines(x=c(0,101),y=c(mytrue_var2[10],mytrue_var2[10]), col="red")#lines(x=c(0,101),y=c(0,0), col="red")
# lines(x=c(0,101),y=c(0,0), col="blue")
# dev.off()

mean(sapply(1:length(s), function(x) 
  cover[[x]][1]<=0&&cover[[x]][2]>=0 )) 

#credible interval, p-value
dim <- 10
vec <- seq(0.01,1,0.01)
mypval <- function(i){
  cout <- rep(NA,length(vec))
  cout <- sapply(1:length(vec), function(j){
    credInt <- hdi(data_list[[i]][-(1:5000),dim], ci=vec[j])
    c1 <- credInt$CI_low
    c2 <- credInt$CI_high
    if(c1*c2>=0) return(min(abs(c1),abs(c2)))
    else return(NA)
  })
  if(sum(cout,na.rm=T)==0){return(0)}
  else{
    index <- which(cout == min(cout,na.rm = T))
    return(1-max(vec[index]))}
}

cov <- mclapply(1:length(s), mypval, mc.cores = 4)
# cov2 <- sapply(1:length(s), function(x) max(cov[[x]]))
# cov2[which(cov2==-Inf)] <- 0
range(cov)
# png("full_h1_10_pvalue.png",height = 350,width = 350)
# par(mfrow= c(1,1),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
# hist(unlist(cov), breaks = seq(0.0,1,0.05),freq = T,xlab = "credible value", 
#      ylab="frequency for 20 simulations",main="")
# dev.off()


# pdf("full_h1_2_10_cov_p.pdf",height = 4,width = 8)
# par(mfrow= c(1,2),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
plot(x=c(1,1), y=c(cover[[1]][1],cover[[1]][2]), type = "l",xlim = c(1,length(s)),
     ylim = c(diff[1]-0.001,diff[2]+0.001),main="(C)",xlab = "Simulation #",ylab = "95% credible interval for beta10 difference")
for(i in 2:length(s)){
  lines(x=c(i,i),y=c(cover[[i]][1],cover[[i]][2]))
}
lines(x=c(0,101),y=c(mytrue_var2[10],mytrue_var2[10]), col="red")#lines(x=c(0,101),y=c(0,0), col="red")
lines(x=c(0,101),y=c(0,0), col="orange")

hist(unlist(cov), breaks = seq(0.0,1,0.05),freq = T,main="(D)",xlab = "credible value", 
     ylab="frequency for 20 simulations")
# dev.off()
# par(mfrow=c(1,1))
######################################################################################



######################################################################################
setwd("~/Desktop/HPCresults/Mac/holly/dim10/Thesis/final/htest_b10_h1_3/")
#full dimension test
mytrue_var <- log(c(0.000461657750477753  ,  0.427503086639990  ,  0.000138933974778344 ,   0.912957765988904 ,   0.000302264234311588  ,
                    0.000264216988034605 ,   0.528433976069209  ,  0.000178629396143585  ,  0.0105686795213842 ,   389.997313257504 ,
                    3.65988124056536 ,   0.000422747180855367 ,   0.0496192767065513   , 7.44289150598270e-05 ,   0.0178629396143585  ,
                    0.0109460124414652 ,   0.0140352064043982 ,   0.00966781987350446  ,  0.0490562874450916  ,  830.998723620212,
                    17.0631986313531 ,   54.1720082557871,    1.25536770067575 ,   1.67337425755250 ,   2.03451173487994,    0.399124270787164 ,
                    1000000.00000000 ,   434.000150742915 ,   1.64330339153493  ,  2.22492836752176 ,   0.600557183706194 ,   2.29538774044506 ,
                    44400.0000000432  ,  89.9889179513641  ,  44000.0000016264 ,   43999.9999979280  ,  0.0266554754467594 ,   0.00151731008995339,
                    0.189003218774087  ,  0.00751136610783774  ,  11.6715973413360  ,  0.495498097191622 ,   2.51038222698186  ,
                    0.0677700081258078  ,  4.57755449879482,    0.579744450285708 ,   42.5463885164964 ,   21.0649543138687 ,   85.0037419701686  ,
                    120.001190410154 ,   65.0980133729533  ,  1.78629396143585 ))

mytrue_var2 <- log(c(0.000465165222015355, 0.242665045280639, 0.00009790512, 0.518326766571717, 0.0002, 0.00009, 0.2,
                     0.000075, 0.075, 274.8265,   0.6, 0.00018 , 0.075, 0.00005, 0.12 , 0.001,0.0008 , 0.00974237918290052,   0.0278492039049790 , 780.0 ,
                     0.6, 54.2492566486518, 0.03 , 0.95 ,   2.00461377104628, 0.398825587055585, 999999.999999978,  434.005733100667,  1.05698640506409,
                     2.24239831194499 , 0.551424353965021, 2.31291778438686, 44400.0000092611,  80.0, 44000.0001051987,  44000.0000615634, 0.0268575860192817,
                     0.000861408841917836,   0.107341967375275, 0.00756915225856982, 11.6618608584811,  0.499324532242705, 2.01627621832653 ,  0.0682897697038028,
                     4.57835999241956, 0.579450741578679, 42.0, 21.0, 120.0, 1.2,  65.0,  0.72))
mytrue_var2[10] <- mytrue_var[10] - mytrue_var2[10] #0.35

s <- seq(1,20,1)
len <- c(24,21,23,21,19,23,21,21,23,23,
         26,21,25,22,19,21,21,28,23,25)
data_list <- list(NA)
for(i in 1:length(s)){
  file_names <- lapply(1:len[i], function(x) {
    paste("htsample2_curr_",s[i],"_",x,".txt", sep = "")
  })
  your_data_frame <- do.call(rbind,lapply(file_names,read.table,header = T))
  data_list[[i]] <- your_data_frame
}

# png("18par2.png",height = 1050,width = 750)
# par(mfrow= c(13,4),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
# for(i in 1:52){
#   plot(data_list[[18]][,i],type = "l", ylab = i,main="")
#   lines(x = c(0,26000),y = c(mytrue_var2[i],mytrue_var2[i]),col='red')
# }
# dev.off()
# par(mfrow= c(1,1))

#credible interval
cover <- list(NA)
for(i in 1:length(s)){
  sample_CI <- hdi(data_list[[i]][-(1:5000),10], ci = 0.95)
  cover[[i]] <- c(sample_CI$CI_low, sample_CI$CI_high)
}

mytrue_var2[10]
diff <- range(sapply(1:length(s), function(x) cover[[x]][1:2]),mytrue_var2[10],0)

# png("full_h1_10_cov.png",height = 350,width = 350)
# par(mfrow= c(1,1),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
# plot(x=c(1,1), y=c(cover[[1]][1],cover[[1]][2]), type = "l",xlim = c(1,length(s)),
#      ylim = c(diff[1]-0.001,diff[2]+0.001),xlab = "",ylab = "95% credible interval for beta10 difference")
# for(i in 2:length(s)){
#   lines(x=c(i,i),y=c(cover[[i]][1],cover[[i]][2]))
# }
# lines(x=c(0,101),y=c(mytrue_var2[10],mytrue_var2[10]), col="red")#lines(x=c(0,101),y=c(0,0), col="red")
# lines(x=c(0,101),y=c(0,0), col="blue")
# dev.off()

mean(sapply(1:length(s), function(x) 
  cover[[x]][1]<=0&&cover[[x]][2]>=0 )) 

#credible interval, p-value
dim <- 10
vec <- seq(0.01,1,0.01)
mypval <- function(i){
  cout <- rep(NA,length(vec))
  cout <- sapply(1:length(vec), function(j){
    credInt <- hdi(data_list[[i]][-(1:5000),dim], ci=vec[j])
    c1 <- credInt$CI_low
    c2 <- credInt$CI_high
    if(c1*c2>=0) return(min(abs(c1),abs(c2)))
    else return(NA)
  })
  if(sum(cout,na.rm=T)==0){return(0)}
  else{
    index <- which(cout == min(cout,na.rm = T))
    return(1-max(vec[index]))}
}

cov <- mclapply(1:length(s), mypval, mc.cores = 4)
# cov2 <- sapply(1:length(s), function(x) max(cov[[x]]))
# cov2[which(cov2==-Inf)] <- 0
range(cov)
# png("full_h1_10_pvalue.png",height = 350,width = 350)
# par(mfrow= c(1,1),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
# hist(unlist(cov), breaks = seq(0.0,1,0.05),freq = T,xlab = "credible value", 
#      ylab="frequency for 20 simulations",main="")
# dev.off()


# pdf("full_h1_10_cov_p.pdf",height = 4,width = 8)
# par(mfrow= c(1,2),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
plot(x=c(1,1), y=c(cover[[1]][1],cover[[1]][2]), type = "l",xlim = c(1,length(s)),
     ylim = c(diff[1]-0.001,diff[2]+0.001),main="(E)",xlab = "Simulation #",ylab = "95% credible interval for beta10 difference")
for(i in 2:length(s)){
  lines(x=c(i,i),y=c(cover[[i]][1],cover[[i]][2]))
}
lines(x=c(0,101),y=c(mytrue_var2[10],mytrue_var2[10]), col="red")#lines(x=c(0,101),y=c(0,0), col="red")
lines(x=c(0,101),y=c(0,0), col="orange")

hist(unlist(cov), breaks = seq(0.0,1,0.05),freq = T,main="(F)",xlab = "credible value", 
     ylab="frequency for 20 simulations")
dev.off()
par(mfrow=c(1,1))
######################################################################################




######################################################################################
setwd("~/Desktop/HPCresults/Mac/holly/dim10/Thesis/final/real_data/test/")

s <- c(2,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,27,28,
       33,48,49,50,51)+1
#                                                              22    
len <- c(24,50,29,41,33,36,42,33,23,29,36,21,27,35,36,27,32,35,20,38,29,
#            
         24,22,22,22,38)
startlen <- c(8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
              8,8,8,8,8)
startlen <- rep(9,length(len))
startlen <- rep(1,length(len)) #for convergence test
startlen <- rep(12,length(len)) #for final convergent choice, start=5501, end=7500, frac1=0.45, frac2=0.50
len <- rep(15,length(len))
data_list <- list(NA)
for(i in 1:length(s)){
    file_names <- lapply(startlen[i]:len[i], function(x) 
      paste("htsample2","_curr_",s[i]-1,"_",x,".txt", sep = ""))
  your_data_frame <- do.call(rbind,lapply(file_names,read.table,header = T))
  data_list[[i]] <- your_data_frame
  
}


write.table(save_data_frame[-(1:5000),1:26],
            '~/Desktop/HPCresults/Mac/holly/dim10/Thesis/final/real_data/cancer/results08dec2022/cppinput.csv',
            sep=',',row.names = FALSE,col.names = FALSE)
# mean(data_list[[1]][,6])

# png("par2.png",height = 750,width = 350)
# par(mfrow= c(3,1),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
# for (i in 1:length(s)){

# plot(data_list[[2]][,s[2]],type = "l", ylab = paste("difference of beta",s[2],sep=" "))
# lines(x=c(0,20000),y=c(0,0),col="red")

# }
# dev.off()
# par(mfrow= c(1,1))


#credible interval
cov <- list(NA)
for(i in 1:length(s)){
  sample_CI <- hdi(data_list[[i]][,s[i]], ci = 0.95)
  cov[[i]] <- c(sample_CI$CI_low, sample_CI$CI_high)
}

b22diff <- range(sapply(1:length(s), function(x) cov[[x]][1:2]))


# p-value
vec <- seq(0.005,1,0.005)
mypval <- function(i){
  cout <- rep(NA,length(vec))
  cout <- sapply(1:length(vec), function(j){
    credInt <- hdi(data_list[[i]][,s[i]], ci=vec[j])
    c1 <- credInt$CI_low
    c2 <- credInt$CI_high
    if(c1*c2>=0) return(min(abs(c1),abs(c2)))
    else return(NA)
  })
  if(sum(cout,na.rm=T)==0){return(1)}
  else{
    index <- which(cout == min(cout,na.rm = T))
    return(1-max(vec[index]))
    }
}

test_p <- sapply(1:length(s), mypval)

# png("real_test_cov.png",height = 350,width = 350)
# par(mfrow= c(1,1),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
# plot(x=c(1,1), y=c(cov[[1]][1],cov[[1]][2]), type = "l",xlim = c(1,length(s)),
#      ylim = c(b22diff[1]-0.001,b22diff[2]+0.001),xlab = "",ylab = "95% credible interval")
# for(i in 2:length(s)){
#   lines(x=c(i,i),y=c(cov[[i]][1],cov[[i]][2]))
# }
# lines(x=c(0,101),y=c(0,0), col="red")
# dev.off()

#1,6,11,14
s[c(1,10,20)]

s[c(10,18,20)]
len[c(10,18,20)]

library(rmeta)
library(grid)

tabletext <- cbind(c("Parameter",beta_formatted[s]),c("Credible Value",round(test_p,3)))

pdf("real_test_cov_horiz0411.pdf",height = 7,width = 8)
par(mfrow= c(1,1),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
forestplot(tabletext,
           c(NA,sapply(1:length(s),function(x) mean(cov[[x]]))),
           c(NA,sapply(1:length(s),function(x) cov[[x]][1])),
           c(NA,sapply(1:length(s),function(x) cov[[x]][2])),
           zero=0,
           boxsize = 0.2,graphwidth = unit(3,"inches"),
           is.summary = FALSE,
           clip=c(-5,4),xlim=c(-4,4),
           col=meta.colors(box="royalblue",line=c("darkblue"),summary = "black")
           )
dev.off()



#converge
stability <- sapply(1:length(s), function(j) geweke.diag(data_list[[j]][-(1:4000),s[j]],0.1,0.5)$z )
pnorm(stability)
stable_check_t <- sum((stability>qnorm(0.025))*(stability<qnorm(0.975)))/length(s)

#final
stability <- sapply(1:length(s), function(j) geweke.diag(data_list[[j]][1:2000,s[j]],0.45,0.5)$z )
pval <- (1-pnorm(abs(stability)))*2
adj_pval <- p.adjust(pval, method = "bonferroni")
sum(adj_pval>0.05)


stability <- sapply(1:length(s), function(j) heidel.diag(data_list[[j]][,s[j]],0.1,0.05)[1,1] )
stable_check_t <- sum(stability)/length(s) 


stack_vec <- NA
flag0 <- seq(6000,10000,500)
flag1 <- seq(4000, 8000,500)
flag2 <- seq(0.1,0.45,0.05)
flag3 <- seq(0.15,0.5,0.05)
for (i1 in 1:length(flag0)){
  for(i2 in 1:length(flag1)){
  for(i3 in 1:length(flag2)){
    for(i4 in 1:length(flag3)){
      if((flag1[i2]+1)<flag0[i1]){
      stability <- sapply(1:length(s), function(j)
        geweke.diag(data_list[[j]][(flag1[i2]+1):flag0[i1],s[j]],flag2[i3],flag3[i4])$z )
      if(sum((pnorm(stability)>0.025)*(pnorm(stability)<1-0.025))>25){
        stack_vec <- c(stack_vec,i1,i2,i3,i4,sum((stability>qnorm(0.025))*(stability<qnorm(1-0.025))))
        print(c(flag1[i2]+1,flag0[i1],flag2[i3],flag3[i4],sum((stability>qnorm(0.025))*(stability<qnorm(1-0.025)))))
      }
      }
    }
  }
 }
}
######################################################################################

