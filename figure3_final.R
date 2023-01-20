find.index<- function(Ti, times){
    return(sum(Ti>=times))
}
aa = read.table("flux_7501_10000.txt", sep='\t')
bb = read.table("flux_10001_13000.txt", sep='\t')
event_time = sort(unique(c(aa[,1], bb[,1])))
N = 100  
selected_time = event_time[seq(1,length(event_time),length.out=N)] ## To make the code run fast, I only selected a subset of 100 time points. But since the curves are very smooth, 100 time points seem to be sufficient
startpos1 = c(which(aa[,1]==0), nrow(aa))
startpos2 = c(which(bb[,1]==0), nrow(bb))
niter1 = length(startpos1)-1
niter2 = length(startpos2)-1
niter = niter1 + niter2
y_C0 = matrix(nrow=niter, ncol=N); y_C1 = matrix(nrow=niter, ncol=N)  # the concentrations at selected_time
for (i in 1:niter){
    if (i<=niter1){
        iteri = startpos1[i]:startpos1[i+1]
        index = sapply(selected_time, find.index, times=aa[iteri,1])
        y_C0[i,] = aa[iteri[index],2]; y_C1[i,] = aa[iteri[index],3]
    } else{
        iteri = startpos2[i-niter1]:startpos2[i-niter1+1]
        index = sapply(selected_time, find.index, times=bb[iteri,1])
        y_C0[i,] = bb[iteri[index],2]; y_C1[i,] = bb[iteri[index],3]
    }
}
all_mcmc_beta <- read.csv("~/Dropbox/XuZhang/manuscripts/BMC Bioinformatics/submit/data_for_figure_3/beta52_mcmc_iterations5500.csv", header=F) #saved under /Dropbox/XuZhang/manuscripts/BMC Bioinformatics/submit/data_for_figure_3

vv5 = matrix(nrow=niter, ncol=N); vv6 = matrix(nrow=niter, ncol=N)
for (i in 1:niter){
beta_i_th_mcmc <- colMeans(as.matrix(all_mcmc_beta))
    
v5 <- beta_i_th_mcmc[3] * y_C0[i,]/6.02214129e14   #Vf_Glc_SerPool1_13C0_13C0D000
v6 <- beta_i_th_mcmc[3] * y_C1[i,]/6.02214129e14   #Vf_Glc_SerPool1_13C1_13C3D000

vv5[i,] <- v5 
vv6[i,] <- v6 
}
upper_C0 = apply(vv5,2,quantile,0.975); lower_C0 = apply(vv5,2,quantile,0.025)
upper_C1 = apply(vv6,2,quantile,0.975); lower_C1 = apply(vv6,2,quantile,0.025)
est_C0 = colMeans(vv5); est_C1 = colMeans(vv6)

#figure for one MCMC iteration result of beta.
pdf("the6flux_est.pdf",height = 5,width = 5)
par(mfrow= c(1,1),mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
plot(x=selected_time, y=est_C0,type="l", ylab = "Estimated flux", # standardized by SD",
     xlab = "Time (min)",main = "",col=3,cex=0.4,xlim=c(0,15),
     lwd=2,lty=1)
lines(x=selected_time,y=est_C1,col=6,cex=0.4,lty=1,lwd=2)
lines(x=selected_time,y=upper_C0,col=3,cex=0.4,lty=2,lwd=1)
lines(x=selected_time,y=lower_C0,col=3,cex=0.4,lty=2,lwd=1)
lines(x=selected_time,y=upper_C1,col=6,cex=0.4,lty=2,lwd=1)
lines(x=selected_time,y=lower_C1,col=6,cex=0.4,lty=2,lwd=1)

legend(5, 0.002,
       legend = c("Vf_Glc_SerPool1_13C0_13C0D000", "Vf_Glc_SerPool1_13C1_13C3D000"), 
       col = c(3,6), 
       lty = c(1,1),
       lwd = c(2,2),
       cex = 0.8, 
       text.col = "black", 
       box.lty = 1,box.col = "black",
       horiz = F , 
       inset = c(1, 1))
dev.off()