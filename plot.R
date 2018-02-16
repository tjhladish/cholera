d = read.table('output', header=F, sep=' ')
#y_max = max(unlist(d))
#d = output;
#rain = c((267.029+45.824+30.495), 26.68579, 169.00367, 275.74283, 418.06171, 91.20199, 140.94853, 140.94853, 456.00995, 39.24, 98.1, 215.82, 307.38, 18.086, 280.333, 253.204, 352.677, 39.537, 237.222, 263.58, 777.561, 68.464, 136.928, 179.718, 470.69, 11.488, 183.808, 333.152, 620.352, 164.065, 94.985, 198.605, 405.845, 30.348, 106.218, 227.61, 394.524, 10.091, 171.547, 272.457, 544.914, 40.552, 172.346, 263.588, 537.314, 16.198, 121.485, 251.069, 421.148, 110.964, 206.076, 467.634, 515.19, 13.714, 143.997, 198.853, 322.279)
obs_cholera = c(3,0,0,0,0,0,0,0,2,18,5,4,0,0,0,2,4,1,15,6,6,2,2,3,0,1,2,3,8,2,12,8,3,2,0,5,5,0,1,6,11,15,0,0,2,1,3,0,0,0,1,5,5,1,6,3,1,6,5,4,1,0,1,0,0,2,1,2,1,1,0,2,2,0,0,0,3,4,11,5,4,2,1,2,8,4,3,1,5,9,1,1,5,2,4,5,0,1,1,2,0,4,3,3,0,1,2,3,1,6,17,8,9,14,15,2,1,4,1,5,6,1,0,4,5,3,1,1,1,0,1,2,2,1,2,2,1,4,1,2,1,0,1,0,1,1,0,1,1,3,0,2,2,0,4,1,1,1,1,2,2,4,2,2,0,0,0,0)
mean_rain = c(74.892,8.407,14.214,28.9,46.386,76.279,43.7,60.8,94.15,102.729,207.964,196.021)
perturbation = 1;
#mean_rain = rep(mean_rain,14)
#rain = t(Vellore_monthly_rainfall);
#rain = t('Vellore_monthly_rainfall');
rain=read.table('Vellore_monthly_rainfall.txt', header=F, sep = ' ')
trans_rain=t(rain);
rain = mean(trans_rain)/(trans_rain+perturbation);
png('cholera.png', width=3000, height=4000, res=400)
par(mfrow=c(6,1))
par(oma=c(2,0,2,0))
par(mar=c(3,5,0,2))
#par(oma=c(0,0,0,0))

plot(seq(1,168,1), obs_cholera, type = 'l', col=1, xlab='', ylab='', axes=F)
mtext('Obs. incidence', 2, line=3)
axis(1, at=seq(0,168,12))
axis(2)

plot(seq(1,168,1), rain, type = 'l', col=1, xlab='', ylab='', axes=F)
mtext('Rainfall Deveiation', 2, line=3)
axis(1, at=seq(0,168,12))
axis(2)

plot(d[,1], type = 'l', col=1, xlab='', ylab='', axes=F)
mtext('S (susceptible)', 2, line=3)
axis(1, at=seq(0,168,12))
axis(2)

plot(d[,2], type = 'l', col=2, xlab='', ylab='', axes=F)
mtext('I (symptomatic)', 2, line=3)
axis(1, at=seq(0,168,12))
axis(2)

plot(d[,3], type = 'l', col=3, xlab='', ylab='', axes=F)
mtext('Y (asymptomatic)', 2, line=3)
axis(1, at=seq(0,168,12))
axis(2)

plot(d[,4], type = 'l', col=4, xlab='', ylab='', axes=F)
mtext('R (resistant)', 2, line=3)
axis(1, at=seq(0,168,12))
axis(2)

mtext('Time (months)', 1, outer=T)
dev.off()
