d = read.table('output', header=F, sep=' ')
#y_max = max(unlist(d))

rain = c((267.029+45.824+30.495), 26.68579, 169.00367, 275.74283, 418.06171, 91.20199, 140.94853, 140.94853, 456.00995, 39.24, 98.1, 215.82, 307.38, 18.086, 280.333, 253.204, 352.677, 39.537, 237.222, 263.58, 777.561, 68.464, 136.928, 179.718, 470.69, 11.488, 183.808, 333.152, 620.352, 164.065, 94.985, 198.605, 405.845, 30.348, 106.218, 227.61, 394.524, 10.091, 171.547, 272.457, 544.914, 40.552, 172.346, 263.588, 537.314, 16.198, 121.485, 251.069, 421.148, 110.964, 206.076, 467.634, 515.19, 13.714, 143.997, 198.853, 322.279)


png('cholera.png', width=3000, height=4000, res=400)
par(mfrow=c(5,1))
par(oma=c(2,0,2,0))
par(mar=c(3,5,0,2))
#par(oma=c(0,0,0,0))

plot(seq(0,168,3), rain, type = 'l', col=1, xlab='', ylab='', axes=F)
mtext('Quartly rainfall', 2, line=3)
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
