setwd("./results/")
rm(list=ls(all=TRUE))
source("../InferenceFunctions.R")

load("K3_L100.RData")
rm(cur.compound, reaction.right.borders, reaction.left.borders, lines, network.file, reaction.borders, i, lines.count)

num.gibbs.samples = length(sample.results$Z)
est.theta = estimate.theta(sample.results$Z, K, 0.01, num.gibbs.samples)

est.theta.list = list()
for (i in 1:num.gibbs.samples) est.theta.list[[i]] = estimate.theta(sample.results$Z, K, 0.01, i)

setwd("plots/")

library("plotrix")
# coloring for separating carnivores, omnivores and herbivores.
# to get numbers for samples use sample.names 
sample.colors = rep("", num.samples)
names(sample.colors) = sample.names
sample.colors[c(1, 4, 5, 8, 11, 14, 16, 17, 18, 19, 21, 22, 23, 26, 27, 28, 30, 34, 36, 37, 38)] = "green"
sample.colors[c(2, 9, 15, 20, 24, 25, 29)] = "red"
sample.colors[c(3, 6, 7, 10, 12, 13, 31, 32, 33, 35)] = "black"

pdf("est-theta-list.pdf")
for (i in 1:num.gibbs.samples) triax.plot(est.theta.list[[i]], col.symbols=sample.colors, pch=20)
dev.off()

pdf(file="est.theta-colored.pdf", width=10, height=10)
triax.plot(est.theta, main="", col.axis="blue", cex.axis=1, cex.ticks=1.5, cex=1.5, align.labels=TRUE, show.grid=TRUE, col.grid="gray", lty.grid=1, pch=19, col.symbols=sample.colors)
dev.off()

est.reaction.membership = estimate.reaction.membership(sample.results$Y, Reactions, L, num.gibbs.samples)
est.reaction.membership.normalized = est.reaction.membership
for (l in 1:L) est.reaction.membership.normalized[,l] = est.reaction.membership[,l] / sum(est.reaction.membership[,l])

est.phi = estimate.phi(sample.results$Z, sample.results$Y, K, L, 0.01, num.gibbs.samples)
for (k in 1:K) est.phi[k,] = est.phi[k,] / sum(est.phi[k,])
colnames(est.phi) = 1:L
rownames(est.phi) = 1:K

setwd("../subnetworks/")
subsystems.reactions = list()
for (l in 1:L) subsystems.reactions[[l]] = est.reaction.membership.normalized[order(-est.reaction.membership.normalized[,l]), l]
for(l in 1:L) {
  sink(paste("subnetwork-", l,"-reactions.txt", sep=""));
  cat(paste( names(subsystems.reactions[[l]]), subsystems.reactions[[l]], "\n", sep=" "));
  sink();
}


major.subsystems.group1 = which(est.phi[1,] > 2.0/L)
major.subsystems.group2 = which(est.phi[2,] > 2.0/L)
major.subsystems.group3 = which(est.phi[3,] > 2.0/L)
major.subsystems = sort(union( union(major.subsystems.group1, major.subsystems.group2), major.subsystems.group3))
major.subsystems.colors = rainbow(length(major.subsystems))
names(major.subsystems.colors) = major.subsystems

setwd("../plots/")

###ribbon plot ####
# preparing phi matrix for plotting (only major subsystems)
est.phi.col.max = apply(est.phi,2,max)
est.phi.for.plotting = matrix(0, nrow=nrow(est.phi[,major.subsystems]), ncol=ncol(est.phi[,major.subsystems])*2)
est.phi.for.plotting[,seq(1,by=2,length.out=ncol(est.phi[,major.subsystems]))] = est.phi[,major.subsystems]
for (i in 1:nrow(est.phi)) est.phi.for.plotting[i,seq(2,by=2,length.out=ncol(est.phi[,major.subsystems]))] = est.phi.col.max[major.subsystems] - est.phi[i,major.subsystems] + 0.01


# plotting phi matrix, all bars are colored
phi.row.colors = rep(0, ncol(est.phi.for.plotting))
phi.row.colors[seq(1,by=2,length.out=ncol(est.phi[,major.subsystems]))] = rainbow(ncol(est.phi[,major.subsystems])) 
phi.row.colors[seq(2,by=2,length.out=ncol(est.phi[,major.subsystems]))] = 0
pdf(file="phi-barplot.pdf", width=80, height=10)
barplot(t(est.phi.for.plotting), horiz=TRUE, density=100,col=phi.row.colors, names.arg=c("Type 1","Type 2","Type 3"), width=c(8,1,1), space=0.2, border=NA, axes=FALSE, las=2)
text(c(0, cumsum(est.phi.col.max[major.subsystems[-length(major.subsystems)]]+0.01))+0.01, 10.0, paste("S", major.subsystems), cex=4.5, pos=3, srt=90)
barplot(t(est.phi.for.plotting), horiz=TRUE, density=100,col=phi.row.colors, names.arg=c("Type 1","Type 2","Type 3"), width=c(1,8,1), space=0.2, border=NA, axes=FALSE, las=2 )
text(c(0, cumsum(est.phi.col.max[major.subsystems[-length(major.subsystems)]]+0.01))+0.01, 10.0, paste("S", major.subsystems), cex=4.5, pos=3, srt=90)
barplot(t(est.phi.for.plotting), horiz=TRUE, density=100,col=phi.row.colors, names.arg=c("Type 1","Type 2","Type 3"), width=c(1,1,8), space=0.2, border=NA, axes=FALSE, las=2 )
text(c(0, cumsum(est.phi.col.max[major.subsystems[-length(major.subsystems)]]+0.01))+0.01, 10.0, paste("S", major.subsystems), cex=4.5, pos=3, srt=90)
dev.off()

# plotting phi matrix with bars only for main type colored
est.phi.col.max = apply(est.phi,2,max)
est.phi.for.plotting = matrix(0, nrow=nrow(est.phi[,major.subsystems]), ncol=ncol(est.phi[,major.subsystems])*2)
est.phi.for.plotting[,seq(1,by=2,length.out=ncol(est.phi[,major.subsystems]))] = est.phi[,major.subsystems]
for (i in 1:nrow(est.phi)) est.phi.for.plotting[i,seq(2,by=2,length.out=ncol(est.phi[,major.subsystems]))] = est.phi.col.max[major.subsystems] - est.phi[i,major.subsystems] + 0.01
est.phi.for.plotting = t(est.phi.for.plotting)
# trick taken from stackoverflow.com/questions/5169265/
tmp = est.phi.for.plotting
est.phi.for.plotting = matrix(0, nrow=nrow(tmp)*ncol(tmp), ncol=ncol(tmp))
est.phi.for.plotting[cbind(1:(nrow(tmp)*ncol(tmp)), rep(1:ncol(tmp), each=nrow(tmp)))] = tmp
phi.row.colors = rep(0, ncol(est.phi.for.plotting))
phi.row.colors[seq(1,by=2,length.out=ncol(est.phi[,major.subsystems]))] = rainbow(ncol(est.phi[,major.subsystems])) 
phi.row.colors[seq(2,by=2,length.out=ncol(est.phi[,major.subsystems]))] = 0
phi.row.colors.inactive = rep(0, ncol(est.phi.for.plotting))
phi.row.colors.inactive[seq(1,by=2,length.out=ncol(est.phi[,major.subsystems]))] = gray(seq(0.6,0.7,length.out=ncol(est.phi[,major.subsystems]))) #"black"
phi.row.colors.inactive[seq(2,by=2,length.out=ncol(est.phi[,major.subsystems]))] = 0
pdf(file="phi-barplot-gray.pdf", width=80, height=10)
barplot(est.phi.for.plotting, horiz=TRUE, density=100,col=c(phi.row.colors,phi.row.colors.inactive,phi.row.colors.inactive), names.arg=c("Type 1","Type 2","Type 3"), cex.names=5, width=c(8,1,1), space=0.2, border=NA, axes=FALSE, las=2, beside=FALSE)
text(c(0, cumsum(est.phi.col.max[major.subsystems[-length(major.subsystems)]]+0.01))+0.01, 10.0, paste("S", major.subsystems), cex=5.0, pos=3, srt=90)
barplot(est.phi.for.plotting, horiz=TRUE, density=100,col=c(phi.row.colors.inactive,phi.row.colors,phi.row.colors.inactive), names.arg=c("Type 1","Type 2","Type 3"), cex.names=5, width=c(1,8,1), space=0.2, border=NA, axes=FALSE, las=2, beside=FALSE )
text(c(0, cumsum(est.phi.col.max[major.subsystems[-length(major.subsystems)]]+0.01))+0.01, 10.0, paste("S", major.subsystems), cex=5.0, pos=3, srt=90)
barplot(est.phi.for.plotting, horiz=TRUE, density=100,col=c(phi.row.colors.inactive,phi.row.colors.inactive,phi.row.colors), names.arg=c("Type 1","Type 2","Type 3"), cex.names=5, width=c(1,1,8), space=0.2, border=NA, axes=FALSE, las=2, beside=FALSE )
text(c(0, cumsum(est.phi.col.max[major.subsystems[-length(major.subsystems)]]+0.01))+0.01, 10.0, paste("S", major.subsystems), cex=5.0, pos=3, srt=90)
dev.off()


# plotting differences between one metabotype and the other two. We include all three comparisons.  
library("shape")
subsystems.to.plot = major.subsystems
difference.range = log( (est.phi[1,subsystems.to.plot] / pmax(est.phi[2,subsystems.to.plot], est.phi[3,subsystems.to.plot]) +1))
names(difference.range) = subsystems.to.plot
difference.range.color = ( max(difference.range) - difference.range ) / (max(difference.range) - min(difference.range))
pdf("subnetwork-differences_1-major.pdf") # , width=12, height=6)
barplot(difference.range, col=gray(difference.range.color), names.arg=paste("S", subsystems.to.plot, sep=" "), cex.names=1, las=2 ) 
#colorlegend(col=gray(seq(max(difference.range.color), min(difference.range.color),by=-0.001)), c(min(difference.range), max(difference.range)), zlevels=3, posy=c(0.6,0.95),posx=c(0.8,0.82), digit=2)
colorlegend(col=gray(seq(max(difference.range.color), min(difference.range.color),length.out=length(difference.range.color))), c(min(difference.range), max(difference.range)), zlevels=3, posy=c(0.6,0.95),posx=c(0.8,0.82), digit=2)
dev.off()

subsystems.to.plot = 1:L
difference.range = log((1+ est.phi[1,subsystems.to.plot] / pmax(est.phi[2,subsystems.to.plot], est.phi[3,subsystems.to.plot])))
names(difference.range) = subsystems.to.plot
difference.range.color = ( max(difference.range) - difference.range ) / (max(difference.range) - min(difference.range))
pdf("subnetwork-differences_1.pdf" , width=12, height=6)
barplot(difference.range, col=gray(difference.range.color), names.arg=paste("S", subsystems.to.plot, sep=" "), cex.names=0.5, las=2 ) 
#colorlegend(col=gray(seq(max(difference.range.color), min(difference.range.color),by=-0.001)), c(min(difference.range), max(difference.range)), zlevels=3, posy=c(0.6,0.95),posx=c(0.8,0.82), digit=2)
colorlegend(col=gray(seq(max(difference.range.color), min(difference.range.color),length.out=length(difference.range.color))), c(min(difference.range), max(difference.range)), zlevels=3, posy=c(0.6,0.95),posx=c(0.8,0.82), digit=2)
dev.off()


subsystems.to.plot = major.subsystems
difference.range = log( (est.phi[2,subsystems.to.plot] / pmax(est.phi[1,subsystems.to.plot], est.phi[3,subsystems.to.plot]) +1))
names(difference.range) = subsystems.to.plot
difference.range.color = ( max(difference.range) - difference.range ) / (max(difference.range) - min(difference.range))
pdf("subnetwork-differences_2-major.pdf") # , width=12, height=6)
barplot(difference.range, col=gray(difference.range.color), names.arg=paste("S", subsystems.to.plot, sep=" "), cex.names=1, las=2 ) 
#colorlegend(col=gray(seq(max(difference.range.color), min(difference.range.color),by=-0.001)), c(min(difference.range), max(difference.range)), zlevels=3, posy=c(0.6,0.95),posx=c(0.8,0.82), digit=2)
colorlegend(col=gray(seq(max(difference.range.color), min(difference.range.color),length.out=length(difference.range.color))), c(min(difference.range), max(difference.range)), zlevels=3, posy=c(0.6,0.95),posx=c(0.8,0.82), digit=2)
dev.off()

subsystems.to.plot = 1:L
difference.range = log((1+ est.phi[2,subsystems.to.plot] / pmax(est.phi[1,subsystems.to.plot], est.phi[3,subsystems.to.plot])))
names(difference.range) = subsystems.to.plot
difference.range.color = ( max(difference.range) - difference.range ) / (max(difference.range) - min(difference.range))
pdf("subnetwork-differences_2.pdf" , width=12, height=6)
barplot(difference.range, col=gray(difference.range.color), names.arg=paste("S", subsystems.to.plot, sep=" "), cex.names=0.5, las=2 ) 
#colorlegend(col=gray(seq(max(difference.range.color), min(difference.range.color),by=-0.001)), c(min(difference.range), max(difference.range)), zlevels=3, posy=c(0.6,0.95),posx=c(0.8,0.82), digit=2)
colorlegend(col=gray(seq(max(difference.range.color), min(difference.range.color),length.out=length(difference.range.color))), c(min(difference.range), max(difference.range)), zlevels=3, posy=c(0.6,0.95),posx=c(0.8,0.82), digit=2)
dev.off()

subsystems.to.plot = major.subsystems
difference.range = log( (est.phi[3,subsystems.to.plot] / pmax(est.phi[2,subsystems.to.plot], est.phi[1,subsystems.to.plot]) +1))
names(difference.range) = subsystems.to.plot
difference.range.color = ( max(difference.range) - difference.range ) / (max(difference.range) - min(difference.range))
pdf("subnetwork-differences_3-major.pdf") # , width=12, height=6)
barplot(difference.range, col=gray(difference.range.color), names.arg=paste("S", subsystems.to.plot, sep=" "), cex.names=1, las=2 ) 
#colorlegend(col=gray(seq(max(difference.range.color), min(difference.range.color),by=-0.001)), c(min(difference.range), max(difference.range)), zlevels=3, posy=c(0.6,0.95),posx=c(0.8,0.82), digit=2)
colorlegend(col=gray(seq(max(difference.range.color), min(difference.range.color),length.out=length(difference.range.color))), c(min(difference.range), max(difference.range)), zlevels=3, posy=c(0.6,0.95),posx=c(0.8,0.82), digit=2)
dev.off()

subsystems.to.plot = 1:L
difference.range = log((1+ est.phi[3,subsystems.to.plot] / pmax(est.phi[2,subsystems.to.plot], est.phi[1,subsystems.to.plot])))
names(difference.range) = subsystems.to.plot
difference.range.color = ( max(difference.range) - difference.range ) / (max(difference.range) - min(difference.range))
pdf("subnetwork-differences_3.pdf" , width=12, height=6)
barplot(difference.range, col=gray(difference.range.color), names.arg=paste("S", subsystems.to.plot, sep=" "), cex.names=0.5, las=2 ) 
#colorlegend(col=gray(seq(max(difference.range.color), min(difference.range.color),by=-0.001)), c(min(difference.range), max(difference.range)), zlevels=3, posy=c(0.6,0.95),posx=c(0.8,0.82), digit=2)
colorlegend(col=gray(seq(max(difference.range.color), min(difference.range.color),length.out=length(difference.range.color))), c(min(difference.range), max(difference.range)), zlevels=3, posy=c(0.6,0.95),posx=c(0.8,0.82), digit=2)
dev.off()

save.image(file="../K3_L100_Analysis.RData")
