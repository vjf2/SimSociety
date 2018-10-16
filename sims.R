#plot random points on 10 x 10 grid

library(colorspace)
library(SocGen)
library(reshape2)
library(fastnet)
library(igraph)
library(adehabitatHR)
library(rgdal)
library(rgeos)

options(stringsAsFactors = FALSE)

set.seed(2)

n<-100
d<-500

xnet<-net.holme.kim(n, 3, 1)
ixnet<-as.undirected(to.igraph(xnet))
layout<-layout.fruchterman.reingold(ixnet)
tkcs<-layout*2
el<-get.edgelist(ixnet)
el_pp<-apply(el, 2, function(x) paste0("id", x)) #add id prefix
colnames(el_pp)<-c("ID1", "ID2")
el_pp<-as.data.frame(el_pp)
el_pp$status<-"preference"

#starting coordinates
x<-x1<-tkcs[,1]
y<-y1<-tkcs[,2]

cols=rainbow(n, alpha=0.3)
windows()
plot(x1, y1, pch=21, bg=cols, xlim=c(-150,150), ylim=c(-150, 150), asp=1)

xl<-list(x1)
yl<-list(y1)

i=2

replicate(d-1, {
  x<<-xl[[i]]<<-x-rnorm(n, 0, abs(rnorm(1, 0.2, 0.5))) #gamma, inverse wishart, distributions
  y<<-yl[[i]]<<-y-rnorm(n, 0, abs(rnorm(1, 0.2, 0.5)))
  points(x, y, col=cols, pch=16)
  i<<-i+1
})

xs<-do.call("rbind", xl)
colnames(xs)<-paste0("id", 1:n)
ys<-do.call("rbind", yl)
colnames(ys)<-paste0("id", 1:n)

xm<-melt(xs)[,-1]
xm[,3]<-as.vector(ys)

names(xm)<-c("id", "x", "y")
xm[,"day"]<-rep(1:d, n)

text(xm[xm$day==1, c("x", "y", "id")], cex=1.2)

dm<-as.matrix(dist(xm[xm$day==1, c("x", "y")]))
diag(dm)<-NA
dimnames(dm)<-list(paste0("id",1:n), paste0("id", 1:n))
dm2<-reduce_pairs(mat2dat(as.matrix(dm), "distance"), "ID1", "ID2")
per5<-round(.05*nrow(dm2))

#limit preferences and avoids to top 1/3 of network
grid_buffer=20
x <- seq(min(xm[,"x"])-grid_buffer,max(xm[,"x"])+grid_buffer,by=1)
y <- seq(min(xm[,"y"])-grid_buffer,max(xm[,"y"])+grid_buffer,by=1)
xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE

hrxydata<-SpatialPointsDataFrame(xm[,c("x","y")],xm["id"])
uds_href<-kernelUD(hrxydata[,1],grid=xy)
VI<-kerneloverlaphr(uds_href, method = "VI", percent = 90)
diag(VI)<-NA
mean(VI, na.rm=TRUE)
vi<-mat2dat(VI, "VI")

###select pref and avoid

dm2<-merge_pairs(dm2, vi, "ID1", "ID2", all.x=FALSE, all.y=FALSE)
dm2<-reduce_pairs(dm2, "ID1", "ID2")
dmt<-merge_pairs(dm2, el_pp, "ID1", "ID2", all.x=TRUE, all.y=FALSE)
dmtr<-dmt[which(is.na(dmt$status)),]
#need to subtract for random
dmta<-dmtr[which(dmtr$VI>0.05),]

#prefs
# wrp<-sample(rownames(dmt), per5, prob=(1/dmt$VI)) #rows containing preferred pairs

#avoids
wra<-sample(rownames(dmta), per5)

#Pair Classifications
# pp<-dm2[rownames(dm2) %in% wrp, 1:2]
pp<-el_pp[,1:2]
ap<-dmta[rownames(dmta) %in% wra, 1:2]
rp<-dmtr[setdiff(rownames(dmtr), wra),1:2]

ppm<-pp
pptogether<-list()

starttime<-Sys.time()

# all_mats<-list()
# 
# for (j in 1:100){

xm2<-xm #original
xm2$group<-NA

for (k in 1:(d-1)){
  dy1<-xm2[xm2$day==k,]
  dy1$group<-NA
  
  distmat<-as.matrix(dist(dy1[,2:3]))
  dimnames(distmat)<-list(paste0("id",1:n), paste0("id", 1:n))  
  # distmat[lower.tri(distmat)]<-NA
  
  #need to fix indexing

  distmat[matrix(c(pp[,1], pp[,2], pp[,2], pp[,1]),ncol=2)]<-distmat[matrix(c(pp[,1], pp[,2], pp[,2], pp[,1]),ncol=2)]*0.1
  
  distmat[matrix(c(rp[,1], rp[,2], rp[,2], rp[,1]),ncol=2)]<-distmat[matrix(c(rp[,1], rp[,2], rp[,2], rp[,1]),ncol=2)]+5  
  
  distmat[matrix(c(ap[,1], ap[,2], ap[,2], ap[,1]),ncol=2)]<-distmat[matrix(c(ap[,1], ap[,2], ap[,2], ap[,1]),ncol=2)]*20
  # distmat[lower.tri(distmat)]<-t(distmat)[lower.tri(distmat)]
  foragers<-paste0("id", sample(1:n, n*0.15))
  distmat<-distmat[setdiff(rownames(distmat),foragers),
                   setdiff(colnames(distmat),foragers)]
  num_clust<-20
  
  #xk<-kmeans(distmat, num_clust, nstart=10)

  xc<-hclust(as.dist(distmat))
  xcc<-cutree(xc, k=num_clust)
  
  # dy1$group<-xk$cluster[match(dy1$id, names(xk$cluster))]
  
  dy1$group<-xcc[match(dy1$id, names(xcc))]
  
  # keep preferences in group
  ppm$group1<-dy1$group[match(ppm$ID1, dy1$id)]
  ppm$group2<-dy1$group[match(ppm$ID2, dy1$id)]
  implode<-which(ppm$group1==ppm$group2)
  pptogether[[k]]<-length(implode)
  
  # loners<-apply(apm[explode,1:2], 1, function(x) sample(c(x), 1))
  # loners<-sample(loners, floor(length(loners)*.5))
  # dy1$group[which(dy1$id %in% loners)]<-(num_clust+1):(num_clust+length(unique(loners)))

  dy2<-xm2[xm2$day==k+1,]
  
  gs<-unique(na.omit(dy1$group))
  
  for (i in gs){
     # dy2[which(dy1$group==i), "x"]<-mean(dy1[which(dy1$group==i), "x"])
     # dy2[which(dy1$group==i), "y"]<-mean(dy1[which(dy1$group==i), "y"])
     nr<-nrow(dy1[which(dy1$group==i),])
     dy1[which(dy1$group==i), "weight"]<-sample(c(2, rep(1, nr-1)), nr)
     
     dy2[which(dy1$group==i), "x"]<-weighted.mean(dy2[which(dy1$group==i), "x"], 
                                                  dy1[which(dy1$group==i), "weight"])
     dy2[which(dy1$group==i), "y"]<-weighted.mean(dy2[which(dy1$group==i), "y"], 
                                                  dy1[which(dy1$group==i), "weight"])
  }
  
  #need to make pp stickier
  
  comp<-cbind(dy1, dy2[,2:3])
  comp$distance<-apply(comp, 1, function(x) dist(rbind(c(x[2], x[3]),
                                                       c(x[7], x[8]))))
  names(comp)[c(7,8)]<-c("x1", "y1")
  
  comp$x1<-ifelse(comp$distance>20, (comp$x+comp$x1)/2, comp$x1)
  comp$y1<-ifelse(comp$distance>20, (comp$y+comp$y1)/2, comp$y1)
  
  #modify xm2
  xm2[xm2$day==k+1, c("x", "y")]<-dy2[,c("x", "y")]
  xm2[xm2$day==k+1, c("group")]<-dy1$group
  

} #end big loop (takes about 20 sec)

# xm2$group[which(is.na(xm2$group))]<-100:(nrow(xm2[which(is.na(xm2$group)),])+99)
# xm2$groupid<-paste0(xm2$day, "_",xm2$group)
# 
# mat<-simple_ratio(sightings=xm2, group_variable = "groupid", dates="day", IDs="id", symmetric = FALSE)
# 
# all_mats[[j]]<-mat
# }

endtime<-Sys.time()
endtime-starttime

xm2$group[which(is.na(xm2$group))]<-100:(nrow(xm2[which(is.na(xm2$group)),])+99)
xm2$groupid<-paste0(xm2$day, "_",xm2$group)
mean(table(xm2$groupid))

#mean association rate
mat<-simple_ratio(sightings=xm2, group_variable = "groupid", dates="day", IDs="id", symmetric = FALSE)
mean(mat, na.rm = TRUE)

AIs<-mat2dat(mat, "SRI")

pp2<-merge_pairs(pp, AIs, "ID1", "ID2", all.x=TRUE, all.y=FALSE)
rp2<-merge_pairs(rp, AIs, "ID1", "ID2", all.x=TRUE, all.y=FALSE)
ap2<-merge_pairs(ap, AIs, "ID1", "ID2", all.x=TRUE, all.y=FALSE)

mean(pp2$SRI)
mean(rp2$SRI)
mean(ap2$SRI)

windows()
par(mfrow=c(1,2))
hist(table(xm2$groupid), nclass=30)
hist(mat, nclass=20)

#time alone
sum(table(xm2$groupid)==1)/length(unique(xm2$groupid)) #not correct?


#create a grid on which to model animal home ranges
grid_buffer=20

x <- seq(min(xm2[,"x"])-grid_buffer,max(xm2[,"x"])+grid_buffer,by=0.5) # where resolution is the pixel size you desire. 100 is the smallest i would go, if you make it larger you'll get coarser resolution, but faster runtimes
y <- seq(min(xm2[,"y"])-grid_buffer,max(xm2[,"y"])+grid_buffer,by=0.5)

xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE

#create UDs for each animal and extract the href smoothing parameter (need to manually select h for boundary method)

hrxydata<-SpatialPointsDataFrame(xm2[,c("x","y")],xm2["id"])

uds_href<-kernelUD(hrxydata[,1],grid=xy)

# windows()
# par(mfrow=c(2,5), mar=c(0,0,0,0))
# for (i in 11:20) {image(uds_href[[i]])}

# hr<-getverticeshr(uds_href, percent=90)

# ka<-kernel.area(uds_href, percent=90, standardize = FALSE, unin="km", unout="km2")

VI<-kerneloverlaphr(uds_href, method = "VI", percent = 90)
mean(VI)

vi<-mat2dat(VI, "VI")

#check categories

pp2$status<-"preference"
rp2$status<-"random"
ap2$status<-"avoidance"

p2<-rbind(pp2, rp2, ap2)

check<-merge_pairs(p2, vi, "ID1", "ID2", all.x = TRUE, all.y=FALSE)

windows()
plot(SRI~VI, data=check, type="n")
points(SRI~VI, data=check[which(check$status=="random"),], col="grey")
points(SRI~VI, data=check[which(check$status=="avoidance"),], col="red")
points(SRI~VI, data=check[which(check$status=="preference"),], col="green")

abline(lm(SRI~VI+0, data=check[which(check$status=="random"),]), col="grey", lwd=2)
abline(lm(SRI~VI+0, data=check[which(check$status=="avoidance"),]), col="red", lwd=2)
abline(lm(SRI~VI+0, data=check[which(check$status=="preference"),]), col="green", lwd=2)


# save.image(file="simulation_output.RData")

#match VI to AI

#mean distance moved in a day

real_prefs<-check

pp2$status<-"preference"
rp2$status<-"random"
ap2$status<-"avoidance"

p2_rand<-rbind(pp2, rp2, ap2)

#compare p2_rand and check

#calculate the proportion of time that each individual is alone

group_size<-table(xm2$groupid)

xm2$group_size<-group_size[match(xm2$groupid, names(group_size))]

#individual group sizes
igs<-aggregate(group_size~id, data=xm2, mean)

ita<-aggregate(group_size~id, data=xm2[which(xm2$group_size==1),], sum)
ita$group_size<-ita$group_size/d

mean(igs$group_size)

# xm2<-xm2[!xm2$day %in% eod,]

#before
windows();plot(y~x, data=xm, col=as.factor(xm$id))
text(xm[xm$day==1, c("x", "y", "id")], cex=1.2)

#after
windows();plot(y~x, data=xm2, col=cols[as.factor(xm2$id)], bg=cols[as.factor(xm2$id)], pch=21, xlim=c(-150, 150), ylim=c(-150,150), asp=1)
text(xm2[xm2$day==1, c("x", "y", "id")], cex=1.2)

#individual home ranges
windows()
par(mfrow=c(2,5))
ids<-unique(xm2$id)
j=1

for (i in j:(j+9)){
  plot(xm$x,xm$y, type="n", xlim=c(-100, 100), ylim=c(-100,100), asp=1)
  points(xm2[xm2$id==ids[i],"x"],xm2[xm2$id==ids[i],"y"], col=NA, bg=cols[i], pch=21)
  points(0, 0, pch=3)
}
j=j+10

#original home ranges
windows()
par(mfrow=c(2,5))
ids<-unique(xm$id)
j=1

for (i in j:(j+9)){
  plot(xm$x,xm$y, type="n", xlim=c(-100, 100), ylim=c(-100,100), asp=1)
  points(xm[xm$id==ids[i],"x"],xm[xm$id==ids[i],"y"], col=NA, bg=cols[i], pch=21)
  points(0, 0, pch=3)
}
j=j+10

#network graph
library(igraph)

g<-graph.adjacency(mat, mode="undirected", diag=FALSE, weighted=TRUE)

windows();
plot(g, edge.width = edge_attr(g)$weight*20,
     edge.curved = rep(-.4,length(edge_attr(g)$weight)),
     vertex.color = cols, vertex.size=3)

i=1
windows()
today=xm2[xm2$day==i,]
plot(y~x, data=today, col=cols[as.factor(today$id)], bg=cols[as.factor(today$id)], pch=21)
text(today[, c("x", "y", "id")], cex=1.2)
i=i+1



#density 
#overlap of home ranges
#group size
#time spent alone
#how fast the groups change

#calculate output parameters

#average amount of home range overlap

#create spatial lines
#create UDs