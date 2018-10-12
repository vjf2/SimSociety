#plot random points on 10 x 10 grid

library(colorspace)
library(SocGen)
library(reshape2)

options(stringsAsFactors = FALSE)

set.seed(2)

n<-200
d<-500

x<-x1<-runif(n, -50, 50)
y<-y1<-runif(n, -20, 20)

cols=rainbow(n, alpha=0.3)
windows()
plot(x1, y1, pch=21, bg=cols, xlim=c(-100, 100), ylim=c(-100, 100), asp=1)

xl<-list(x1)
yl<-list(y1)

i=2

replicate(d-1, {
  x<<-xl[[i]]<<-x-rnorm(n, 0, abs(rnorm(1, 0.5, 1))) #gamma, inverse wishart, distributions
  y<<-yl[[i]]<<-y-rnorm(n, 0, abs(rnorm(1, 0.5, 1)))
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

#preferences
#default attraction
#preference attraction

dm<-as.matrix(dist(xm[xm$day==1, c("x", "y")]))
diag(dm)<-NA
dimnames(dm)<-list(paste0("id",1:n), paste0("id", 1:n))
dm2<-reduce_pairs(mat2dat(as.matrix(dm), "distance"), "ID1", "ID2")
per5<-(.05*nrow(dm2))
rownames(dm2)<-1:nrow(dm2)

#limit preferences and avoids to top 1/3 of network

thresholds<-quantile(dm2$distance, probs=c(0.1, 0.33))

dmt<-dm2[dm2$distance<=thresholds[1],]
dmt2<-dm2[dm2$distance>thresholds[1] & dm2$distance<=thresholds[2],]

#prefs
wrp<-as.numeric(sample(rownames(dmt), per5, prob=(1/dmt$distance))) #rows containing preferred pairs

#avoids
wra<-as.numeric(sample(setdiff(rownames(dmt2), wrp), per5, prob=(1/dmt2[setdiff(rownames(dmt2), wrp), "distance"])))

#Pair Classifications
pp<-dm2[wrp, 1:2]
ap<-dm2[wra, 1:2]
rp<-dm2[-c(wrp, wra), 1:2]

apm<-ap


#now make it recursive
#change locations and groups on xm

xm2<-xm #backup
xm2$group<-NA

starttime<-Sys.time()

for (k in 1:(d-1)){
  dy1<-xm2[xm2$day==k,]
  dy1$group<-NA
  
  distmat<-as.matrix(dist(dy1[,2:3]))
  dimnames(distmat)<-list(paste0("id",1:n), paste0("id", 1:n))  
  distmat[upper.tri(distmat)]<-NA
  distmat[as.matrix(pp)]<-distmat[as.matrix(pp)]*0.01
  distmat[as.matrix(ap)]<-distmat[as.matrix(ap)]*10
  distmat[upper.tri(distmat)]<-t(distmat)[upper.tri(distmat)]  
  foragers<-paste0("id", sample(1:n, n*0.2))
  distmat<-distmat[setdiff(rownames(distmat),foragers),
                   setdiff(colnames(distmat),foragers)]
  num_clust<-35
  
# xk<-kmeans(distmat, num_clust, nstart=10)

  xc<-hclust(as.dist(distmat))
  xcc<-cutree(xc, k=num_clust)
  
  # dy1$group<-xk$cluster[match(dy1$id, names(xk$cluster))]
  
  dy1$group<-xcc[match(dy1$id, names(xcc))]
  
  #remove avoiders from group
  # apm$group1<-dy1$group[match(apm$ID1, dy1$id)]
  # apm$group2<-dy1$group[match(apm$ID2, dy1$id)]
  # explode<-which(apm$group1==apm$group2)
  # loners<-apply(apm[explode,1:2], 1, function(x) sample(c(x), 1))
  # loners<-sample(loners, floor(length(loners)*.5))
  # dy1$group[which(dy1$id %in% loners)]<-(num_clust+1):(num_clust+length(unique(loners)))
  # 
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
  
  comp$x1<-ifelse(comp$distance>50, (comp$x+comp$x1)/2, comp$x1)
  comp$y1<-ifelse(comp$distance>50, (comp$y+comp$y1)/2, comp$y1)
  
  
  #modify xm2
  xm2[xm2$day==k+1, c("x", "y")]<-dy2[,c("x", "y")]
  xm2[xm2$day==k+1, c("group")]<-dy1$group
} #end big loop (takes about 20 sec)

endtime<-Sys.time()
endtime-starttime

xm2$group[which(is.na(xm2$group))]<-100:(nrow(xm2[which(is.na(xm2$group)),])+99)
xm2$groupid<-paste0(xm2$day, "_",xm2$group)
mean(table(xm2$groupid))
windows();hist(table(xm2$groupid), probability = TRUE, nclass=30)

#mean association rate
mat<-simple_ratio(sightings=xm2, group_variable = "groupid", dates="day", IDs="id")
mat[upper.tri(mat)]<-NA
mean(mat, na.rm = TRUE)
windows();hist(mat, nclass=20)

AIs<-mat2dat(mat, "SRI")

pp2<-merge_pairs(pp, AIs, "ID1", "ID2", all.x=TRUE, all.y=FALSE)
rp2<-merge_pairs(rp, AIs, "ID1", "ID2", all.x=TRUE, all.y=FALSE)
ap2<-merge_pairs(ap, AIs, "ID1", "ID2", all.x=TRUE, all.y=FALSE)

mean(pp2$SRI)
mean(rp2$SRI)
mean(ap2$SRI)

#time alone
sum(table(xm2$groupid)==1)/length(unique(xm2$groupid)) #not correct?

#calculate the proportion of time that each individual is alone

group_size<-table(xm2$groupid)

xm2$group_size<-group_size[match(xm2$groupid, names(group_size))]

#individual group sizes
igs<-aggregate(group_size~id, data=xm2, mean)

ita<-aggregate(group_size~id, data=xm2[which(xm2$group_size==1),], sum)
ita$group_size<-ita$group_size/d

mean(igs$group_size)

#maybe need to prevent chaining
#have dolphins be social on even days?
#put a hard limit on distance moved in one step

#check before and after plots

#only take social days

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

g<-graph.adjacency(t(mat), mode="undirected", diag=FALSE, weighted=TRUE)

windows();
plot(g, edge.width = edge_attr(g)$weight*20,
     edge.curved = rep(-.4,length(edge_attr(g)$weight)),
     vertex.color = cols)

i=1
windows()
today=xm2[xm2$day==i,]
plot(y~x, data=today, col=cols[as.factor(today$id)], bg=cols[as.factor(today$id)], pch=21)
text(today[, c("x", "y", "id")], cex=1.2)
i=i+1

library(adehabitatHR)
library(rgdal)
library(rgeos)

#create a grid on which to model animal home ranges
grid_buffer=50

x <- seq(min(xm2[,"x"])-grid_buffer,max(xm2[,"x"])+grid_buffer,by=0.5) # where resolution is the pixel size you desire. 100 is the smallest i would go, if you make it larger you'll get coarser resolution, but faster runtimes
y <- seq(min(xm2[,"y"])-grid_buffer,max(xm2[,"y"])+grid_buffer,by=0.5)

xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE

#create UDs for each animal and extract the href smoothing parameter (need to manually select h for boundary method)

hrxydata<-SpatialPointsDataFrame(xm2[,c("x","y")],xm2["id"])

uds_href<-kernelUD(hrxydata[,1],grid=xy)

windows()
par(mfrow=c(2,5), mar=c(0,0,0,0))
for (i in 11:20) {image(uds_href[[i]])}

hr<-getverticeshr(uds_href, percent=90)

ka<-kernel.area(uds_href, percent=90, standardize = FALSE, unin="km", unout="km2")

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


save.image(file="simulation_output.RData")

#match VI to AI

#mean distance moved in a day


#density 
#overlap of home ranges
#group size
#time spent alone
#how fast the groups change

#calculate output parameters

#average amount of home range overlap

#create spatial lines
#create UDs

#mean group size

#need to fix this
