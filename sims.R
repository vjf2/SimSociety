#plot random points on 10 x 10 grid

library(colorspace)
library(SocGen)
library(reshape2)

options(stringsAsFactors = FALSE)

set.seed(2)

n<-50
d<-100

x<-x1<-runif(n, -50, 50)
y<-y1<-runif(n, -50, 50)

cols=rainbow(n, alpha=0.3)
windows()
plot(x1, y1, pch=21, bg=cols, xlim=c(-100, 100), ylim=c(-100, 100), asp=1)

xl<-list(x1)
yl<-list(y1)

i=2

replicate(d-1, {
  x<<-xl[[i]]<<-x-rnorm(n, 0, abs(rnorm(1, 0.1, 0.5))) #gamma, inverse wishart, distributions
  y<<-yl[[i]]<<-y-rnorm(n, 0, abs(rnorm(1, 0.1, 0.5)))
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
per5<-(.1*nrow(dm2))
rownames(dm2)<-1:nrow(dm2)

#limit preferences and avoids to top 1/3 of network

threshold<-quantile(dm2$distance, 0.33)

dm2<-dm2[dm2$distance<=threshold,]

#prefs
wrp<-as.numeric(sample(rownames(dm2), per5, prob=(1/dm2$distance))) #rows containing preferred pairs

#avoids
wra<-as.numeric(sample(setdiff(rownames(dm2), wrp), per5, prob=(1/dm2[setdiff(rownames(dm2), wrp), "distance"])))

dsts<-lapply(1:d, function(i) as.matrix(dist(xm[xm$day==i, c("x", "y")], diag=FALSE, upper=FALSE)))

q05<-quantile(unlist(dsts), 0.05)

#Pair Classifications
pp<-dm2[wrp, 1:2]
ap<-dm2[wra, 1:2]
rp<-dm2[-c(wrp, wra), 1:2]

#now make it recursive
#change locations and groups on xm

xm2<-xm #backup
xm2$group<-0

distrb<-ecdf(unlist(dsts))

#weird because no one groups on odd days
starttime<-Sys.time()

for (k in 1:(d-1)){

  dm<-as.matrix(dist(xm2[xm2$day==k, c("x", "y")]))
  diag(dm)<-NA
  dimnames(dm)<-list(paste0("id",1:n), paste0("id", 1:n))
  dm2<-reduce_pairs(mat2dat(as.matrix(dm), "distance"), "ID1", "ID2")
  rownames(dm2)<-1:nrow(dm2)
  
  #join probs for rand and prefer
  dm2$join_prob<-1-distrb(dm2$distance)
  dm2$weight<-rep(0.01, nrow(dm2))
  dm2$weight[wrp]<-0.4
  dm2$weight[wra]<-0
  dm2$weight[dm2$distance>=100]<-0
  
  dm2$grouped<-rbinom(nrow(dm2), 1, dm2$join_prob*dm2$weight)
  
  #remove 30% of individuals from foraging decisions
  
  foragers<-paste0("id", sample(1:n, n*.3))
  
  dm2$grouped<-ifelse((dm2$ID1 %in% foragers | dm2$ID2 %in% foragers), 
                      0, dm2$grouped)
  #extract groups
  
  g1<-dm2[which(dm2$grouped==1),]
  
  dy1<-xm2[xm2$day==k,]
  
  dy1$group<-0
  #group connected animals
  
  counter<-1
  
  apply(g1, 1, function(x){
    cstatus<-dy1[dy1$id %in% c(x[1], x[2]),"group"]
    dy1[dy1$id %in% c(x[1], x[2]),"group"]<<-ifelse(sum(cstatus)==0, counter, min(cstatus[cstatus > 0]))
    counter<<-max(dy1$group)+1
  })
  
  #move grouped animals to centroid
  
  dy2<-xm2[xm2$day==k+1,]
  
  gs<-unique(dy1$group[dy1$group>0])
  
  dy1$weight<-1
  
  for (i in gs){
    
    nr<-nrow(dy1[which(dy1$group==i),])
    
    dy1[which(dy1$group==i), "weight"]<-sample(c(2, rep(1, nr-1)), nr)
    
    dy2[which(dy1$group==i), "x"]<-weighted.mean(dy1[which(dy1$group==i), "x"], 
                                                 dy1[which(dy1$group==i), "weight"])
    dy2[which(dy1$group==i), "y"]<-weighted.mean(dy1[which(dy1$group==i), "y"], 
                                                 dy1[which(dy1$group==i), "weight"])
  }
  
  #modify xm2
  xm2[xm2$day==k+1, c("x", "y")]<-dy2[,c("x", "y")]
  xm2[xm2$day==k+1, c("group")]<-dy1$group
  
} #end big loop (takes about 20 sec)

endtime<-Sys.time()

endtime-starttime

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

for (i in 1:10){
  plot(xm$x,xm$y, type="n", xlim=c(-150, 150), ylim=c(-150,150), asp=1)
  points(xm2[xm2$id==ids[i],"x"],xm2[xm2$id==ids[i],"y"], col=NA, bg=cols[i], pch=21)
}

#network graph
library(igraph)

g<-graph.adjacency(mat, mode="undirected", diag=FALSE, weighted=TRUE)

windows();
plot(g, edge.width = edge_attr(g)$weight*20,
     edge.curved = rep(-.4,length(edge_attr(g)$weight)),
     vertex.color = cols)

i=1
windows()
today=xm2[xm2$day==i,]
plot(y~x, data=today, col=cols[as.factor(today$id)], bg=cols[as.factor(today$id)], pch=21, xlim=c(-50, 50), ylim=c(-50,50), asp=1)
text(today[, c("x", "y", "id")], cex=1.2)
i=i+1


#flip the timesteps
#some part of dolphins aren't part of social decisions each day


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
xm2$group[xm2$group==0]<-11:(nrow(xm2[which(xm2$group==0),])+10)
xm2$groupid<-paste0(xm2$day, "_",xm2$group)
mean(table(xm2$groupid))
windows();hist(table(xm2$groupid), probability = TRUE)

#mean association rate
mat<-simple_ratio(sightings=xm2, group_variable = "groupid", dates="day", IDs="id")
mat[upper.tri(mat)]<-NA
mean(mat, na.rm = TRUE)
windows();hist(mat, nclass=20)


AIs<-mat2dat(mat, "SRI")

pp<-dm2[wrp, 1:2]
ap<-dm2[wra, 1:2]
rp<-dm2[-c(wrp, wra), 1:2]

pp<-merge_pairs(pp, AIs, "ID1", "ID2", all.x=TRUE, all.y=FALSE)
rp<-merge_pairs(rp, AIs, "ID1", "ID2", all.x=TRUE, all.y=FALSE)
ap<-merge_pairs(ap, AIs, "ID1", "ID2", all.x=TRUE, all.y=FALSE)

mean(pp$SRI)
mean(rp$SRI)
mean(ap$SRI)


#modify the distance matrix directly 
windows()
plot(y~x, data=dy1, col=cols)
distmat<-as.matrix(dist(dy1[,2:3]))
dimnames(distmat)<-list(paste0("id",1:n), paste0("id", 1:n))
distmat[upper.tri(distmat)]<-NA

distmat[as.matrix(pp)]<-distmat[as.matrix(pp)]*0.1
distmat[as.matrix(ap)]<-distmat[as.matrix(ap)]*10

distmat[upper.tri(distmat)]<-t(distmat)[upper.tri(distmat)]
 


num_clust<-floor(n/4.5)
xk<-kmeans(distmat, num_clust, nstart=3)
text(dy1$x, dy1$y, xk$cluster)



apply(pp, 1, function(x)
{
  lines(dy1[dy1$id %in% c(x[1], x[2]),"x"], dy1[dy1$id %in% c(x[1], x[2]), c("y")],
        lty=2)
})

apply(ap, 1, function(x)
{
  lines(dy1[dy1$id %in% c(x[1], x[2]),"x"], dy1[dy1$id %in% c(x[1], x[2]), c("y")],
        lty=2, col="red")
})

