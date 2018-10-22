#Some code for figure plotting

vcols<-sapply(degree(ixnet), function(x) adjustcolor("purple", alpha.f=(1-((1/x)*3))+0.1))

windows()

# pdf(file="preference_network.pdf")
par(mar=c(0,0,0,0))
plot(ixnet, 
     layout = layout,
     vertex.size = c(degree(ixnet)/1.5)+1,
     edge.arrow.size = 0, 
     edge.width = 2,
     edge.curved = rep(-.4,length(E(ixnet))),
     vertex.label = NA,
     vertex.color = vcols, 
     margin=c(0,0,0,0))
dev.off()

windows()
plot(xm[,1:2], col=cols[as.factor(xm[,3])], pch=20)
points(x1,y1, pch=15) #starting points 


#Core home ranges
ids <- unique(xm$id)
windows()
#pdf(file="core_home_ranges.pdf", width=7, height=3.5)
par(mar = c(5, 4, 2, 1))
layout(matrix(
  c(1, 1, 1, 1, 2, 4, 3, 5),
  nrow = 2,
  ncol = 4,
  byrow = FALSE
))
plot(
  xm[, 2:3],
  col = cols[as.factor(xm[, 1])],
  pch = 16,
  xlab = "Core Home Ranges - 500 steps",
  ylab = NA,
  yaxt = "n",
  cex.lab=1.25
)
axis(2, las = 1)
points(x1, y1, pch = 15)
legend(-32,
       47,
       pch = 15,
       legend = "Starting Locations",
       bty = "n",
       cex=1.25)
i = 1
for (i in i:(i + 3)) {
  plot(
    xm$x,
    xm$y,
    type = "n",
    xlim = c(-30, 30),
    ylim = c(-30, 30),
    asp = 1,
    ylab = NA,
    xlab = paste0("ID", i),
    yaxt = "n",
    cex.lab=1.25
  )
  axis(2, las = 1)
  points(xm[xm$id == ids[i], "x"],
         xm[xm$id == ids[i], "y"],
         col = NA,
         bg = rainbow(4, alpha = 0.4)[i],
         pch = 21)
  points(0, 0, pch = 3)
}

#match VI to AI

#mean distance moved in a day

real_prefs<-check

pp2$status<-"preference"
rp2$status<-"random"
ap2$status<-"avoidance"

p2_rand<-rbind(pp2, rp2, ap2)

#compare p2_rand and check

#calculate the proportion of time that each individual is alone

group_size <- table(xm2$groupid)

xm2$group_size <- group_size[match(xm2$groupid, names(group_size))]

#individual group sizes
igs <- aggregate(group_size ~ id, data = xm2, mean)

ita <-
  aggregate(group_size ~ id, data = xm2[which(xm2$group_size == 1), ], sum)
ita$group_size <- ita$group_size / d

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
  plot(xm$x,xm$y, type="n", xlim=c(-20, 20), ylim=c(-20,20), asp=1)
  points(xm2[xm2$id==ids[i],"x"],xm2[xm2$id==ids[i],"y"], col=NA, bg=cols[i], pch=21)
  points(0, 0, pch=3)
}
j=j+10

#original home ranges
windows()
par(mfrow = c(2, 5))
ids <- unique(xm$id)
j = 1

for (i in j:(j + 9)) {
  plot(
    xm$x,
    xm$y,
    type = "n",
    xlim = c(-20, 20),
    ylim = c(-20, 20),
    asp = 1
  )
  points(xm[xm$id == ids[i], "x"],
         xm[xm$id == ids[i], "y"],
         col = NA,
         bg = cols[i],
         pch = 21)
  points(0, 0, pch = 3)
}
j = j + 10

#network graph

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
