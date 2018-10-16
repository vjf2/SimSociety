#assign preferences in a more modular way

windows();plot(test, vertex.size=3, layout=lm)

erg<-erdos.renyi.game(100, 0.05, type = c("gnp"), directed = FALSE)

windows();plot(erg, vertex.size=3)

layout<-xm[xm$day==1,]
layout<-layout[match(V(test)$name, layout$id),]
lm<-as.matrix(layout[,c("x", "y")])


tkplot(test, vertex.size=3, layout=lm)

library(ergm)

g.sim <- simulate(network(100) ~ edges + mutual, coef=c(0, 0))
windows();plot(g.sim)

install.packages("fastnet")

xnet<-net.holme.kim(100, 3, 1)

mean(unlist(lapply(xnet, length)))

ixnet<-to.igraph(xnet)

ixnet<-as.undirected(ixnet)

tkplot(ixnet)

lec<-leading.eigenvector.community(ixnet)

modularity(ixnet, lec$membership)

tkc<-tk_coords(6)

tkcs<-scale(tkc)

tkcs<-tkcs*10

el<-get.edgelist(ixnet)

windows();plot(ixnet, vertex.size=3)

el_pp<-apply(el, 2, function(x) paste0("id", x)) #add id prefix

colnames(el_pp)<-c("ID1", "ID2")

el_pp<-as.data.frame(el_pp)

el_pp$status<-"preference"

