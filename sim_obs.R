fulldata<-read.csv("dolphin_sim.csv", row.names = 1)

#observation areas
#approx 50 x 50 study site
#2500 study site

windows()
plot(fulldata$x, fulldata$y)
rect(-10, -10, 20, 10, border="red", lwd=2)

#xleft, ybottom, xright, ytop 

t1<-c(-10, 0, 0, 10) 
t2<-c(0, 0, 10, 10)
t3<-c(10, 0, 20, 10)
t4<-c(-10, -10, 0, 0)
t5<-c(0, -10, 10, 0)
t6<-c(10, -10, 20, 0)

rect(t1[1], t1[2], t1[3], t1[4], col=adjustcolor("blue", alpha.f=0.4))
rect(t2[1], t2[2], t2[3], t2[4], col=adjustcolor("green", alpha.f=0.4))
rect(t3[1], t3[2], t3[3], t3[4], col=adjustcolor("yellow", alpha.f=0.4))
rect(t4[1], t4[2], t4[3], t4[4], col=adjustcolor("orange", alpha.f=0.4))
rect(t5[1], t5[2], t5[3], t5[4], col=adjustcolor("red", alpha.f=0.4))
rect(t6[1], t6[2], t6[3], t6[4], col=adjustcolor("purple", alpha.f=0.4))



sample_days<-seq(1, d, 5)

all_obs<-list()

for(i in 1:length(sample_days)){
  
  day<-fulldata[fulldata$day==i,]
  squares<-sample(paste0("t",1:6), 3)
  day_obs<-list()
  
  for (j in 1:length(squares)){
    cs<-get(squares[j])
    obs<-day[day$x>cs[1] & day$x<cs[3] &
               day$y>cs[2] & day$y<cs[4],]
    day_obs[[j]]<-obs
  }
  all_obs[[i]]<-do.call("rbind", day_obs)
}

all_obs<-do.call("rbind", all_obs)

points(all_obs$x, all_obs$y, pch=16, col="red")
length(unique(all_obs$groupid))

points(all_obs[all_obs$day==6,"x"], all_obs[all_obs$day==6,"y"], pch=16, col="blue")

#make grid and home ranges

#make daily mcps

daily_xydata<-SpatialPointsDataFrame(all_obs[,c("x","y")],all_obs["day"])
mcps<-mcp(daily_xydata[,1], percent=100, unin=c("km"), unout=c("km2"))

#Add buffer, make sure whole area is covered

buff_days<-gBuffer(mcps, byid=TRUE,width=1)

#Number of animals in study
n<-length(unique(all_obs$id))

#Number of survey days
d<-length(unique(all_obs$day))
dates<-sort(unique(all_obs$day))

#make fake schedule 
dcounts<-table(all_obs$id)
dolphins<-names(dcounts[dcounts>=5])

fast_avail<-data.frame(dolphin_id=dolphins, entry=1, depart=500)

alive<-Vectorize(FUN=function(r,c) 
  isTRUE(r>=fast_avail$entry[which(fast_avail$dolphin_id==c)] 
         & r<=fast_avail$depart[which(fast_avail$dolphin_id==c)]))

schedule<-outer(as.numeric(dates), dolphins, FUN=alive) #takes 2 min to run

matnames<-list(as.character(dates),dolphins)
dimnames(schedule)<-matnames

#need to think about area

dolphin_density_per_km<-dim(all_obs)[1]/gArea(buff_days)

gArea(buff_days)

areakm<-gArea(buff_days, byid=TRUE)
numdol<-round(areakm*dolphin_density_per_km)
numdol<-ifelse(numdol<=1, 2, numdol) 


num_sim=100 #number of simulations to run

udsgdf <- as(estUDm2spixdf(uds_href),"SpatialGridDataFrame")
fullgrid(udsgdf)<-FALSE

#Set up cluster for parallelization
#Timing depends on number of cores available
#3 cores - 1000 sims in 50 min

library(parallel)
library(pbapply)

cl<-makeCluster(detectCores()-1)
clusterEvalQ(cl, library(sp))
clusterEvalQ(cl, library(SocGen))
clusterExport(cl, c("d", "buff_days", "udsgdf", "schedule", "num_sim", "numdol"))

starttime<-Sys.time()

nest_days<-pblapply(seq_len(d), FUN=function(i){
  
  bound<-buff_days[i,]
  nd<-numdol[i]
  dailygrid<-udsgdf[bound,,drop=TRUE] 
  cellsize<-udsgdf@grid@cellsize
  probweights<-colSums(dailygrid@data, na.rm=TRUE)
  probweights<-probweights[names(probweights) %in% colnames(schedule)[schedule[i,]==TRUE]]
  dc<-coordinates(dailygrid)
  dgdf<-dailygrid@data
  holder<-replicate(num_sim, fast_random_points(probweights = probweights, 
                                                nd = nd, 
                                                dc = dc,
                                                dgdf = dgdf,
                                                gridrad = cellsize/2), 
                    simplify=FALSE)
  return(holder) }, cl=cl)

# })

endtime<-Sys.time()

stopCluster(cl)

endtime-starttime #check run time #last run was 1.8 hours

sim_surveys<-sapply(1:num_sim, function(i) lapply(nest_days, "[[", i), simplify = FALSE)

#save(sim_surveys, file="sim_surveys100_20180621.RData")
#load("sim_surveys100_20180621.RData")

rm(nest_days)

#Calculate mean group size in real data (xydata3)

mean_group_size<-mean(table(all_obs$groupid))
sd(table(all_obs$groupid))


dpd<-sapply(unique(all_obs$day), function(i) length(all_obs$day[which(all_obs$day==i)]))

gpd<-sapply(unique(all_obs$day), function(i) length(unique(all_obs$groupid[which(all_obs$day==i)])))

prob_per_day<-data.frame(gpd,dpd)

prob_per_day$gs<-prob_per_day$dpd/prob_per_day$gpd

mod<-lm(gpd~dpd+areakm, data=prob_per_day) 

#Assign groups, (should take about 35 min for 1000 sims of full data)

bd_id<-buff_days$id

kfinal<-lapply(seq_len(num_sim), function(w) {
  single_sim<-sim_surveys[[w]]
  counter=1
  each_days_assoc<-lapply(single_sim, function(x1) {
    
    num_clust<-round(mod$coefficients[1]+(nrow(x1)*mod$coefficients[2]))
    distmat<-dist(x1[,1:2])
    size<-attr(distmat, "Size")*(attr(distmat, "Size")-1)/2
    num_clust<-ifelse(size < num_clust, size, num_clust)
    
    xc<-hclust(as.dist(distmat))
    xcc<-cutree(xc, k=num_clust)
    dayAssocK<-data.frame(IDs=x1[,3],Group=as.numeric(xcc))
    
    #assign by gprox
    
    
    dayAssocK$Permutation<-c(rep(as.character(bd_id[counter]),
                                 nrow(dayAssocK)))
    dayAssocK$Permutation<-dayAssocK$Permutation
    counter <<- counter + 1
    return(dayAssocK)
  })
  eda<-do.call("rbind", each_days_assoc)
  eda$id<-paste0(eda$Permutation,"_",eda$Group)
  return(eda)
}) #end kfinal *apply

endtime<-Sys.time()

endtime-starttime #check run time #2.7 hours not in parallel

random_group_sizes<-lapply(kfinal, function(x) mean(table(x$id)))

random_group_sizes_sd<-lapply(kfinal, function(x) sd(table(x$id)))
mean(unlist(random_group_sizes))


kfinal<-lapply(kfinal, function(x) {names(x)<-c("dolphin_id", "DayGroup", "Date", "observation_id");x})
kfinal<-lapply(kfinal, function(x) {x[,"observation_id"]<-as.numeric(as.factor(x[,"observation_id"]));x})
# kfinal<-lapply(kfinal, function(x) {x[,"Date"]<-as.numeric(x[,"Date"]);x})
kfinal<-lapply(kfinal, function(x) {x[,"dolphin_id"]<-as.character(x[,"dolphin_id"]);x})

masked_randoms<-lapply(kfinal, function(x){
  simple_ratio(sightings=x,
               group_variable="observation_id", 
               dates="Date", 
               IDs="dolphin_id", 
               symmetric=FALSE)})

obs_network<-simple_ratio(sightings=all_obs,
                            group_variable="groupid", 
                            dates="day", 
                            IDs="id", 
                            symmetric=FALSE)


mn<-mat2dat(obs_network, "realSRI")

mr<-lapply(masked_randoms, mat2dat)

for(i in 1:length(mr)) {names(mr[[i]])[3]<-paste0("randSRI", i)}

mr<-c(list(mn), mr)

fdata<-Reduce(function(x, y) merge(x, y, all.x = T), mr) #about 2 minutes

#merge removes anything missing, individuals that didn't overlap have NaN and get removed in merge

fdata$quantile975<-apply(fdata, 1, function(x) quantile(as.numeric(x[4:13]), 0.975, na.rm=TRUE))

fdata$meanSRI<-apply(fdata[,4:13],1, mean, na.rm=TRUE)

fdata$affiliation975<-ifelse(fdata$realSRI>fdata$quantile975, TRUE, FALSE)
fdata<-fdata[,grep("rand", names(fdata), invert=TRUE)]

#merge in real preferences

results<-merge_pairs(check, fdata, "ID1", "ID2", all.x=TRUE, all.y=FALSE)


windows()
plot(realSRI~VI, data=results, col=as.factor(affiliation975), bg=as.factor(status), pch=21)


table(results$status, results$affiliation975)

###ok, try timestamp swapping method

#make in format of sim_surveys

sim_surveys<-list()

for (i in 1:num_sim) {

swapped_obs<-split(all_obs[,1:4], as.factor(all_obs$id)) 

swapped_obs<-lapply(swapped_obs, function(x){ x$day<-sample(x$day) 
                                              return(x)})

swapped_obs<-do.call("rbind", swapped_obs)

sim_surveys1<-split(swapped_obs, f=as.factor(swapped_obs$day))

sim_surveys[[i]]<-sim_surveys1

}


kfinal<-lapply(seq_len(num_sim), function(w) {
  single_sim<-sim_surveys[[w]]
  counter=1
  each_days_assoc<-lapply(single_sim, function(x1) {
    
    num_clust<-round(mod$coefficients[1]+(nrow(x1)*mod$coefficients[2]))
    distmat<-dist(x1[,2:3])
    size<-attr(distmat, "Size")*(attr(distmat, "Size")-1)/2
    num_clust<-ifelse(size < num_clust, size, num_clust)
    
    xc<-hclust(as.dist(distmat))
    xcc<-cutree(xc, k=num_clust)
    dayAssocK<-data.frame(IDs=x1[,1],Group=as.numeric(xcc))
    
    
    dayAssocK$Permutation<-c(rep(as.character(bd_id[counter]),
                                 nrow(dayAssocK)))
    dayAssocK$Permutation<-dayAssocK$Permutation
    counter <<- counter + 1
    return(dayAssocK)
  })
  eda<-do.call("rbind", each_days_assoc)
  eda$id<-paste0(eda$Permutation,"_",eda$Group)
  return(eda)
}) #end kfinal *apply

endtime<-Sys.time()

endtime-starttime #check run time #2.7 hours not in parallel

random_group_sizes<-lapply(kfinal, function(x) mean(table(x$id)))

random_group_sizes_sd<-lapply(kfinal, function(x) sd(table(x$id)))
mean(unlist(random_group_sizes))


kfinal<-lapply(kfinal, function(x) {names(x)<-c("dolphin_id", "DayGroup", "Date", "observation_id");x})
kfinal<-lapply(kfinal, function(x) {x[,"observation_id"]<-as.numeric(as.factor(x[,"observation_id"]));x})
# kfinal<-lapply(kfinal, function(x) {x[,"Date"]<-as.numeric(x[,"Date"]);x})
kfinal<-lapply(kfinal, function(x) {x[,"dolphin_id"]<-as.character(x[,"dolphin_id"]);x})

masked_randoms<-lapply(kfinal, function(x){
  simple_ratio(sightings=x,
               group_variable="observation_id", 
               dates="Date", 
               IDs="dolphin_id", 
               symmetric=FALSE)})

obs_network<-simple_ratio(sightings=all_obs,
                          group_variable="groupid", 
                          dates="day", 
                          IDs="id", 
                          symmetric=FALSE)


mn<-mat2dat(obs_network, "realSRI")

mr<-lapply(masked_randoms, mat2dat)

for(i in 1:length(mr)) {names(mr[[i]])[3]<-paste0("randSRI", i)}

mr<-c(list(mn), mr)

fdata<-Reduce(function(x, y) merge(x, y, all.x = T), mr) #about 2 minutes

#merge removes anything missing, individuals that didn't overlap have NaN and get removed in merge

fdata$quantile975<-apply(fdata, 1, function(x) quantile(as.numeric(x[4:13]), 0.975, na.rm=TRUE))

fdata$meanSRI<-apply(fdata[,4:13],1, mean, na.rm=TRUE)

fdata$affiliation975<-ifelse(fdata$realSRI>fdata$quantile975, TRUE, FALSE)
fdata<-fdata[,grep("rand", names(fdata), invert=TRUE)]

results_swap<-merge_pairs(check, fdata, "ID1", "ID2", all.x=TRUE, all.y=FALSE)


windows()
plot(realSRI~VI, data=results_swap, col=as.factor(affiliation975), bg=as.factor(status), pch=21)


table(results_swap$status, results_swap$affiliation975)

table(results$status, results$affiliation975)

results$numID1<-dcounts[match(results$ID1, names(dcounts))]
results$numID2<-dcounts[match(results$ID2, names(dcounts))]

falsepos<-results[which(results$status=="random" & results$affiliation975==TRUE),c(1:5, 110:112)]

mean(falsepos$numID1+falsepos$numID2)

mean(results$numID1+results$numID2)
