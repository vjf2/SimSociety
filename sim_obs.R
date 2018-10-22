library(SocGen)
library(adehabitatHR)
library(rgdal)
library(rgeos)

options(stringsAsFactors = FALSE)

setwd("C:/Users/froug/Desktop/PizaRocaReview")

set.seed(4)

n<-100
d<-2000

fulldata<-read.csv("sim_res_to_test2000.csv")

answer_key<-read.csv("sim_cats2000.csv")

#observation areas
#approx 50 x 50 study site
#2500 study site

fullsite<-c(-30, -20, 30, 20)

windows()
plot(fulldata$x, fulldata$y, col=adjustcolor("black", alpha.f=0.3), yaxt="n", ylab=NA, xlab="Subsampled Observations", pch=16)
axis(2, las=1)
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


sample_days<-seq(2, 2000, 20)

sample_days<-1:d

all_obs<-list()
all_coords<-list()

for(i in sample_days){
  
  day<-fulldata[fulldata$day==i,]
  squares<-c("t1", sample(paste0("t",2:6), 3)) #bias toward square 1
  day_obs<-list()
  day_coords<-list()
    
  for (j in 1:length(squares)){
    cs<-get(squares[j])
    obs<-day[day$x>cs[1] & day$x<cs[3] &
               day$y>cs[2] & day$y<cs[4],]
    day_obs[[j]]<-obs
    day_coords[[j]]<-matrix(c(cs[1], cs[1], cs[3], cs[3],
                              cs[2], cs[4], cs[2], cs[4],
                              i,i,i,i), ncol=3, byrow=FALSE)
  }
  all_obs[[i]]<-do.call("rbind", day_obs)
  all_coords[[i]]<-do.call("rbind", day_coords)
}

all_obs<-do.call("rbind", all_obs)
all_coords<-do.call("rbind", all_coords)

points(all_obs$x, all_obs$y, pch=16, col="red")
length(unique(all_obs$groupid))

points(all_obs[all_obs$day==11,"x"], all_obs[all_obs$day==11,"y"], pch=16, col="blue")

#make grid and home ranges

#test out with full data set

# all_obs<-fulldata[which(fulldata$x<fullsite[3] & fulldata$x>fullsite[1] &
                          # fulldata$y<fullsite[4] & fulldata$y>fullsite[2]),]

#make daily mcps

daily_xydata<-SpatialPointsDataFrame(all_coords[,1:2],data.frame(day=all_coords[,3]))

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

all_obs<-all_obs[which(all_obs$id %in% dolphins),]

fast_avail<-data.frame(dolphin_id=dolphins, entry=1, depart=2000)

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
numdol<-ifelse(numdol>length(dolphins), length(dolphins), numdol)

num_sim=1000 #number of simulations to run

grid_buffer=5
x <- seq(min(all_obs[,"x"])-grid_buffer,max(all_obs[,"x"])+grid_buffer,by=0.5) 
y <- seq(min(all_obs[,"y"])-grid_buffer,max(all_obs[,"y"])+grid_buffer,by=0.5)
xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE

#create UDs for each animal and extract the href smoothing parameter (need to manually select h for boundary method)

hrxydata<-SpatialPointsDataFrame(all_obs[,c("x","y")],all_obs["id"])

uds_href<-kernelUD(hrxydata[,1],grid=xy)

udsgdf <- as(estUDm2spixdf(uds_href),"SpatialGridDataFrame")
fullgrid(udsgdf)<-FALSE

gridrad<-udsgdf@grid@cellsize[1]/2

#Set up cluster for parallelization
#Timing depends on number of cores available
#3 cores - 1000 sims in 50 min

library(parallel)
library(pbapply)

cl<-makeCluster(detectCores()-1)
clusterEvalQ(cl, library(sp))
clusterEvalQ(cl, library(SocGen))
clusterExport(cl, c("d", "buff_days", "udsgdf", "schedule", "num_sim", "numdol", "gridrad"))

starttime<-Sys.time()

nest_days<-pblapply(seq_len(d), FUN=function(i){
  
  bound<-buff_days[i,]
  nd<-numdol[i]
  dailygrid<-udsgdf[bound,,drop=TRUE] 
  probweights<-colSums(dailygrid@data, na.rm=TRUE)
  probweights<-probweights[names(probweights) %in% colnames(schedule)[schedule[i,]==TRUE]]
  dc<-coordinates(dailygrid)
  dgdf<-dailygrid@data
  holder<-replicate(num_sim, fast_random_points(probweights = probweights, 
                                                nd = nd, 
                                                dc = dc,
                                                dgdf = dgdf,
                                                gridrad = gridrad), 
                    simplify=FALSE)
  return(holder) }, cl=cl)

endtime<-Sys.time()

stopCluster(cl)

endtime-starttime #check run time #last run was 1.8 hours

sim_surveys<-sapply(1:num_sim, function(i) lapply(nest_days, "[[", i), simplify = FALSE)

rm(nest_days)

#save(sim_surveys, file="sim_surveysUDd1000_2000.RData")

#Assign groups, (should take about 35 min for 1000 sims of full data)

mean_group_size<-mean(table(all_obs$groupid))

bd_id<-buff_days$id

starttime<-Sys.time()

kfinal_swap<-group_assign(data=sim_surveys, id="id", xcoord ="x", ycoord="y", time="day", group_size = mean_group_size)

endtime<-Sys.time()

endtime-starttime #check run time #2.7 hours not in parallel

random_group_sizes<-lapply(kfinal, function(x) mean(table(x$id)))
mean(unlist(random_group_sizes))


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

fdata$quantile975<-apply(fdata, 1, function(x) quantile(as.numeric(x[4:(num_sim+3)]), 1, na.rm=TRUE))

fdata$meanSRI<-apply(fdata[,4:(num_sim+3)],1, mean, na.rm=TRUE)

fdata$maxSRI<-apply(fdata[,4:num_sim+3],1, max, na.rm=TRUE)

fdata$affiliation975<-ifelse(fdata$realSRI>fdata$quantile975, TRUE, FALSE)
fdata<-fdata[,grep("rand", names(fdata), invert=TRUE)]

#merge in real preferences

results<-merge_pairs(answer_key, fdata, "ID1", "ID2", all.x=TRUE, all.y=FALSE)

results$dcounts1<-dcounts[match(results$ID1, names(dcounts))]
results$dcounts2<-dcounts[match(results$ID2, names(dcounts))]
results$joint_counts<-results$dcounts1+results$dcounts2

windows()
plot(realSRI~VI, data=results, col=as.factor(affiliation975), bg=as.factor(status), pch=21)


table(results[results$dcounts1>=35 
              & results$dcounts2>=35,"status"], results[results$dcounts1>=35 
                                                        & results$dcounts2>=35,"affiliation975"])

###ok, try timestamp swapping method

#make in format of sim_surveys

sim_surveys_swap<-list()

for (i in 1:num_sim) {
  swapped_obs<-split(all_obs[,1:4], as.factor(all_obs$id)) 
  swapped_obs<-lapply(swapped_obs, function(x){ x$day<-sample(x$day) 
                                                return(x)})
  swapped_obs<-do.call("rbind", swapped_obs)
  sim_surveys1<-split(swapped_obs, f=as.factor(swapped_obs$day))
  sim_surveys_swap[[i]]<-sim_surveys1
}

kfinal_swap<-group_assign(data=sim_surveys_swap, id="id", xcoord ="x", ycoord="y", time="day", group_size = mean_group_size)

endtime<-Sys.time()

endtime-starttime #check run time #2.7 hours not in parallel

random_group_sizes<-lapply(kfinal_swap, function(x) mean(table(x$id)))
mean(unlist(random_group_sizes))


kfinal_swap<-lapply(kfinal_swap, function(x) {names(x)<-c("dolphin_id", "DayGroup", "Date", "observation_id");x})
kfinal_swap<-lapply(kfinal_swap, function(x) {x[,"observation_id"]<-as.numeric(as.factor(x[,"observation_id"]));x})
kfinal_swap<-lapply(kfinal_swap, function(x) {x[,"dolphin_id"]<-as.character(x[,"dolphin_id"]);x})

masked_randoms_swap<-lapply(kfinal_swap, function(x){
  simple_ratio(sightings=x,
               group_variable="observation_id", 
               dates="day", 
               IDs="id", 
               symmetric=FALSE)})

obs_network<-simple_ratio(sightings=all_obs,
                          group_variable="groupid", 
                          dates="day", 
                          IDs="id", 
                          symmetric=FALSE)


mn<-mat2dat(obs_network, "realSRI")

mrs<-lapply(masked_randoms_swap, mat2dat)

for(i in 1:length(mrs)) {names(mrs[[i]])[3]<-paste0("randSRI", i)}

mrs<-c(list(mn), mrs)

fdatas<-Reduce(function(x, y) merge(x, y, all.x = T), mrs) #about 2 minutes

#merge removes anything missing, individuals that didn't overlap have NaN and get removed in merge

fdatas$quantile975<-apply(fdatas, 1, function(x) quantile(as.numeric(x[4:(num_sim+3)]), 1, na.rm=TRUE))

fdatas$meanSRI<-apply(fdatas[,4:(num_sim+3)],1, mean, na.rm=TRUE)

fdatas$affiliation975<-ifelse(fdatas$realSRI>fdatas$quantile975, TRUE, FALSE)
fdatas<-fdatas[,grep("rand", names(fdatas), invert=TRUE)]

answer_key<-reduce_pairs(answer_key, "ID1", "ID2")

results_swap<-merge_pairs(answer_key, fdatas, "ID1", "ID2", all.x=TRUE, all.y=FALSE)


windows()
plot(realSRI~VI, data=results_swap, col=as.factor(affiliation975), bg=as.factor(status), pch=21)


table(results_swap$status, results_swap$affiliation975)

table(results$status, results$affiliation975)





