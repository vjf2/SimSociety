save(all_mats, file="all_mats.RData")



rand_mats<-lapply(all_mats, function(x) { 
AIs<-mat2dat(x, "SRI")

pp2<-merge_pairs(pp, AIs, "ID1", "ID2", all.x=TRUE, all.y=FALSE)
rp2<-merge_pairs(rp, AIs, "ID1", "ID2", all.x=TRUE, all.y=FALSE)
ap2<-merge_pairs(ap, AIs, "ID1", "ID2", all.x=TRUE, all.y=FALSE)


pp2$status<-"preference"
rp2$status<-"random"
ap2$status<-"avoidance"

p2_rand<-rbind(pp2, rp2, ap2)
return(p2_rand)
})

rms<-do.call("cbind", lapply(rand_mats, function(x) x$SRI))

qx<-apply(rms, 1, function(x) quantile(x, probs = c(0.025, 0.975)))

rms<-cbind(t(qx), rms)          

all_res<-cbind(p2_rand[,1:2], rms)

all_res<-merge_pairs(all_res[,1:4], real_prefs, "ID1", "ID2", all.x=TRUE, all.y=FALSE)

all_res<-reduce_pairs(all_res, "ID1", "ID2")

windows()
plot(all_res[all_res$status=="preference",4], 
     all_res[all_res$status=="preference",5], 
     xlab="97.5 percentile", ylab="SRI",
     xlim=c(0, 0.5), ylim=c(0,0.5))
abline(0, 1, col="red")

sum(all_res[which(all_res$status=="avoidance"),3]<all_res[which(all_res$status=="avoidance"),5])

check_rand<-all_res

check_rand$`2.5%`<-sample(check_rand$`2.5%`)
check_rand$`97.5%`<-sample(check_rand$`97.5%`)

sum(check_rand[which(check_rand$status=="random"),4]>check_rand[which(check_rand$status=="random"),5])

#57% of prefs correct
#19% of randoms are false positive pref
#11% of randoms are false positive avoid
#42% of prefs are false negative
#8% of avoids are false negative 


