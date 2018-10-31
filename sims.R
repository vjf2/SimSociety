#Generated simulated society

#library(devtools)
#install_github("vjf2/SocGen")

library(colorspace)
library(SocGen)
library(igraph)
library(adehabitatHR)
library(rgdal)
library(rgeos)
library(reticulate)

options(stringsAsFactors = FALSE)

use_python("C:/Python27") #add your python 2.7 directory

rmg_dir <- "../random-modular-network-generator-master/" #add directory for Sah 2014 random graph modules

#import python modules and script
nx <- import("networkx")
sg <- import_from_path("sequence_generator", path=rmg_dir)
rmg <- import_from_path("random_modular_generator_variable_modules", path=rmg_dir)
script <- readLines("generator_base_script.py")

#add graph properties
n <- 100 #number of nodes
dg <- 6 #average node degree
m <- 9 #number of modules
Q <- 0.75 #modularity

script<-sub("N= [0-9]+", paste0("N= ", n), script)
script<-sub("d= [0-9]+", paste0("d= ", dg), script)
script<-sub("m= [0-9]+", paste0("m= ", m), script)
script<-sub("Q= [0-9]+\\.[0-9]+", paste0("Q= ", Q), script)

writeLines(script, "generator_input_script.py")

#Note that the script uses a random number generator and not all
#starting configurations will converge on a final graph, so if you
#receive an error message simply re-run py_run_file until the graph 
#converges

py_run_file("generator_input_script.py")

g<-read_graph(paste0("output_files/random_graph_poisson_N", n, "_d", dg, ".graphml"), 
                     format="graphml")

set.seed(4)

layout <- layout.fruchterman.reingold(g)
tkcs <- layout * 2
el_pp <- get.edgelist(g)
el_pp[,1:2] <- apply(el_pp[,1:2], 2, function(x) paste0("id", x))
colnames(el_pp) <- c("ID1", "ID2")
el_pp <- as.data.frame(el_pp)
el_pp$status <- "preference"

#starting coordinates
x <- x1 <- tkcs[, 1]
y <- y1 <- tkcs[, 2]

cols = rainbow(n, alpha = 0.3)

maxDist <- 10
d <- 2000 #number of days (time-steps to run simulation)

#each individual takes a semi-random (Levy flight) bounded walk
walks <- lapply(1:n, function(i) {
  bounded_walk(n.times = d, maxDist, start = c(x[i], y[i]))
})

xm <- as.data.frame(do.call("rbind", walks))

xm$id <- paste0("id", rep(1:n, each = d))
xm$day <- rep(1:d, n)

names(xm)[1:2] <- c("x", "y")

xm <- xm[, c(3, 1:2, 4)] #rearrange columns

#limit avoidances to pairs with some minimum home range overlap
grid_buffer = 20
x <- seq(min(xm[, "x"]) - grid_buffer, max(xm[, "x"]) + grid_buffer, by = 1)
y <- seq(min(xm[, "y"]) - grid_buffer, max(xm[, "y"]) + grid_buffer, by = 1)
xy <- expand.grid(x = x, y = y)
coordinates(xy) <- ~ x + y
gridded(xy) <- TRUE

hrxydata <- SpatialPointsDataFrame(xm[, c("x", "y")], xm["id"])
uds_href <- kernelUD(hrxydata[, 1], grid = xy)
VI <- kerneloverlaphr(uds_href, method = "VI", percent = 90)
VI[upper.tri(VI)] <- NA
diag(VI) <- NA
mean(VI, na.rm = TRUE)
vi <- mat2dat(VI, "VI")

#select avoidances
dmt <- merge_pairs(vi, el_pp, "ID1", "ID2", all.x = TRUE, all.y = FALSE)
dmtr <- dmt[which(is.na(dmt$status)), ]
dmta <- dmtr[which(dmtr$VI > 0.10), ]

navoid <- nrow(el_pp)
wra <- sample(rownames(dmta), navoid)

#all pair classifications
pp <- el_pp[, 1:2]
ap <- dmta[rownames(dmta) %in% wra, 1:2]
rp <- dmtr[setdiff(rownames(dmtr), wra), 1:2]

#costs
pref_jw_cost <- 0.1
pref_mcost <- 0
rand_jw_cost <- 1
rand_mcost <- 1
avoid_jw_cost <- 10
avoid_mcost <- 2

pref_index <- matrix(c(pp[, 1], pp[, 2], pp[, 2], pp[, 1]), ncol = 2)
rand_index <- matrix(c(rp[, 1], rp[, 2], rp[, 2], rp[, 1]), ncol = 2)
avoid_index <- matrix(c(ap[, 1], ap[, 2], ap[, 2], ap[, 1]), ncol = 2)

#Start Simulation

xm2 <- xm #original backup
xm2$group <- NA


for (k in 1:(d - 1)) {
  dy1 <- xm2[xm2$day == k, ]
  dy1$group <- NA
  
  distmat <- as.matrix(dist(dy1[, 2:3]))
  dimnames(distmat) <- list(paste0("id", 1:n), paste0("id", 1:n))
  
  distmat[pref_index] <- distmat[pref_index] * pref_jw_cost + pref_mcost
  
  distmat[rand_index] <- distmat[rand_index] * rand_jw_cost + rand_mcost
  
  distmat[avoid_index] <- distmat[avoid_index] * avoid_jw_cost + avoid_mcost

  foragers <- paste0("id", sample(1:n, n * 0.15))
  distmat <- distmat[setdiff(rownames(distmat), foragers),
                     setdiff(colnames(distmat), foragers)]
  num_clust <- 20
  
  xc <- hclust(as.dist(distmat))
  xcc <- cutree(xc, k = num_clust)
  
  dy1$group <- xcc[match(dy1$id, names(xcc))]
  
  dy2 <- xm2[xm2$day == k + 1, ]
  
  gs <- unique(na.omit(dy1$group))
  
  
  for (i in gs) {
    nr <- nrow(dy1[which(dy1$group == i), ])
    dy1[which(dy1$group == i), "weight"] <-
      sample(c(2, rep(1, nr - 1)), nr)
    dy2[which(dy1$group == i), "x"] <-
      weighted.mean(dy2[which(dy1$group == i), "x"],
                    dy1[which(dy1$group ==
                                i), "weight"])
    dy2[which(dy1$group == i), "y"] <-
      weighted.mean(dy2[which(dy1$group == i), "y"],
                    dy1[which(dy1$group ==
                                i), "weight"])
  }
  
  comp <- cbind(dy1, dy2[, 2:3])
  
  comp$distance <-
    apply(comp, 1, function(x)
      dist(rbind(c(x[2], x[3]),
                 c(x[7], x[8]))))
  
  names(comp)[c(7, 8)] <- c("x1", "y1")
  
  comp$x1 <- ifelse(comp$distance > 20, (comp$x + comp$x1) / 2, comp$x1)
  comp$y1 <- ifelse(comp$distance > 20, (comp$y + comp$y1) / 2, comp$y1)
  
  xm2[xm2$day == k + 1, c("x", "y")] <- dy2[, c("x", "y")]
  xm2[xm2$day == k + 1, c("group")] <- dy1$group
} 

#Assign lone individuals unique group numbers
xm2$group[which(is.na(xm2$group))] <- n:(nrow(xm2[which(is.na(xm2$group)), ]) + (n - 1))
xm2$groupid <- paste0(xm2$day, "_", xm2$group)
mean(table(xm2$groupid))

#mean association rate
mat <-
  simple_ratio(
    sightings = xm2,
    group_variable = "groupid",
    dates = "day",
    IDs = "id",
    symmetric = FALSE
  )
mean(mat, na.rm = TRUE)

AIs <- mat2dat(mat, "SRI")

pp2 <- merge_pairs(pp, AIs, "ID1", "ID2", all.x = TRUE, all.y = FALSE)
rp2 <- merge_pairs(rp, AIs, "ID1", "ID2", all.x = TRUE, all.y = FALSE)
ap2 <- merge_pairs(ap, AIs, "ID1", "ID2", all.x = TRUE, all.y = FALSE)

mean(pp2$SRI)
mean(rp2$SRI)
mean(ap2$SRI)

#distribution of association indices and group sizes
dev.new()
par(mfrow = c(1, 2))
hist(table(xm2$groupid), nclass = 30)
hist(mat, nclass = 20)

#check percent groups size=1
sum(table(xm2$groupid) == 1) / length(unique(xm2$groupid))

#create UDs for each animal
hrxydata <- SpatialPointsDataFrame(xm2[, c("x", "y")], xm2["id"])
uds_href <- kernelUD(hrxydata[, 1], grid = xy)
VI <- kerneloverlaphr(uds_href, method = "VI", percent = 90)
mean(VI)
vi <- mat2dat(VI, "VI")

#combine category assignments with results

pp2$status <- "preference"
rp2$status <- "random"
ap2$status <- "avoidance"

p2 <- rbind(pp2, rp2, ap2)

results <-
  merge_pairs(p2, vi, "ID1", "ID2", all.x = TRUE, all.y = FALSE)

#plot results
dev.new()
# pdf(file="sim_res2000.pdf")
plot(
  SRI ~ VI,
  data = results,
  type = "n",
  xlab = "Home Range Overlap",
  ylab = "Association Index",
  yaxt = "n",
  cex.lab = 1.25
)
points(
  SRI ~ VI,
  data = results[which(results$status == "random"), ],
  col = adjustcolor("grey", alpha.f = 0.4),
  pch = 16
)
points(
  SRI ~ VI,
  data = results[which(results$status == "avoidance"), ],
  col = adjustcolor("red", alpha.f = 0.4),
  pch = 16
)
points(
  SRI ~ VI,
  data = results[which(results$status == "preference"), ],
  col = adjustcolor("green", alpha.f = 0.4),
  pch = 16
)
axis(2, las = 1)

abline(lm(SRI ~ VI + 0, data = results[which(results$status == "random"), ]), col =
         "grey", lwd = 1)
abline(lm(SRI ~ VI + 0, data = results[which(results$status == "avoidance"), ]), col =
         "red", lwd = 1)
abline(lm(SRI ~ VI + 0, data = results[which(results$status == "preference"), ]), col =
         "green", lwd = 1)

legend(
  "topleft",
  pch = 16,
  legend = c("Preference", "Random", "Avoidance"),
  col = c("green", "grey", "red"),
  cex = 1.25
)

# dev.off()

# write.csv(xm2, "output_files/sim_results_2000steps.csv", row.names = FALSE)
# write.csv(results, "output_files/sim_categories_2000steps.csv", row.names = FALSE)

