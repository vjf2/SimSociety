#Ben Bolker's Bounded Random Walk Code 

walk <- function(n.times=100,
                 xlim=c(-20,20),
                 ylim=c(-20,20),
                 start=c(0,0),
                 stepsize=c(1,1)) {
  ## blank plot of arena
  # plot(c(0,0),type="n",xlim=xlim,ylim=ylim,
       # xlab="Easting",ylab="Northing") 
  ## extract starting point
  x <- start[1]
  y <- start[2]
  xl <- list()
  yl <- list()
  # ## define potential step sizes
  # steps <- 1/c(1,2,4,8,12,16)
  # ## all possible positive or negative steps for N-S movement
  # steps.y <- c(steps,-steps,0)
  # ## bias E-W movement by leaving out some positive steps
  # steps.x <- c(steps,-steps,0)
  #   # c(steps[c(1,5,6)],-steps,0)
  # ## plot starting location
  # points(x,y,pch=16,col="red",cex=1)
  for (i in 1:n.times) {
    repeat {
      ## pick jump sizes
      xi <- rnorm(1, 0, abs(rnorm(1, 0.2, 0.5)))
      yi <- rnorm(1, 0, abs(rnorm(1, 0.2, 0.5)))
      ## new candidate locations
      newx <- x+xi
      newy <- y+yi
      ## IF new locations are within bounds, then
      ##    break out of the repeat{} loop (otherwise
      ##    try again)
      if (newx>xlim[1] && newx<xlim[2] &&
          newy>ylim[1] && newy<ylim[2]) break
    }
    # lines(c(x,newx),c(y,newy),col="blue") ## draw move
    ## set new location to candidate location
    x<- xl[[i]] <- newx
    y<- yl[[i]] <-newy
  }

  mvmt<-cbind(unlist(xl), unlist(yl))
  
  return(mvmt)
  
  }

set.seed(101)
test<-walk(100)

#starting coordinates

x<-x1<-tkcs[,1]
y<-y1<-tkcs[,2]

n=2

w2<-lapply(1:n, function(i){              
                        walk(n.times=10000,
                        xlim=c(x[i]-20,x[i]+20),
                        ylim=c(y[i]-20,y[i]+20),
                        start=c(x[i],y[i]),
                        stepsize=c(1,1))
          })

plot(w2[[1]]) #[[1]],w2[[2]][[2]])

points(x1[1], y1[1], col="blue", pch=16) #starting point

w3<-do.call("rbind", w2)

