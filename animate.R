library(tweenr)

#list of data frames

tiny_xm2<-xm2[xm2$day %in% 400:450,]

#jitter points
tiny_xm2$x<-jitter(tiny_xm2$x, amount=3)
tiny_xm2$y<-jitter(tiny_xm2$y, amount=3)

ldf<-split(tiny_xm2, as.factor(tiny_xm2$day))
# ldf<-c(ldf, ldf[1])

#original is sine-out

tween_data <- tween_states(ldf, tweenlength = 1,
                           statelength = 0, 
                           ease = "back-in-out", 
                           nframe = 500)

gcols=rainbow(n, alpha=0.4)

windows()
for(i in 1:max(tween_data$.frame)) {
  # png(file=paste0("pngs/testcase", sprintf("%04d", i), ".png"))
  plot(y~x, data=tween_data[tween_data$.frame==i,], 
     col=gcols, pch=15, cex=1,
     xlim=c(-80, 80), ylim=c(-60,60),
     xlab=NA, ylab=NA, xaxt="n", yaxt="n")
   Sys.sleep(.1)
  # dev.off()
}


# setwd("./pngs")
# shell("magick convert -delay 40 *.png example_1.gif")

#magick convert -delay 40 *.png example_4.gif


