#Create gif of simulation

library(tweenr)

xm2 <- read.csv("output_files/sim_results_2000steps.csv")

#list of data frames

tiny_xm2 <- xm2[xm2$day %in% 400:450, ]

#jitter points so individuals in same group are still visible
tiny_xm2$x <- jitter(tiny_xm2$x, amount = 0.5)
tiny_xm2$y <- jitter(tiny_xm2$y, amount = 0.5)

ldf <- split(tiny_xm2, as.factor(tiny_xm2$day))

tween_data <- tween_states(
  ldf,
  tweenlength = 1,
  statelength = 0,
  ease = "linear",
  nframe = 500
)

gcols = rainbow(length(unique(xm2$id)), alpha = 0.4)

windows()
for (i in 1:max(tween_data$.frame)) {
  # png(file=paste0("pngs/frames", sprintf("%04d", i), ".png"))
  plot(
    y ~ x,
    data = tween_data[tween_data$.frame == i, ],
    col = gcols,
    pch = 16,
    cex = 1,
    xlim = c(-20, 20),
    ylim = c(-20, 20),
    xlab = NA,
    ylab = NA,
    xaxt = "n",
    yaxt = "n"
  )
  Sys.sleep(.1)
  # dev.off()
}

#compiled into gif with imagemagick
#magick convert -delay 15 *.png simdays400_450.gif
