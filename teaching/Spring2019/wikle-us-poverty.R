## Visualizing uncertainty in areal data with bivariate choropleth maps, map pixelation and glyph rotation, STAT

library(dplyr)
library(ggplot2)
library(grDevices)
library(grid)
library(gridExtra)
library(maptools)
library(plyr)
library(rgdal)
library(rgeos)


load('us_data.rda')
load('us_geo.rda')
us <- us_geo
us_pov <- us_data

qnt1_r <- quantile(us_pov$pov_rate, c(.3333))
qnt2_r <- quantile(us_pov$pov_rate, c(.6666))
qnt3_r <- quantile(us_pov$pov_rate, c(1))
qnt1_m <- quantile(us_pov$pov_moe, c(.3333))
qnt2_m <- quantile(us_pov$pov_moe, c(.6666))
qnt3_m <- quantile(us_pov$pov_moe, c(1))

f <- colorRamp(c("#CCCCFF", "#0000FF"))
g <- colorRamp(c("#FFFFCC", "#FFFF00"))
vec1 <- c(0, .5, 1, 0, .5, 1, 0, .5, 1)
vec2 <- c(0, 0, 0, .5, .5, .5, 1, 1, 1)
xy <- cbind(vec1,vec2)
xy <- as.data.frame(xy)

one <- f(xy$vec1)
two <- g(xy$vec2)
both <- as.data.frame(cbind(one, two))
colnames(both) <- c("r1", "g1", "b1", "r2", "g2", "b2")

both$r.ave <- round((both$r1 + both$r2)/2)
both$g.ave <- round((both$g1 + both$g2)/2)
both$b.ave <- round((both$b1 + both$b2)/2)
both$ave <- paste(both$r.ave, both$g.ave, both$b.ave)
colors <- both$ave
colors <- sapply(strsplit(colors, " "),
                 function(colors) rgb(colors[1], colors[2], colors[3], maxColorValue = 255))

xy <- expand.grid(x=seq(1, 3, by=1), y=seq(1, 3, by=1))
grid <- cbind(xy, colors)
grid$colors <- as.character(grid$colors)

us_pov$hex_code <-
    ifelse(us_pov$pov_rate<=qnt1_r & us_pov$pov_moe<=qnt1_m, colors[1],
           ifelse(us_pov$pov_rate<=qnt1_r & us_pov$pov_moe<=qnt2_m, colors[4],
                  ifelse(us_pov$pov_rate<=qnt1_r & us_pov$pov_moe<=qnt3_m, colors[7],
                         ifelse(us_pov$pov_rate<=qnt2_r & us_pov$pov_moe<=qnt1_m, colors[2],
                                ifelse(us_pov$pov_rate<=qnt2_r & us_pov$pov_moe<=qnt2_m, colors[5],
                                       ifelse(us_pov$pov_rate<=qnt2_r & us_pov$pov_moe<=qnt3_m, colors[8],
                                              ifelse(us_pov$pov_rate<=qnt3_r & us_pov$pov_moe<=qnt1_m, colors[3],
                                                     ifelse(us_pov$pov_rate<=qnt3_r & us_pov$pov_moe<=qnt2_m, colors[6],
                                                            ifelse(us_pov$pov_rate<=qnt3_r & us_pov$pov_moe<=qnt3_m, colors[9], "NA")))))))))

colnames(us_pov)[1] <- "GEO_ID"
us@data <- left_join(us@data, us_pov, by="GEO_ID")

us@data$id <- rownames(us@data)
points <- fortify(us, region="id")
final_pov <- join(points, us@data, by="id")

terciles <- c(min(us_pov$pov_rate), qnt1_r, qnt2_r, qnt3_r,
              min(us_pov$pov_moe), qnt1_m, qnt2_m, qnt3_m)
x <- c(.5, 1.5, 2.5, 3.5, 0, 0, 0, 0)
y <- c(0, 0, 0, 0, .5, 1.5, 2.5, 3.5)
angle <- c(-90, -90, -90, -90, 180, 180, 180, 180)
labels <- data.frame(x, y, terciles, angle)
labels$terciles <- as.character(labels$terciles)

color.grid <- ggplot() +
    geom_tile(data=grid, aes(x=x, y=y, fill=colors), colour="black", size=.1) +
    scale_fill_identity() +
    scale_y_continuous(trans="reverse") +
    coord_equal(ratio=1) +
    geom_text(aes(x=-1, y=2, label="Margin of error \n (terciles)", angle=-90), size=2.5) +
    geom_text(aes(x=2, y=-1, label="Poverty rate \n (terciles)", angle=180), size=2.5) +
    geom_text(data=labels, aes(x=x, y=y, label=terciles, angle=angle), size=2.5) +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),axis.title.y=element_blank(),
          panel.background=element_blank())

states <- map_data("state")

map <- ggplot() +
    geom_polygon(data=final_pov, aes(x=long, y=lat, group=group, fill=hex_code),
                 colour="black", size=.1) +
    scale_fill_identity() +
    labs(title="United States 2015 Poverty Map", subtitle = "Percentage of families whose income
was below the poverty level") +
    geom_path(data=states, aes(x=long, y=lat, group=group), size=.3) +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),axis.title.y=element_blank(),
          panel.background=element_blank(), plot.title = element_text(size = 20),
          plot.subtitle = element_text(size=15)) +
    coord_equal(ratio=1.3)

cg <- editGrob(ggplotGrob(color.grid), vp=viewport(angle=135))

lay <- rbind(c(1,1,1,1,1),
             c(1,1,1,1,1),
             c(1,1,1,1,2),
             c(1,1,1,1,1))
dev.new(dev.new())
plot <- grid.arrange(map, cg, layout_matrix=lay)
