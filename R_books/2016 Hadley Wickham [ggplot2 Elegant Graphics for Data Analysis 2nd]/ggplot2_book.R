library(ggplot2)
library(dplyr)
# =======================================================================================
# 1&2. Getting Started with ggplot2
# =======================================================================================
ggplot2_built_data <- data(package = "ggplot2")
ggplot2_built_data$results[,"Item"]

ggplot(select(mpg,manufacturer,model), aes(manufacturer)) +
  geom_bar()

mpg
ggplot(mpg,aes(x = cty, y = hwy)) + 
  geom_point()

ggplot(mpg, aes(model, manufacturer)) + geom_point()

ggplot(mpg, aes(cty, hwy)) + geom_point()
ggplot(diamonds, aes(carat, price)) + geom_point()
ggplot(economics, aes(date, unemploy)) + geom_line()
ggplot(mpg, aes(cty)) + geom_histogram(binwidth = 1)

ggplot(mpg, aes(displ, cty, colour = class)) +
  geom_point()

ggplot(mpg, aes(displ, hwy)) + geom_point(aes(colour = "blue"))
ggplot(mpg, aes(displ, hwy)) + geom_point(colour = "blue")

ggplot(mpg, aes(displ, hwy)) +
  geom_point(shape = 3) +
  facet_wrap(~cyl)

ggplot(mpg, aes(displ, hwy)) +
  geom_point() +
  geom_smooth(method = "lm")

library(mgcv)
ggplot(mpg, aes(displ, hwy)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x))

ggplot(mpg, aes(drv, hwy)) +
  geom_point()

ggplot(mpg, aes(drv, hwy)) + geom_jitter()
ggplot(mpg, aes(drv, hwy)) + geom_boxplot(fill = "pink")
ggplot(mpg, aes(drv, hwy)) + geom_violin(fill = "pink")

ggplot(mpg, aes(hwy)) + geom_histogram()
ggplot(mpg, aes(hwy)) + geom_freqpoly()

ggplot(mpg, aes(displ, colour = drv)) +
  geom_freqpoly(binwidth = 0.5)

ggplot(mpg, aes(displ, fill = drv)) +
  geom_histogram(binwidth = 0.5) +
  facet_wrap(~drv, ncol = 1)

ggplot(economics, aes(date, unemploy / pop)) +
  geom_line()
ggplot(economics, aes(date, uempmed)) +
  geom_line()

year <- function(x) as.POSIXlt(x)$year + 1900
ggplot(economics, aes(unemploy / pop, uempmed)) +
  geom_path(colour = "grey50") +
  geom_point(aes(colour = year(date)))

ggplot(mpg, aes(class, hwy)) + geom_boxplot()

ggplot(mpg, aes(reorder(class, hwy), hwy)) +
  geom_boxplot()

ggplot(diamonds,aes(carat)) + 
  geom_histogram(binwidth = 0.05)

ggplot(mpg, aes(cty, hwy)) +
  geom_point(alpha = 1 / 3)

ggplot(mpg, aes(cty, hwy)) +
  geom_point(alpha = 1 / 3) +
  xlab("city driving (mpg)") +
  ylab("highway driving (mpg)")

# Remove the axis labels with NULL
ggplot(mpg, aes(cty, hwy)) +
  geom_point(alpha = 1 / 3) +
  xlab(NULL) +
  ylab(NULL)


ggplot(mpg, aes(drv, hwy)) +
  geom_jitter(width = 0.1)

ggplot(mpg, aes(drv, hwy)) +
  geom_jitter(width = 0.25) +
  xlim("f", "r") +
  ylim(20, 30)

#> Warning: Removed 138 rows containing missing values (geom_point).
# For continuous scales, use NA to set only one limit
ggplot(mpg, aes(drv, hwy)) +
  geom_jitter(width = 0.25, na.rm = TRUE) +
  ylim(NA, 30)

p <- ggplot(mpg, aes(displ, hwy, colour = factor(cyl))) +
  geom_point()
p

# Save png to disk
ggsave("plot.pdf", width = 5, height = 5)
summary(p)


# =======================================================================================
# 3. Toolbox
# =======================================================================================
df <- data.frame(
  x = c(3, 1, 5),
  y = c(2, 4, 6),
  label = c("a","b","c")
)
p <- ggplot(df, aes(x, y, label = label)) +
  labs(x = NULL, y = NULL) + # Hide axis label
  theme(plot.title = element_text(size = 12)) # Shrink plot title
p + geom_point() + ggtitle("point")
p + geom_text() + ggtitle("text")
p + geom_bar(stat = "identity") + ggtitle("bar")
p + geom_tile() + ggtitle("raster")

p + geom_line() + ggtitle("line")
p + geom_area() + ggtitle("area")
p + geom_path() + ggtitle("path")
p + geom_polygon() + ggtitle("polygon")

df <- data.frame(x = 1, y = 3:1, family = c("sans", "serif", "mono"))
ggplot(df, aes(x, y)) +
  geom_text(aes(label = family, family = family))

df <- data.frame(x = 1, y = 3:1, face = c("plain", "bold", "italic"))
ggplot(df, aes(x, y)) +
  geom_text(aes(label = face, fontface = face))

df <- data.frame(trt = c("a", "b", "c"), resp = c(1.2, 3.4, 2.5))
ggplot(df, aes(resp, trt)) +
  geom_point() +
  geom_text(aes(label = paste0("(", resp, ")")), nudge_y = -0.25) +
  xlim(1, 3.6)

ggplot(mpg, aes(displ, hwy)) +
  geom_text(aes(label = model)) +
  xlim(1, 8)

ggplot(mpg, aes(displ, hwy)) +
  geom_text(aes(label = model), check_overlap = TRUE) +
  xlim(1, 8)

label <- data.frame(
  waiting = c(55, 80),
  eruptions = c(2, 4.3),
  label = c("peak one", "peak two")
)
ggplot(faithfuld, aes(waiting, eruptions)) +
  geom_tile(aes(fill = density)) +
  geom_label(data = label, aes(label = label))

ggplot(mpg, aes(displ, hwy, colour = class)) +
  geom_point()
ggplot(mpg, aes(displ, hwy, colour = class)) +
  geom_point(show.legend = FALSE) +
  directlabels::geom_dl(aes(label = class), method = "smart.grid")

ggplot(economics, aes(date, unemploy)) +
  geom_line()

ggplot(economics) +
  geom_rect(
    aes(xmin = start, xmax = end, fill = party),
    ymin = -Inf, ymax = Inf, alpha = 0.2,
    data = presidential
  ) +
  geom_vline(
    aes(xintercept = as.numeric(start)),
    data = presidential,
    colour = "grey50", alpha = 0.5
  ) +
  geom_text(
    aes(x = start, y = 2500, label = name),
    data = presidential,
    size = 3, vjust = 0, hjust = 0, nudge_x = 50
  ) +
  geom_line(aes(date, unemploy)) +
  scale_fill_manual(values = c("blue", "red"))

yrng <- range(economics$unemploy)
xrng <- range(economics$date)
caption <- paste(strwrap("Unemployment rates in the US have
                           varied a lot over the years", 40), collapse = "\n")

ggplot(economics, aes(date, unemploy)) +
  geom_line() +
  geom_text(aes(x, y, label = caption),
  data = data.frame(x = xrng[1], y = yrng[2], caption = caption),
  hjust = 0, vjust = 1, size = 4
  )

ggplot(economics, aes(date, unemploy)) +
  geom_line() +
  annotate("text", x = xrng[1], y = yrng[2], label = caption,
           hjust = 0, vjust = 1, size = 4
  )

ggplot(diamonds, aes(log10(carat), log10(price))) +
  geom_bin2d() +
  facet_wrap(~cut, nrow = 1)


mod_coef <- coef(lm(log10(price) ~log10(carat), data = diamonds))
ggplot(diamonds, aes(log10(carat), log10(price))) +
  geom_bin2d() +
  geom_abline(intercept = mod_coef[1], slope = mod_coef[2],
              colour = "white", size = 1) +
  facet_wrap(~cut, nrow = 1)

data(Oxboys, package = "nlme")
head(Oxboys)

ggplot(Oxboys, aes(age, height, group = Subject)) +
  geom_point() +
  geom_line()

ggplot(Oxboys, aes(age, height, group = Subject)) +
  geom_line() +
  geom_smooth(method = "lm", se = FALSE)

ggplot(Oxboys, aes(age, height)) +
  geom_line(aes(group = Subject)) +
  geom_smooth(method = "lm", size = 2, se = FALSE)

ggplot(Oxboys, aes(Occasion, height)) +
  geom_boxplot()

ggplot(Oxboys, aes(Occasion, height)) +
  geom_boxplot() +
  geom_line(aes(group = Subject), colour = "#3366FF", alpha = 0.5)

df <- data.frame(x = 1:3, y = 1:3, colour = c(1,3,5))
ggplot(df, aes(x, y, colour = factor(colour))) +
  geom_line(aes(group = 1), size = 2) +
  geom_point(size = 5)
ggplot(df, aes(x, y, colour = colour)) +
  geom_line(aes(group = 1), size = 2) +
  geom_point(size = 5)

ggplot(mpg, aes(drv)) +
  geom_bar()
ggplot(mpg, aes(drv, fill = hwy, group = hwy)) +
  geom_bar()
library(dplyr)
mpg2 <- mpg %>% arrange(hwy) %>% mutate(id = seq_along(hwy))
ggplot(mpg2, aes(drv, fill = hwy, group = id)) +
  geom_bar()


ggplot(faithfuld, aes(eruptions, waiting)) +
  geom_contour(aes(z = density, colour = ..level..))

ggplot(faithfuld, aes(eruptions, waiting)) +
  geom_raster(aes(fill = density))

# Bubble plots work better with fewer observations
small <- faithfuld[seq(1, nrow(faithfuld), by = 10), ]
ggplot(small, aes(eruptions, waiting)) +
  geom_point(aes(size = density), alpha = 1/3) +
  scale_size_area()

library(dplyr)
mi_counties <- map_data("county", "michigan") %>%
  dplyr::select(lon = long, lat, group, id = subregion)


mi_cities <- maps::us.cities %>%
  tbl_df() %>%
  filter(country.etc =="MI") %>%
  dplyr::select(-country.etc, lon = long) %>%
  arrange(desc(pop))
mi_cities


if (file.exists("mi_raster.rds")) {
  mi_raster <- readRDS("mi_raster.rds")
} else {
  bbox <- c(
    min(mi_counties$lon), min(mi_counties$lat),
    max(mi_counties$lon), max(mi_counties$lat)
  )
  mi_raster <- ggmap::get_openstreetmap(bbox, scale = 8735660)
  saveRDS(mi_raster, "mi_raster.rds")
}
#(Finding the appropriate scale required a lot of manual tweaking.) You can then plot it with:
  ggmap::ggmap(mi_raster)
ggmap::ggmap(mi_raster) +
  geom_point(aes(size = pop), mi_cities, colour = "red") +
  scale_size_area()
#If you have raster data from the raster package, you can convert it to the form needed by ggplot2 with the following code:
  df <- as.data.frame(raster::rasterToPoints(x))
names(df) <- c("lon", "lat", "x")
ggplot(df, aes(lon, lat)) +
  geom_raster(aes(fill = x))

mi_census <- midwest %>%
  tbl_df() %>%
  filter(state =="MI") %>%
  mutate(county = tolower(county)) %>%
  dplyr::select(county, area, poptotal, percwhite, percblack)
mi_census

census_counties <- left_join(mi_census, mi_counties, by = c("county" ="id"))
census_counties

ggplot(census_counties, aes(lon, lat, group = county)) +
  geom_polygon(aes(fill = poptotal)) +
  coord_quickmap()
ggplot(census_counties, aes(lon, lat, group = county)) +
  geom_polygon(aes(fill = percwhite)) +
  coord_quickmap()


y <- c(18, 11, 16)
df <- data.frame(x = 1:3, y = y, se = c(1.2, 0.5, 1.0))
base <- ggplot(df, aes(x, y, ymin = y - se, ymax = y + se))
base + geom_crossbar()
base + geom_pointrange()
base + geom_smooth(stat = "identity")
base + geom_errorbar()
base + geom_linerange()
base + geom_ribbon()

# Unweighted
ggplot(midwest, aes(percwhite, percbelowpoverty)) +
  geom_point()
# Weight by population
ggplot(midwest, aes(percwhite, percbelowpoverty)) +
  geom_point(aes(size = poptotal / 1e6)) +
  scale_size_area("Population\n(millions)", breaks = c(0.5, 1, 2, 4))

# Unweighted
ggplot(midwest, aes(percwhite, percbelowpoverty)) +
  geom_point() +
  geom_smooth(method = lm, size = 1)
# Weighted by population
ggplot(midwest, aes(percwhite, percbelowpoverty)) +
  geom_point(aes(size = poptotal / 1e6)) +
  geom_smooth(aes(weight = poptotal), method = lm, size = 1) +
  scale_size_area(guide = "none")

ggplot(midwest, aes(percbelowpoverty)) +
  geom_histogram(binwidth = 1) +
  ylab("Counties")
ggplot(midwest, aes(percbelowpoverty)) +
  geom_histogram(aes(weight = poptotal), binwidth = 1) +
  ylab("Population (1000s)")

ggplot(diamonds, aes(depth)) +
  geom_histogram()
ggplot(diamonds, aes(depth)) +
  geom_histogram(binwidth = 0.1) +
  xlim(55, 70)

ggplot(diamonds, aes(depth)) +
  geom_freqpoly(aes(colour = cut), binwidth = 0.1, na.rm = TRUE) +
  xlim(58, 68) +
  theme(legend.position = "none")

ggplot(diamonds, aes(depth)) +
  geom_histogram(aes(fill = cut), binwidth = 0.1, position = "fill",na.rm = TRUE) +
  xlim(58, 68) +
  theme(legend.position = "none")

ggplot(diamonds, aes(depth)) +
  geom_density(na.rm = TRUE) +
  xlim(58, 68) +
  theme(legend.position = "none")

ggplot(diamonds, aes(depth, fill = cut, colour = cut)) +
  geom_density(alpha = 0.2, na.rm = TRUE) +
  xlim(58, 68) +
  theme(legend.position = "none")

df <- data.frame(x = rnorm(2000), y = rnorm(2000))
norm <- ggplot(df, aes(x, y)) + xlab(NULL) + ylab(NULL)
norm + geom_point()
norm + geom_point(shape = 1) # Hollow circles
norm + geom_point(shape = ".") # Pixel sized

norm + geom_point(alpha = 1 / 3)
norm + geom_point(alpha = 1 / 5)
norm + geom_point(alpha = 1 / 10)

norm + geom_bin2d()
norm + geom_bin2d(bins = 10)
norm + geom_hex()
norm + geom_hex(bins = 10)

ggplot(diamonds, aes(color)) +
  geom_bar()
ggplot(diamonds, aes(color, price)) +
  geom_bar(stat = "summary_bin", fun.y = mean)

ggplot(diamonds, aes(table, depth)) +
  geom_bin2d(binwidth = 1, na.rm = TRUE) +
  xlim(50, 70) +
  ylim(50, 70)

ggplot(diamonds, aes(table, depth, z = price)) +
  geom_raster(binwidth = 1, stat = "summary_2d", fun = mean,na.rm = TRUE) +
  xlim(50, 70) +
  ylim(50, 70)


# =======================================================================================
# 4. Mastering the Grammar
# =======================================================================================

ggplot(mpg, aes(displ, hwy, colour = factor(cyl))) +
  geom_point()

ggplot(mpg, aes(displ, hwy, colour = factor(cyl))) +
  geom_line() +
  theme(legend.position = "none")
ggplot(mpg, aes(displ, hwy, colour = factor(cyl))) +
  geom_bar(stat = "identity", position = "identity", fill = NA) +
  theme(legend.position = "none")

ggplot(mpg, aes(displ, hwy, colour = factor(cyl))) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(mpg, aes(displ, hwy)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~year)

# =======================================================================================
# 5. Build a Plot Layer by Layer
# =======================================================================================
p <- ggplot(mpg, aes(displ, hwy))
p
p + geom_point()

p + layer(
  mapping = NULL,
  data = NULL,
  geom = "point",
  stat = "identity",
  position = "identity"
)

mod <- loess(hwy ~ displ, data = mpg)
grid <- data_frame(displ = seq(min(mpg$displ), max(mpg$displ), length = 50))
grid$hwy <- predict(mod, newdata = grid)
grid

std_resid <- resid(mod) / mod$s
outlier <- filter(mpg, abs(std_resid) > 2)
outlier

ggplot(mpg, aes(displ, hwy)) +
  geom_point() +
  geom_line(data = grid, colour = "blue", size = 1.5) +
  geom_text(data = outlier, aes(label = model))

ggplot(mapping = aes(displ, hwy)) +
  geom_point(data = mpg) +
  geom_line(data = grid) +
  geom_text(data = outlier, aes(label = model))

library(dplyr)
class <- mpg %>%
  group_by(class) %>%
  summarise(n = n(), hwy = mean(hwy))

ggplot(mpg, aes(class, hwy)) +
  geom_jitter(width = 0.05, size = 2) +
  geom_point(aes(y = hwy), data = class, size = 4, color = "red") +
  geom_text(aes(y = 10, label = paste0("n = ", n)), data = class)

ggplot(mpg, aes(displ, hwy, colour = class)) +
  geom_point()
ggplot(mpg, aes(displ, hwy)) +
  geom_point(aes(colour = class))
ggplot(mpg, aes(displ)) +
  geom_point(aes(y = hwy, colour = class))
ggplot(mpg) +
  geom_point(aes(displ, hwy, colour = class))


ggplot(mpg, aes(displ, hwy, colour = class)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme(legend.position = "none")

ggplot(mpg, aes(displ, hwy)) +
  geom_point(aes(colour = class)) +
  geom_smooth(method = "lm", se = FALSE) +
  theme(legend.position = "none")

ggplot(mpg, aes(cty, hwy)) +
  geom_point(colour = "darkblue")
ggplot(mpg, aes(cty, hwy)) +
  geom_point(aes(colour = "darkblue"))

ggplot(mpg, aes(cty, hwy)) +
  geom_point(aes(colour = "darkblue")) +
  scale_colour_identity()

ggplot(mpg, aes(displ, hwy)) +
  geom_point() +
  geom_smooth(aes(colour = "loess"), method = "loess", se = FALSE) +
  geom_smooth(aes(colour = "lm"), method = "lm", se = FALSE) +
  labs(colour = "Method")

ggplot(mpg) +
  geom_point(aes(mpg$disp, mpg$hwy))
ggplot() +
  geom_point(mapping = aes(y = hwy, x = cty), data = mpg) +
  geom_smooth(data = mpg, mapping = aes(cty, hwy))
ggplot(diamonds, aes(carat, price)) +
  geom_point(aes(log(brainwt), log(bodywt)), data = msleep)

ggplot(mpg) +
  geom_point(aes(class, cty)) +
  geom_boxplot(aes(trans, hwy))


mod <- loess(hwy ~displ, data = mpg)
smoothed <- data.frame(displ = seq(1.6, 7, length = 50))
pred <- predict(mod, newdata = smoothed, se = TRUE)
smoothed$hwy <- pred$fit
smoothed$hwy_lwr <- pred$fit - 1.96 * pred$se.fit
smoothed$hwy_upr <- pred$fit + 1.96 * pred$se.fit

ggplot(mpg, aes(displ, hwy)) +
  geom_ribbon(data = smoothed, aes(ymin = hwy_lwr, ymax = hwy_upr, color = "red"), fill = "grey") +
  geom_point(position = "jitter") +
  geom_line(data = smoothed, aes(y = hwy), size = 1, color = "red")


mpg %>% 
  group_by(drv, trans) %>% 
  summarise(n = n()) %>% 
  mutate(drvtrans = paste0(drv, "-", trans)) %>% 
  ggplot() +
  geom_bar(aes(drvtrans, n), stat = "identity") +
  coord_flip()

dplot <- ggplot(diamonds, aes(color, fill = cut)) +
  xlab(NULL) + ylab(NULL) + theme(legend.position = "none")

dplot + geom_bar()
dplot + geom_bar(position = "fill")
dplot + geom_bar(position = "dodge")

dplot + geom_bar(position = "identity", alpha = 1 / 2, colour = "grey50")
ggplot(diamonds, aes(color, colour = cut)) +
  geom_line(aes(group = cut), stat = "count") +
  xlab(NULL) + ylab(NULL) +
  theme(legend.position = "none")

ggplot(mpg, aes(displ, hwy)) +
  geom_point(position = "jitter")
ggplot(mpg, aes(displ, hwy)) +
  geom_point(position = position_jitter(width = 0.05, height = 0.5))
ggplot(mpg, aes(displ, hwy)) +
  geom_jitter(width = 0.05, height = 0.5)


# =======================================================================================
# 6. Scales, Axes and Legends
# =======================================================================================

































