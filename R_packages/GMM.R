#https://tinyheero.github.io/2015/10/13/mixture-model.html
library(ggplot2)
library(dplyr)
library(mixtools)

options(scipen = 999)

ggplot(faithful, aes(x = waiting)) +
  geom_density()


plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

set.seed(1)
wait <- faithful$waiting
mixmdl <- normalmixEM(wait, k = 2)

data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "blue", lwd = 1.5) +
  ylab("Density")

mixmdl$mu
mixmdl$sigma
mixmdl$lambda
post.df <- as.data.frame(cbind(x = mixmdl$x, mixmdl$posterior))
head(post.df, 10)  # Retrieve first 10 rows
post.df %>%
  filter(x > 66, x < 68)

post.df %>%
  mutate(label = ifelse(comp.1 > 0.3, 1, 2)) %>% 
  ggplot(aes(x = factor(label))) +
  geom_bar() +
  xlab("Component") +
  ylab("Number of Data Points")

post.df %>%
  mutate(label = ifelse(comp.1 > 0.8, 1, 2)) %>%
  ggplot(aes(x = factor(label))) +
  geom_bar() +
  xlab("Component") +
  ylab("Number of Data Points")


