library(ggplot2)

Tw = 40
Tp = 25

x = 1:(Tw + Tp)
print(x)

uniform <- function() {
  y = rep(1, 65)
}

threetwoone <- function() {
  y = c(
    rep(3, 30),
    rep(2, 10),
    rep(1, 25)
  )
}

opt_af <- function(beta) {
  y = exp(beta[1] * ((Tw + Tp - 1):0) + beta[2] * c((Tw-1):0, rep(0,25)) + beta[3] * c((Tw-1):0, rep(0,Tp))^2)
}

df = data.frame(x, uniform(), threetwoone(), opt_af(c(0.0037, -0.0265, 0.0018)))
colnames(df) = c('x', 'y1', 'y2', 'y3')

ggplot(df, aes(x=x)) + 
  geom_line(aes(y=y1)) + 
  geom_line(aes(y=y2)) +
  geom_line(aes(y=y3)) +
  ylim(0,5) +
  xlim(0,70)
