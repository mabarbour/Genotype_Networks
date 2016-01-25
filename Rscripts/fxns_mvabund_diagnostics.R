# Curve plots that are meant to be plotted on top of meanvar.plots in the 'mvabund' package or any other mean-variance plots. 
poisson_curve <- function(from, to, add = TRUE, color = "blue", line.width = 3, line.type = 1){
  curve(expr = 1*x, # poisson assumes that mean = variance
        from = from, to = to, add = add, col = color, lwd = line.width, lty = line.type)
}

quasipoisson_curve <- function(from, to, quasi.scalar, add = TRUE, color = "red", line.width = 3, line.type = 1){
  curve(expr = quasi.scalar*x, # quaispoisson assumes that variance scales linearly with mean according to a constant.
        from = from, to = to, add = add, col = color, lwd = line.width, lty = line.type)
}

neg.binomial_curve <- function(from, to, theta.negbin, add = TRUE, color = "black", line.width = 3, line.type = 1){
  curve(expr = x + theta.negbin*x^2, # negative binomial assumes this mean-variance relationship (quadratic)
        from = from, to = to, add = add, col = color, lwd = line.width, lty = line.type)
}