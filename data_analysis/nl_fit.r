library(minpack.lm)
args = commandArgs(trailingOnly=TRUE)
data <- read_csv(args[2],show_col_types = FALSE)

start_values <- c(1.0,1.0)
if (length(args) ==4){
    start_values <- c(a=args[3], b=args[4])
}

fit <- nls(y ~ a * exp(b * x),
           start = start_values,
           algorithm = "port",
           control = nls.control(maxiter = 1000))
print(coef(fm1DNase1))
print(confint(fm1DNase1))
