get_language_data <- function(file_path) {
  language = read.table(file_path, header = FALSE)
  colnames(language) = c("vertices","degree_2nd_moment", "mean_length")
  language = language[order(language$vertices), ]
  return(language)
}

fit_model <- function(data, a_initial, b_initial) {
  nonlinear_model = nls(degree_2nd_moment~a*vertices^b, data=data,start = list(a = a_initial, b = b_initial), trace = TRUE)
}

print_results <- function(model) {
  cat("Deviance: ", deviance(nonlinear_model), '\n')
  
  cat("AIC: ", AIC(nonlinear_model), '\n')
  
  cat("s: ", sqrt(deviance(nonlinear_model)/df.residual(nonlinear_model)), '\n')
  
  cat("Best parameters: ", coef(nonlinear_model), '\n')
}

inital_plots <- function(data) {
  cat(language, length(data$vertices), mean(data$vertices), sd(data$vertices), mean(data$degree_2nd_moment), sd(data$degree_2nd_moment), '\n')
  
  plot(data$vertices, data$degree_2nd_moment,xlab = "vertices", ylab = paste(language, 'degree second moment'))
  
  plot(log(data$vertices), log(data$degree_2nd_moment), xlab = "log(vertices)", ylab = paste(language, "log(degree second moment)"))
  
  mean_data = aggregate(data, list(data$vertices), mean)
  plot(mean_data$vertices, mean_data$degree_2nd_moment, xlab = "vertices", ylab = paste(language, "mean degree second moment"))
  plot(log(mean_data$vertices), log(mean_data$degree_2nd_moment), xlab = "vertices", ylab = paste(language, "mean degree second moment"))
  
  plot(data$vertices, data$degree_2nd_moment,xlab = "vertices", ylab = "degree 2nd moment")
  lines(mean_data$vertices,mean_data$degree_2nd_moment, col = "green")
  lines(data$vertices,(1 - 1/data$vertices)*(5 - 6/data$vertices), col = "red")
  lines(data$vertices,4-6/data$vertices, col = "blue")
  lines(data$vertices,data$vertices-1, col = "blue")
}
  

# MAIN
source = read.table("catalan_only.txt",  header = TRUE, as.is = c("language","file"))

for (x in 1:nrow(source)) {
  language <- source$language[x]
  data = get_language_data(source$file[x])

  linear_model = lm(log(mean_length)~log(vertices), data)
    
  a_initial = exp(coef(linear_model)[1])
  b_initial = coef(linear_model)[2]
  
  mean_data = aggregate(data, list(data$vertices), mean)
  
  nonlinear_model <- fit_model(mean_data, a_initial, b_initial)
  print_results(nonlinear_model)
}

