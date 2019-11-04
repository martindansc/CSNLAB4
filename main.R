get_language_data <- function(file_path) {
  language = read.table(file_path, header = FALSE)
  colnames(language) = c("vertices","degree_2nd_moment", "mean_length")
  language = language[order(language$vertices), ]
  return(language)
}

print_null_results <- function(data) {
  RSS <- sum((data$mean_length-(data$vertices+1)/3)^2)
  n <- length(data$vertices)
  p <- 0
  s <- sqrt(RSS/(n - p))
  
  AIC <- n*log(2*pi) + n*log(RSS/n) + n + 2*(p + 1)
  
  cat("Deviance: ", RSS, '\n')
  
  cat("AIC: ", AIC, '\n')
  
  cat("s: ", s, '\n')
}

print_results <- function(model) {
  cat("Deviance: ", deviance(model), '\n')
  
  cat("AIC: ", AIC(model), '\n')
  
  cat("s: ", sqrt(deviance(model)/df.residual(model)), '\n')
  
  cat("Best parameters: ", coef(model), '\n')
}

save_plot <- function(language, name, data, mean_data, model, xlab="vertices", ylab="degree second moment") {
  directory <- paste("plots/", language, sep="")
  if(!dir.exists(directory)) {
    dir.create(directory)
  }
  
  png(filename=paste(directory, "/", name, ".png", sep = ""))
  plot(log(data$vertices), log(data$degree_2nd_moment),xlab = xlab, ylab = ylab)
  lines(log(mean_data$vertices), log(fitted(model)), col = "green")
  dev.off()
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

initial_check <- function(data) {
  
  v <- length(data$vertices)
  secDeg <- mean(data$degree_2nd_moment)
  meanL <- mean(data$mean_length)
  
  cat("Checking metrics satisfaction\n")
  if ((4 - (6/v)) > secDeg) {
    cat("First check failed\n")
    return(FALSE)
  }
  if(secDeg > v) {
    cat("Second check failed\n")
    return(FALSE)
  }
  if(((v/(8*(v-1)))*secDeg+(1/2)) > meanL) {
    cat("Third check failed\n")
    return(FALSE)
  }
  if(meanL > (v-1)) {
    cat("Forth check failed\n")
    return(FALSE)
  }
  cat("All checks passed\n")
  return(TRUE)
  
}

fit_model_0 <- function(language, data, mean_data) {
  print_null_results(data)
  print_null_results(mean_data)
  
  plot(data$vertices, data$degree_2nd_moment,xlab = "vertices", ylab = "degree 2nd moment")
  lines(mean_data$vertices,mean_data$degree_2nd_moment, col = "green")
  lines(data$vertices,(1 - 1/data$vertices)*(5 - 6/data$vertices), col = "red")
  lines(data$vertices,4-6/data$vertices, col = "blue")
  lines(data$vertices,data$vertices-1, col = "blue")
}

fit_model_1 <- function(language, data, mean_data) {
  linear_model = lm(log(degree_2nd_moment)~log(vertices), data)
  b_initial = coef(linear_model)[2]
  
  nonlinear_model <- nls(degree_2nd_moment~(vertices/2)^b, data=mean_data, start = list(b = b_initial), trace = TRUE)
  print_results(nonlinear_model)
  
  save_plot(language, "model_1", data, mean_data, nonlinear_model)
}

fit_model_1_plus <- function(language, data, mean_data) {
  linear_model = lm(log(degree_2nd_moment)~log(vertices), data)
  b_initial = coef(linear_model)[2]
  
  nonlinear_model <- nls(degree_2nd_moment~(vertices/2)^b+d, data=mean_data, start = list(b = b_initial, d = 0), trace = TRUE)
  print_results(nonlinear_model)
  
  save_plot(language, "model_1_+", data, mean_data, nonlinear_model)
}

fit_model_2 <- function(language, data, mean_data) {
  linear_model = lm(log(degree_2nd_moment)~log(vertices), data)
  
  a_initial = exp(coef(linear_model)[1])
  b_initial = coef(linear_model)[2]
  
  nonlinear_model <- nls(degree_2nd_moment~a*vertices^b, data=mean_data,start = list(a = a_initial, b = b_initial), trace = TRUE)
  print_results(nonlinear_model)
  
  save_plot(language, "model_2", data, mean_data, nonlinear_model)
}

fit_model_3 <- function(language, data, mean_data) {
  linear_model = lm(log(degree_2nd_moment)~vertices, data)
  
  a_initial = exp(coef(linear_model)[1]) 
  c_initial = coef(linear_model)[2]
  
  nonlinear_model <- nls(degree_2nd_moment~a*exp(c*vertices), data=mean_data,start = list(a = a_initial, c = c_initial), trace = TRUE)
  print_results(nonlinear_model)
  
  save_plot(language, "model_3", data, mean_data, nonlinear_model)
}

fit_model_4 <- function(language, data, mean_data) {
  linear_model = lm(log(degree_2nd_moment)~log(vertices), data)
  
  a_initial = coef(linear_model)[2]
  
  nonlinear_model <- nls(degree_2nd_moment~a*log(vertices), data=mean_data,start = list(a = a_initial), trace = TRUE)
  print_results(nonlinear_model)
  
  save_plot(language, "model_4", data, mean_data, nonlinear_model)
}


# MAIN
source = read.table("catalan_only.txt",  header = TRUE, as.is = c("language","file"))

for (x in 1:nrow(source)) {
  language <- source$language[x]
  data = get_language_data(source$file[x])
  
  if(initial_check(data) == TRUE) {
    mean_data = aggregate(data, list(data$vertices), mean)
    
    fit_model_0(language, data, mean_data)
    fit_model_1(language, data, mean_data)
    fit_model_2(language, data, mean_data)
    fit_model_3(language, data, mean_data)
    fit_model_4(language, data, mean_data)
    #fit_model_1_plus(language, data, mean_data)
  }
}

