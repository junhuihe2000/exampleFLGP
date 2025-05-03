# multinomial classification for RBF Gaussian process

library(kernlab)
library(FLGP)


train_multi_classification <- function(x, y) {
  J = max(y) + 1
  aug_y = multi_train_split(y)
  multi_model = list()
  for(j in 1:J) {
    model = gausspr(x=x, y=as.factor(aug_y[,j]), type="classification", scaled=FALSE)
    multi_model[[j]] = model
  }
  return(list("aug_y"=aug_y,
              "models"=multi_model,
              "x_train"=x))
}

predict_multi_classification <- function(x, multi_classifier) {
  aug_y = multi_classifier$aug_y
  models = multi_classifier$models
  x_train = multi_classifier$x_train
  J = length(models)
  n = nrow(x)
  probs = matrix(NA, n, J)
  mean = matrix(NA, nrow(x), J)
  cov = matrix(NA, nrow(x), J)
  for(j in 1:J) {
    model = models[[j]]
    output = predict(model, x, type="probabilities")
    probs[,j] = output[,2]

    kernel = model@kernelf
    C11 = kernelMatrix(kernel, x_train)
    C21 = kernelMatrix(kernel, x, x_train)
    C22 = rep(1, nrow(x))
    post = posterior_distribution_classification(C11, C21, cbind(C22), aug_y[,j])
    mean[,j] = post$mean
    cov[,j] = post$cov
  }
  y_hat = apply(probs, 1, function(x) {return(which.max(x)-1)})
  return(list("pred"=y_hat,
              "posterior"=list("mean"=mean, "cov"=cov)))
}
