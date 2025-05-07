pathpre = "/nfs/turbo/jiankanggroup/gxma/hejunhui_HCP/HCP data"
pathsuf = "_2bk-0bk.nii"
library(oro.nifti)
library(neurobase)
library(readxl)
library(FLGP)
library(kernlab)

if(!file.exists(file.path(pathpre, "tables", "brain"))) {dir.create(file.path(pathpre, "tables", "brain"))}

task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
idx = task_id
m_ratios = c(0.2, 0.3, 0.4)
num_exps = 10
a = 0.01
mset = array(NA, c(4, num_exps, length(m_ratios)))
mse = array(NA, c(4, num_exps, length(m_ratios)))
nll = array(NA, c(4, num_exps, length(m_ratios)))
time = array(NA, c(4, num_exps, length(m_ratios)))
s_ratio = 0.1; r = 3; K_ratio = 0.3
models = list(subsample="kmeans", kernel="lae", gl="cluster-normalized", root=TRUE)
index.subject = read.csv(file.path(pathpre, "index2bk-0bk.csv"))
subject = index.subject$subject[idx]
# read data
data = read.csv(file.path(pathpre, "brain", paste0(subject,"_2bk-0bk.csv")))
# Activation Voxel
I = dim(data)[1]

load("/home/gxma/FLAG/act_count.rdata")
act_idx = (sapply(c(1:I), function(i) {return(act_count[data[i,1], data[i,2], data[i,3]] >= 2)}) & (abs(data[,4]) > qnorm(0.75)))
act_idx[is.na(act_idx)] = FALSE

if (sum(act_idx) > 5000) {
  # Process data
  data_act = data[act_idx,]
  X = as.matrix(data_act[,1:3])
  X = apply(X,2,function(x) {return((x-min(x))/(max(x)-min(x)))}) # coordinate transform
  y = data_act[,4]
  n = nrow(X)
  # Run Regression
  set.seed(1234)
  s = round(n*s_ratio); K = round(s*K_ratio)
  for(i in 1:length(m_ratios)) {
    for(num_exp in 1:num_exps) {
      # regression
      m = round(n*m_ratios[i])
      train.index = sample.int(n, m)
      test.index = c(1:n)[-train.index]
      X.train = X[train.index,]; X.test = X[test.index,]
      y.train = y[train.index]; y.test = y[test.index]

      # radial basis kernel function
      score_var <- function(log_var) {
        gp = gausspr(x=X.train, y=y.train, scaled=FALSE, var=exp(log_var))
        # compute rbf kernel matrix
        K_noise = kernelMatrix(gp@kernelf, X.train) + diag(exp(log_var), m)
        R = chol(K_noise)
        alpha = forwardsolve(t(R), y.train)
        alpha = backsolve(R, alpha)
        logdet = 2 * sum(log(diag(R)))

        nmll_gp = 0.5 * (crossprod(y.train, alpha) + logdet)
        return(nmll_gp)
      }

      opt = optimize(score_var, interval = log(c(1e-2, 3)), tol = 1e-2)
      var <- exp(opt$minimum)

      t1 = Sys.time()
      # var = res_lkflag.hcp$pars[2] + 5e-3
      rbf.hcp = gausspr(x=X.train, y=y.train, scaled=FALSE, var=var, variance.model=TRUE)
      train_pred.hcp = predict(rbf.hcp, X.train)
      test_pred.hcp = predict(rbf.hcp, X.test)
      t2 = Sys.time()
      mset.rbf = mean((y.train-train_pred.hcp)^2)
      mset[1,num_exp,i] = mset.rbf
      idx.rbf = order((y.test-test_pred.hcp)^2)[round((n-m)*a):round((n-m)*(1-a))]
      mse.rbf.win = mean(((y.test-test_pred.hcp)^2)[idx.rbf])
      # mse.rbf.win = mean(sort((y.test-test_pred.hcp)^2)[round((n-m)*a):round((n-m)*(1-a))])
      mse[1,num_exp,i] = mse.rbf.win

      mean = cbind(test_pred.hcp[idx.rbf])
      # cov = cbind(predict(rbf.hcp, X.test, type="variance")[idx.rbf] + var)

      kernel = rbf.hcp@kernelf
      C11 = kernelMatrix(kernel, X.train) + diag(rep(var,m))
      C21 = kernelMatrix(kernel, X.test, X.train)
      C22 = rep(1+var,n-m)
      cov = cbind((C22 - rowSums((C21%*%solve(C11))*C21))[idx.rbf])

      nll.rbf = negative_log_likelihood(mean, cov, y.test[idx.rbf], "regression")
      nll[1,num_exp,i] = nll.rbf

      time[1,num_exp,i] = as.numeric(t2-t1, units="hours")*3600

      # Nystrom extension
      t1 = Sys.time()
      res_nystrom.hcp = fit_nystrom_regression_gp_rcpp(X.train, y.train, X.test, s, K, sigma = 5e-3)
      t2 = Sys.time()
      y_nystrom.hcp = res_nystrom.hcp$Y_pred
      train_pred_nystrom.hcp = y_nystrom.hcp$train
      test_pred_nystrom.hcp = y_nystrom.hcp$test
      mset.nystrom = mean((y.train-train_pred_nystrom.hcp)^2)
      mset[2,num_exp,i] = mset.nystrom
      idx.nystrom = order((y.test-test_pred_nystrom.hcp)^2)[round((n-m)*a):round((n-m)*(1-a))]
      mse.nystrom.win = mean(((y.test-test_pred_nystrom.hcp)^2)[idx.nystrom])
      # mse.nystrom.win = mean(sort((y.test-test_pred_nystrom.hcp)^2)[round((n-m)*a):round((n-m)*(1-a))])
      mse[2,num_exp,i] = mse.nystrom.win

      mean = cbind(res_nystrom.hcp$posterior$mean[idx.nystrom])
      cov = cbind(res_nystrom.hcp$posterior$cov[idx.nystrom])
      nll.nystrom = negative_log_likelihood(mean, cov, y.test[idx.nystrom], "regression")
      nll[2,num_exp,i] = nll.nystrom

      time[2,num_exp,i] = as.numeric(t2-t1, units="hours")*3600

      # SKFLGP
      models$subsample = "kmeans"; models$gl = "cluster-normalized"
      t1 = Sys.time()
      res_skflag.hcp = fit_se_regression_gp_rcpp(X.train, y.train, X.test, s, r, K, models = models, sigma = 5e-3)
      t2 = Sys.time()
      y_skflag.hcp = res_skflag.hcp$Y_pred
      train_pred_skflag.hcp = y_skflag.hcp$train
      test_pred_skflag.hcp = y_skflag.hcp$test
      mset.skflag = mean((y.train-train_pred_skflag.hcp)^2)
      mset[3,num_exp,i] = mset.skflag
      idx.skflag = order((y.test-test_pred_skflag.hcp)^2)[round((n-m)*a):round((n-m)*(1-a))]
      mse.skflag.win = mean(((y.test-test_pred_skflag.hcp)^2)[idx.skflag])
      # mse.lkflag.win = mean(sort((y.test-test_pred_lkflag.hcp)^2)[round((n-m)*a):round((n-m)*(1-a))])
      mse[3,num_exp,i] = mse.skflag.win

      mean = cbind(res_skflag.hcp$posterior$mean[idx.skflag])
      cov = cbind(res_skflag.hcp$posterior$cov[idx.skflag])
      nll.skflag = negative_log_likelihood(mean, cov, y.test[idx.skflag], "regression")
      nll[3,num_exp,i] = nll.skflag

      time[3,num_exp,i] = as.numeric(t2-t1, units="hours")*3600

      # LKFLGP
      models$subsample = "kmeans"; models$gl = "cluster-normalized"
      t1 = Sys.time()
      res_lkflag.hcp = fit_lae_regression_gp_rcpp(X.train, y.train, X.test, s, r, K, models = models, sigma = 5e-3)
      t2 = Sys.time()
      y_lkflag.hcp = res_lkflag.hcp$Y_pred
      train_pred_lkflag.hcp = y_lkflag.hcp$train
      test_pred_lkflag.hcp = y_lkflag.hcp$test
      mset.lkflag = mean((y.train-train_pred_lkflag.hcp)^2)
      mset[4,num_exp,i] = mset.lkflag
      idx.lkflag = order((y.test-test_pred_lkflag.hcp)^2)[round((n-m)*a):round((n-m)*(1-a))]
      mse.lkflag.win = mean(((y.test-test_pred_lkflag.hcp)^2)[idx.lkflag])
      # mse.lkflag.win = mean(sort((y.test-test_pred_lkflag.hcp)^2)[round((n-m)*a):round((n-m)*(1-a))])
      mse[4,num_exp,i] = mse.lkflag.win

      mean = cbind(res_lkflag.hcp$posterior$mean[idx.lkflag])
      cov = cbind(res_lkflag.hcp$posterior$cov[idx.lkflag])
      nll.lkflag = negative_log_likelihood(mean, cov, y.test[idx.lkflag], "regression")
      nll[4,num_exp,i] = nll.lkflag

      time[4,num_exp,i] = as.numeric(t2-t1, units="hours")*3600


      gc()
    }
  }
  save(list = c("mset", "mse", "nll", "time", "n"), file = paste("/home/gxma/FLAG/results0.01_common_5k/", subject, ".Rdata", sep = ""))
}
