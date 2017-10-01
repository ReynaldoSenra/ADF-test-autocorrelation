MIT License

Copyright (c) 2017 Reynaldo Senra

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


library(urca)
# order accounts for upto wich order of autocorrelation you want to account in the ADF test
# q is the maximun lag that you are willing to reach when trying to correct the autocorrelation
# pvalu is the significance level at which you want to perform the Breush Godfrey test
adfnocorr<-function(x, type = "drift", lags = 7, selectlags = "AIC", order =5, order.by = NULL, q = 10, pvalu = 0.01) {
  breush<-c()
  x<-na.omit(x)
  adfresult<-ur.df(x, type, lags, selectlags)
  z<-diff(x)
  n <- length(z)
  if (type == "drift") {
    b<-(nrow(adfresult@testreg$coefficients)-1) # here I need to set one lag more than the one of the ADF 
    w <- embed(z, b)
    z.diff <- w[, 1]
    z.lag.1 <- x[b:n]
    trend <- b:n
    z.diff.lag = w[, 2:b]
    regre<- lm(z.diff ~ z.lag.1 + 1 + z.diff.lag)
    for (i in 1:order){
      godfrey<-bgtest(regre, order = i, order.by, type = "Chisq")
      breush[i]<-godfrey$p.value
    }
    if (any(breush <= pvalu)) {
      repeat {
        adfresult<-ur.df(x, type, lags = b, selectlags = "Fixed")
        b<-b+1
        w <- embed(z, b)
        z.diff <- w[, 1]
        z.lag.1 <- x[b:n]
        z.diff.lag = w[, 2:b]
        regre<- lm(z.diff ~ z.lag.1 + 1 + z.diff.lag)
        for (a in 1:order){
          godfrey<-bgtest(regre, order = a, type = "Chisq")
          breush[a]<-godfrey$p.value
        }
        if ((b == (q+1)) | (all(breush > pvalu))){
          if ((any(breush <= pvalu))) {
            warning("Significant autocorrelation problems the ADF test")
          }						
          break
        }
      }
    } 
  } else if (type == "none") {
    b<-nrow(adfresult@testreg$coefficients) # here I need to set one lag more than the one of the ADF 
    w <- embed(z, b)
    z.diff <- w[, 1]
    z.lag.1 <- x[b:n]
    trend <- b:n
    z.diff.lag = w[, 2:b]
    regre<- lm(z.diff ~ z.lag.1 - 1 + z.diff.lag)
    for (i in 1:order){
      godfrey<-bgtest(regre, order = i, order.by, type = "Chisq")
      breush[i]<-godfrey$p.value
    }
    if (any(breush <= pvalu)) {
      repeat {
        adfresult<-ur.df(x, type, lags = b, selectlags = "Fixed")
        b<-b+1
        w <- embed(z, b)
        z.diff <- w[, 1]
        z.lag.1 <- x[b:n]
        z.diff.lag = w[, 2:b]
        regre<- lm(z.diff ~ z.lag.1 - 1 + z.diff.lag)
        for (a in 1:order){
          godfrey<-bgtest(regre, order = a, type = "Chisq")
          breush[a]<-godfrey$p.value
        }
        if ((b == (q+1)) | (all(breush > pvalu))){
          if ((any(breush <= pvalu))) {
            warning("Significant autocorrelation problems the ADF test")
          }						
          break
        }
      }
    }
  } else {
    b<-(nrow(adfresult@testreg$coefficients)-2) # here I need to set one lag more than the one of the ADF 
    w <- embed(z, b)
    z.diff <- w[, 1]
    z.lag.1 <- x[b:n]
    trend <- b:n
    z.diff.lag = w[, 2:b]
    regre<- lm(z.diff ~ z.lag.1 + 1 + trend + z.diff.lag)
    for (i in 1:order){
      godfrey<-bgtest(regre, order = i, order.by, type = "Chisq")
      breush[i]<-godfrey$p.value
    }
    if (any(breush <= pvalu)) {
      repeat {
        adfresult<-ur.df(x, type, lags = b, selectlags = "Fixed")
        b<-b+1
        w <- embed(z, b)
        z.diff <- w[, 1]
        z.lag.1 <- x[b:n]
        trend <- b:n
        z.diff.lag = w[, 2:b]
        regre<- lm(z.diff ~ z.lag.1 + 1 + trend + z.diff.lag)
        for (a in 1:order){
          godfrey<-bgtest(regre, order = a, type = "Chisq")
          breush[a]<-godfrey$p.value
        }
        if ((b == (q+1)) | (all(breush > pvalu))){
          if ((any(breush <= pvalu))) {
            warning("Significant autocorrelation problems the ADF test")
          }						
          break
        }
      }
    }    
  }
  print(breush)
  summary(adfresult)
}
