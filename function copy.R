## functions
library(tidyverse)
library(dplyr)
library(shiny)
library(devtools)
library(clinfun)
library(Rcpp)
library(RcppArmadillo)
library(shinycssloaders)

## posterior probability table
posterior_matrix <- read.csv("https://raw.githubusercontent.com/Wanni-Lei/data/main/posterior_matrix.csv")

## create scenarios
## use rcpp generate each scenarios
library(Rcpp)

# Define the Rcpp function in the same R script
cppFunction('
  DataFrame generateResultsRcpp(int nmax) {
    std::vector<int> n_values;
    
    std::vector<int> n1_values;
    std::vector<int> r_values;
    std::vector<int> r1_values;
  
    // Fix n
    for (int n = 6; n <= nmax; ++n) {
      // Fix n1
      for (int n1 = 3; n1 <= (n - 3); ++n1) {
       
        
        // Fix r1
        for (int r1 = (n1 - 1); r1 >= 0; --r1) {
          // Fix r
          for (int r = (r1 + (n-n1) -1); r >= (r1 + 1); --r) {
            n_values.push_back(n);
            n1_values.push_back(n1);
            r_values.push_back(r);
            r1_values.push_back(r1);
          }
        }
      }
    }
  
    // Create a DataFrame from the vectors
    DataFrame df = DataFrame::create(_["n"] = n_values, _["n1"] = n1_values, _["r"] = r_values, _["r1"] = r1_values);
  
    return df;
  }
')

# Call the Rcpp function, the maximum sample size for single arm is set to 100
results_rcpp <- generateResultsRcpp(100)

# Remove the first row (initialized with zeros)
results_rcpp <- results_rcpp[-1,]
rownames(results_rcpp) <- NULL

library(tidyverse)
results_rcpp <- as.data.frame(results_rcpp)
a <- c(6, 3, 4, 2) ## the scenarios table lacks this row
results_rcpp <- rbind(results_rcpp,a)
scenarios <- results_rcpp %>% ## this is the scenarios function
  arrange(n,n1,r1) 




optimal_minimax <- function(pa1,pb1, alpha_input, beta_input, pa0, pb0, delta, scenarios, posterior_matrix){
  
  simon <-ph2simon(pu=pa1, pa=pb1, alpha_input, beta_input)
  df <- as.data.frame(simon$out)
  df_optimal<- df[which.min(df$`EN(p0)`), ]
  n_simon <- as.numeric(df_optimal[4]) ## n of Simon's optimal design
  n1_simon <- max(10, n_simon) ## restrict upper bound of n
## chose maximum sample size for single arm
  if((pb1-pa1)>= 0.3){
  choice<- scenarios %>% 
    filter(0.8*n_simon <= n & n <= 1.4*n1_simon ## restrict n
           & n1/(n-n1) >= 0.55 & n1/(n-n1) <= 1 ## restrict n1/n2
           & pb0 >= (r1-2)/n1 & r1>= n1*(pa0-0.3) ## restrict r1
           & pa1 <= (r-r1)/(n-n1) & (r-r1)/(n-n1) <= (pb1+0.4) ) ## restrict r2
  } else{
    choice<- scenarios %>% 
      filter(0.9*n_simon <= n & n <= 1.2*n1_simon ## restrict n
             & n1/(n-n1) >= 0.65 & n1/(n-n1) <= 1.5 ## restrict n1/n2
             & pb0 >= (r1-2)/n1 & r1>= n1*(pa0-0.1) ## restrict r1
             & pa1 <= (r-r1)/(n-n1) & (r-r1)/(n-n1) <= (pb1+0.1)
             )
  }
  
  
  
  
 
  ## apply the simonBayesian function to the dataframe
  results <- choice %>% 
    mutate(EN_p0 = apply(choice, 1, function(x) simon_bayesian(x[1], x[2], x[3], x[4],
                                                               pa1=pa1, pb1=pb1, pa0=pa0, pb0=pb0, delta=delta, posterior_matrix=posterior_matrix )$ess),
           power = apply(choice, 1, function(x) simon_bayesian(x[1], x[2], x[3], x[4],
                                                               pa1=pa1, pb1=pb1, pa0=pa0, pb0=pb0, delta=delta,posterior_matrix=posterior_matrix )$power),
           alpha = apply(choice, 1, function(x) simon_bayesian(x[1], x[2], x[3], x[4],
                                                               pa1=pa1, pb1=pb1, pa0=pa0, pb0=pb0,delta=delta, posterior_matrix=posterior_matrix )$alpha))
  
  
  valid_scenario <- results %>% 
    filter(EN_p0 <= 2*n) %>% 
    filter((1 - power) <= beta_input & alpha <= alpha_input)
  
  if(nrow(valid_scenario)==0){
    stop("No feasible solution found") 
  }
  
   
  optimal <- valid_scenario %>%
    filter(EN_p0 == min(EN_p0))%>% ## filter minimum EN_p0
    ## filter power closest to power_input
    filter(abs(power - (1-beta_input)) == min(abs(power - (1-beta_input))))%>%
    filter(abs(alpha - alpha_input) == min(abs(alpha - alpha_input))) ## alpha power closest to alpha_input
  
  
  
  minimax <- valid_scenario %>%
    filter(n == min(n))%>% ## filter minimum EN_p0
    ## filter power closest to power_input
    filter(abs(power - (1-beta_input)) == min(abs(power - (1-beta_input))))%>%
    filter(abs(alpha - alpha_input) == min(abs(alpha - alpha_input)))%>%
    ## alpha power closest to alpha_input
    filter(EN_p0==min(EN_p0)) ## filter minimum expected sample size
  
  
  
  
  df <- as.data.frame(rbind(optimal,minimax))
  n_optimal <- nrow(optimal)
  n_minimax <- nrow(minimax)
  
  # Generate unique row names
  optimal_names <- paste("optimal", seq_len(n_optimal), sep = "_")
  minimax_names <- paste("minimax", seq_len(n_minimax), sep = "_")
  
  row.names(df) <- c(optimal_names, minimax_names)
  return(df)
  
}



simon_bayesian <- function(n, n1, r, r1, pa1, pb1, pa0, pb0, delta, posterior_matrix){
 
  
  if(pa1>pb1){
    stop("Error: pa1 must be less than pb1.")
  }
  
  n2 <- n-n1
  x12 <- data.frame(x1=rep(0:n1, each=n2+1), 
                    x2=rep(0:n2, n1+1)) 
  x12 <- x12 %>% filter(!(x1<=r1&x2>0)) %>%
    mutate(result=case_when(x1<=r1 ~ "F1", 
                            x1>r1&(x1+x2)<=r ~ "F2", 
                            x1>r1 & (x1+x2)>r ~ "P"),
           pa1=ifelse(x1<=r1, dbinom(x1, n1, p=pa1), 
                      dbinom(x1,n1, p=pa1)*dbinom(x2,n2, p=pa1)), 
           pb1=ifelse(x1<=r1, dbinom(x1, n1, p=pb1),
                      dbinom(x1,n1, p=pb1)*dbinom(x2,n2, p=pb1)),
           pa0=ifelse(x1<=r1, dbinom(x1, n1, p=pa0),
                      dbinom(x1,n1, p=pa0)*dbinom(x2,n2, p=pa0)),
           pb0=ifelse(x1<=r1, dbinom(x1, n1, p=pb0),
                      dbinom(x1,n1, p=pb0)*dbinom(x2,n2, p=pb0))
    )
  
  x12.smry <- x12[,-(1:2)] %>%
    group_by(result) %>%
    summarize(p.pa1=sum(pa1),
              p.pb1=sum(pb1),
              p.pa0=sum(pa0),
              p.pb0=sum(pb0)) %>% data.frame()
  
  row.names(x12.smry) <-x12.smry$result
  
  
  table1 <- outer(x12.smry$p.pa1,x12.smry$p.pb1) 
  
  table0 <- outer(x12.smry$p.pa0,x12.smry$p.pb0)
  
  
  ## when both 2 doses pase stage 1 
  nr.x12 <- dim(x12)[1] 
  ab <- data.frame(a1=rep(x12$x1,nr.x12),
                   a2=rep(x12$x2,nr.x12),
                   b1=rep(x12$x1,each=nr.x12),
                   b2=rep(x12$x2,each=nr.x12))
  
  
  ab.pf <- ab %>% 
    mutate(a.out=case_when(a1<=r1 ~ "aF1",
                           a1>r1&(a1+a2)<=r ~ "aF2",
                           a1>r1 & (a1+a2)>r ~ "aP"),
           b.out=case_when(b1<=r1 ~ "bF1",
                           b1>r1&(b1+b2)<=r ~ "bF2",
                           b1>r1 & (b1+b2)>r ~ "bP")) %>% 
    mutate(bH.aPbP=ifelse(a.out=="aP"& b.out=="bP",1,0))
  
  ## filter only (a1+a2)<(b1+b2)
  bW.aPbP <- ab.pf %>% 
    filter(bH.aPbP==1)  %>% 
    filter((a1+a2)<(b1+b2))
  
  ## Indicator matrix: 
  indicator <- data.frame(a1 = bW.aPbP$a1,
                          a2 = bW.aPbP$a2,
                          b1 = bW.aPbP$b1,
                          b2 = bW.aPbP$b2,
                          a_total = bW.aPbP$a1+bW.aPbP$a2,
                          b_total = bW.aPbP$b1+bW.aPbP$b2,
                          posterior = 0)
  
  ## Sort the data frame by ascending order of a_total and b_total, 
  ## create groups for unique a_total, b_total
  indicator <- indicator %>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total) %>%
    mutate(sample_size = n)  %>%
    mutate(group = cur_group_id()) 
  
  
  ## load posterior matrix
  ## changed simulation in p.rbH to 10000 times
  posterior_matrix_filter <- posterior_matrix %>%
    filter(sample_size==n &a_total>r & b_total>r &a_total<b_total)%>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total)
  
  ## add the "posterior" values from the "posterior_matrix_filter" to 
  ## the "indicator" based on matching values of
  ## "a_total", "b_total", and "sample_size"
  cppFunction('
  DataFrame add_posterior_values(DataFrame indicator, DataFrame posterior_matrix_filter, int n) {
  int nrows = indicator.nrows();
  NumericVector group = indicator["group"];
  NumericVector a_total = indicator["a_total"];
  NumericVector b_total = indicator["b_total"];
  NumericVector posterior = indicator["posterior"];
  
  int num_groups = posterior_matrix_filter.nrows();
  NumericVector filter_a_total = posterior_matrix_filter["a_total"];
  NumericVector filter_b_total = posterior_matrix_filter["b_total"];
  NumericVector filter_sample_size = posterior_matrix_filter["sample_size"];
  NumericVector filter_posterior = posterior_matrix_filter["posterior"];
  
  for (int i = 0; i < num_groups; ++i) {
    double cur_posterior = filter_posterior[i];
    for (int j = 0; j < nrows; ++j) {
      if (group[j] == (i + 1) && a_total[j] == filter_a_total[i] && b_total[j] == filter_b_total[i] && n == filter_sample_size[i]) {
        posterior[j] = cur_posterior;
      }
    }
  }
  
  return DataFrame::create(_["a1"] = indicator["a1"],
                           _["a2"] = indicator["a2"],
                           _["b1"] = indicator["b1"],
                           _["b2"] = indicator["b2"],
                           _["a_total"] = indicator["a_total"],
                           _["b_total"] = indicator["b_total"],
                           _["posterior"] = posterior,    
                           _["group"] = indicator["group"]);
}

')
  indicator <- indicator %>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total) %>%
    mutate(sample_size = n)  %>%
    mutate(group = cur_group_id()) 
  
  posterior_matrix_filter <- posterior_matrix %>%
    filter(sample_size==n &a_total>r & b_total>r &a_total<b_total)%>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total)
  
  indicator <- add_posterior_values(indicator, posterior_matrix_filter, n)
  
  
  ## 
  indicator_p <- indicator %>% 
    mutate(p.ra1rb1=dbinom(a1, n1, p=pa1)*dbinom(a2, n2, p=pa1)*
             dbinom(b1, n1, p=pb1)*dbinom(b2, n2, p=pb1),
           p.ra0rb0=dbinom(a1, n1, p=pa0)*dbinom(a2, n2, p=pa0)*
             dbinom(b1, n1, p=pb0)*dbinom(b2, n2, p=pb0))
  
  ##   
  bW.aPbP.smry <- apply(indicator_p[indicator_p$posterior>delta,
                                    c("p.ra1rb1","p.ra0rb0")],2,sum)
  
  
  ## expected total sample size under p0
  ess <- n1*2*table0[1,1]+              
    (n1*2+n2*2)*(table0[2,2] +table0[2,3]+table0[3,2] +table0[3,3])+    
    (n1*2+n2)*(table0[1,2] +table0[1,3]+table0[2,1] +table0[3,1])  
  
  
  
  ## power and type I error of pick the true winner
  bpower_1 <- bW.aPbP.smry["p.ra1rb1"]
  bpower_0 <- bW.aPbP.smry["p.ra0rb0"]
  
  ## when response rate of 2 doses are not equal
  power <-  table1[1,3]+table1[2,3]+bpower_1 ## pa_1 pb_1 
  alpha <-  table0[1,3]+table0[2,3]+bpower_0 ## pa_0 pb_0 
  
  #my_list <- round(list(n = n,n1 = n1, r = r, r1 = r1, ess=ess, power=power,alpha=alpha),3)
  vector <- round(c(ess, power, alpha), 3)
  vector_list <- as.list(vector)
  names(vector_list) <- c("ess", "power", "alpha")
  #vector_list
  return(vector_list)
}











## operating characteristics function
oc <- function(n, n1, r, r1, pa1, pb1, pa0, pb0, delta, posterior_matrix){
  
  if(pa1>pb1){
    stop("Error: pa1 must be less than pb1.")
  }
  
  n2 <- n-n1
  x12 <- data.frame(x1=rep(0:n1, each=n2+1), 
                    x2=rep(0:n2, n1+1)) 
  x12 <- x12 %>% filter(!(x1<=r1&x2>0)) %>%
    mutate(result=case_when(x1<=r1 ~ "F1", 
                            x1>r1&(x1+x2)<=r ~ "F2", 
                            x1>r1 & (x1+x2)>r ~ "P"),
           pa1=ifelse(x1<=r1, dbinom(x1, n1, p=pa1), 
                      dbinom(x1,n1, p=pa1)*dbinom(x2,n2, p=pa1)), 
           pb1=ifelse(x1<=r1, dbinom(x1, n1, p=pb1),
                      dbinom(x1,n1, p=pb1)*dbinom(x2,n2, p=pb1)),
           pa0=ifelse(x1<=r1, dbinom(x1, n1, p=pa0),
                      dbinom(x1,n1, p=pa0)*dbinom(x2,n2, p=pa0)),
           pb0=ifelse(x1<=r1, dbinom(x1, n1, p=pb0),
                      dbinom(x1,n1, p=pb0)*dbinom(x2,n2, p=pb0))
    )
  
  x12.smry <- x12[,-(1:2)] %>%
    group_by(result) %>%
    summarize(p.pa1=sum(pa1),
              p.pb1=sum(pb1),
              p.pa0=sum(pa0),
              p.pb0=sum(pb0)) %>% data.frame()
  
  row.names(x12.smry) <-x12.smry$result
  
  
  table1 <- outer(x12.smry$p.pa1,x12.smry$p.pb1) 
  
  table0 <- outer(x12.smry$p.pa0,x12.smry$p.pb0)
  
  
  ## when both 2 doses pase stage 1 
  nr.x12 <- dim(x12)[1] 
  ab <- data.frame(a1=rep(x12$x1,nr.x12),
                   a2=rep(x12$x2,nr.x12),
                   b1=rep(x12$x1,each=nr.x12),
                   b2=rep(x12$x2,each=nr.x12))
  
  
  ab.pf <- ab %>% 
    mutate(a.out=case_when(a1<=r1 ~ "aF1",
                           a1>r1&(a1+a2)<=r ~ "aF2",
                           a1>r1 & (a1+a2)>r ~ "aP"),
           b.out=case_when(b1<=r1 ~ "bF1",
                           b1>r1&(b1+b2)<=r ~ "bF2",
                           b1>r1 & (b1+b2)>r ~ "bP")) %>% 
    mutate(bH.aPbP=ifelse(a.out=="aP"& b.out=="bP",1,0))
  
  ## filter only (a1+a2)<(b1+b2)
  bW.aPbP <- ab.pf %>% 
    filter(bH.aPbP==1)  %>% 
    filter((a1+a2)<(b1+b2))
  
  ## Indicator matrix: 
  indicator <- data.frame(a1 = bW.aPbP$a1,
                          a2 = bW.aPbP$a2,
                          b1 = bW.aPbP$b1,
                          b2 = bW.aPbP$b2,
                          a_total = bW.aPbP$a1+bW.aPbP$a2,
                          b_total = bW.aPbP$b1+bW.aPbP$b2,
                          posterior = 0)
  
  ## Sort the data frame by ascending order of a_total and b_total, 
  ## create groups for unique a_total, b_total
  indicator <- indicator %>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total) %>%
    mutate(sample_size = n)  %>%
    mutate(group = cur_group_id()) 
  
  
  ## load posterior matrix
  ## changed simulation in p.rbH to 10000 times
  posterior_matrix_filter <- posterior_matrix %>%
    filter(sample_size==n &a_total>r & b_total>r &a_total<b_total)%>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total)
  
  ## add the "posterior" values from the "posterior_matrix_filter" to 
  ## the "indicator" based on matching values of
  ##"a_total", "b_total", and "sample_size"
  cppFunction('
  DataFrame add_posterior_values(DataFrame indicator, DataFrame posterior_matrix_filter, int n) {
  int nrows = indicator.nrows();
  NumericVector group = indicator["group"];
  NumericVector a_total = indicator["a_total"];
  NumericVector b_total = indicator["b_total"];
  NumericVector posterior = indicator["posterior"];
  
  int num_groups = posterior_matrix_filter.nrows();
  NumericVector filter_a_total = posterior_matrix_filter["a_total"];
  NumericVector filter_b_total = posterior_matrix_filter["b_total"];
  NumericVector filter_sample_size = posterior_matrix_filter["sample_size"];
  NumericVector filter_posterior = posterior_matrix_filter["posterior"];
  
  for (int i = 0; i < num_groups; ++i) {
    double cur_posterior = filter_posterior[i];
    for (int j = 0; j < nrows; ++j) {
      if (group[j] == (i + 1) && a_total[j] == filter_a_total[i] && b_total[j] == filter_b_total[i] && n == filter_sample_size[i]) {
        posterior[j] = cur_posterior;
      }
    }
  }
  
  return DataFrame::create(_["a1"] = indicator["a1"],
                           _["a2"] = indicator["a2"],
                           _["b1"] = indicator["b1"],
                           _["b2"] = indicator["b2"],
                           _["a_total"] = indicator["a_total"],
                           _["b_total"] = indicator["b_total"],
                           _["posterior"] = posterior,    
                           _["group"] = indicator["group"]);
}

')
  indicator <- indicator %>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total) %>%
    mutate(sample_size = n)  %>%
    mutate(group = cur_group_id()) 
  
  posterior_matrix_filter <- posterior_matrix %>%
    filter(sample_size==n &a_total>r & b_total>r &a_total<b_total)%>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total)
  
  indicator <- add_posterior_values(indicator, posterior_matrix_filter, n)
  
  
   
  indicator_p <- indicator %>% 
    mutate(p.ra1rb1=dbinom(a1, n1, p=pa1)*dbinom(a2, n2, p=pa1)*
             dbinom(b1, n1, p=pb1)*dbinom(b2, n2, p=pb1),
           p.ra0rb0=dbinom(a1, n1, p=pa0)*dbinom(a2, n2, p=pa0)*
             dbinom(b1, n1, p=pb0)*dbinom(b2, n2, p=pb0))
  
     
  bW.aPbP.smry <- apply(indicator_p[indicator_p$posterior> delta,
                                    c("p.ra1rb1","p.ra0rb0")],2,sum)
  
  
  ## expected total sample size under p0
  ess <- n1*2*table0[1,1]+              
    (n1*2+n2*2)*(table0[2,2] +table0[2,3]+table0[3,2] +table0[3,3])+    
    (n1*2+n2)*(table0[1,2] +table0[1,3]+table0[2,1] +table0[3,1])  
  
  
  
  ## power and type I error of pick the true winner
  bpower_1 <- bW.aPbP.smry["p.ra1rb1"]
  bpower_0 <- bW.aPbP.smry["p.ra0rb0"]
  
  ## when response rate of 2 doses are not equal
  power <-  table1[1,3]+table1[2,3]+bpower_1 ## pa_1 pb_1 
  alpha <-  table0[1,3]+table0[2,3]+bpower_0 ## pa_0 pb_0 
  
  A.ES.H0 <- sum(table0[1,1:3])
  B.ES.H0 <- colSums(table0)[1]
  AB.ES.H0 <- A.ES.H0 + B.ES.H0 - table0[1,1]
  A.P1.H1 <- 1- sum(table1[1,1:3])
  B.P1.H1 <- 1- colSums(table1)[1]
  AB.P1.H1 <- table1[2,2]+table1[2,3]+table1[3,2]+table1[3,3]
  
  
  vector <- list("n" = n,"n1" = n1, "r" = r, "r1" = r1, "ess"=ess, "power"=power,"alpha"=alpha, "A.ES.H0"=A.ES.H0, "B.ES.H0"=B.ES.H0,"AB.ES.H0"=AB.ES.H0, "A.P1.H1"=A.P1.H1,"B.P1.H1"=B.P1.H1, "AB.P1.H1"=AB.P1.H1)
  vector <- as.data.frame(vector)
  colnames(vector) <- c("n","n1","r","r1","ess", "power", "alpha","A.ES.H0", "B.ES.H0","AB.ES.H0", "A.P1.H1","B.P1.H1", "AB.P1.H1")
  rownames(vector) <- "Results"
  
  return(round(vector,3))
}


## return table1
octable1 <- function(n, n1, r, r1, pa1, pb1, pa0, pb0, delta, posterior_matrix){
  
  if(pa1>pb1){
    stop("Error: pa1 must be less than pb1.")
  }
  
  n2 <- n-n1
  x12 <- data.frame(x1=rep(0:n1, each=n2+1), 
                    x2=rep(0:n2, n1+1)) 
  x12 <- x12 %>% filter(!(x1<=r1&x2>0)) %>%
    mutate(result=case_when(x1<=r1 ~ "F1", 
                            x1>r1&(x1+x2)<=r ~ "F2", 
                            x1>r1 & (x1+x2)>r ~ "P"),
           pa1=ifelse(x1<=r1, dbinom(x1, n1, p=pa1), 
                      dbinom(x1,n1, p=pa1)*dbinom(x2,n2, p=pa1)), 
           pb1=ifelse(x1<=r1, dbinom(x1, n1, p=pb1),
                      dbinom(x1,n1, p=pb1)*dbinom(x2,n2, p=pb1)),
           pa0=ifelse(x1<=r1, dbinom(x1, n1, p=pa0),
                      dbinom(x1,n1, p=pa0)*dbinom(x2,n2, p=pa0)),
           pb0=ifelse(x1<=r1, dbinom(x1, n1, p=pb0),
                      dbinom(x1,n1, p=pb0)*dbinom(x2,n2, p=pb0))
    )
  
  x12.smry <- x12[,-(1:2)] %>%
    group_by(result) %>%
    summarize(p.pa1=sum(pa1),
              p.pb1=sum(pb1),
              p.pa0=sum(pa0),
              p.pb0=sum(pb0)) %>% data.frame()
  
  row.names(x12.smry) <-x12.smry$result
  
  
  table1 <- outer(x12.smry$p.pa1,x12.smry$p.pb1) 
  
  table0 <- outer(x12.smry$p.pa0,x12.smry$p.pb0)
  
 
  #colnames(table1) <- c("B.fail.stage1", "B.fail.stage2", "B.pass")
  #rownames(table1) <- c("A.fail.stage1", "A.fail.stage2", "A.pass")
  
  ## when both 2 doses pase stage 1 
  nr.x12 <- dim(x12)[1] 
  ab <- data.frame(a1=rep(x12$x1,nr.x12),
                   a2=rep(x12$x2,nr.x12),
                   b1=rep(x12$x1,each=nr.x12),
                   b2=rep(x12$x2,each=nr.x12))
  
  
  ab.pf <- ab %>% 
    mutate(a.out=case_when(a1<=r1 ~ "aF1",
                           a1>r1&(a1+a2)<=r ~ "aF2",
                           a1>r1 & (a1+a2)>r ~ "aP"),
           b.out=case_when(b1<=r1 ~ "bF1",
                           b1>r1&(b1+b2)<=r ~ "bF2",
                           b1>r1 & (b1+b2)>r ~ "bP")) %>% 
    mutate(bH.aPbP=ifelse(a.out=="aP"& b.out=="bP",1,0))
  
  ## filter only (a1+a2)<(b1+b2)
  bW.aPbP <- ab.pf %>% 
    filter(bH.aPbP==1)  %>% 
    filter((a1+a2)<(b1+b2))
  
  ## Indicator matrix: 
  indicator <- data.frame(a1 = bW.aPbP$a1,
                          a2 = bW.aPbP$a2,
                          b1 = bW.aPbP$b1,
                          b2 = bW.aPbP$b2,
                          a_total = bW.aPbP$a1+bW.aPbP$a2,
                          b_total = bW.aPbP$b1+bW.aPbP$b2,
                          posterior = 0)
  
  ## Sort the data frame by ascending order of a_total and b_total, 
  ## create groups for unique a_total, b_total
  indicator <- indicator %>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total) %>%
    mutate(sample_size = n)  %>%
    mutate(group = cur_group_id()) 
  
  
  ## load posterior matrix
  ## changed simulation in p.rbH to 10000 times
  posterior_matrix_filter <- posterior_matrix %>%
    filter(sample_size==n &a_total>r & b_total>r &a_total<b_total)%>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total)
  
  ## add the "posterior" values from the "posterior_matrix_filter" to 
  ## the "indicator" based on matching values of
  ##"a_total", "b_total", and "sample_size"
  cppFunction('
  DataFrame add_posterior_values(DataFrame indicator, DataFrame posterior_matrix_filter, int n) {
  int nrows = indicator.nrows();
  NumericVector group = indicator["group"];
  NumericVector a_total = indicator["a_total"];
  NumericVector b_total = indicator["b_total"];
  NumericVector posterior = indicator["posterior"];
  
  int num_groups = posterior_matrix_filter.nrows();
  NumericVector filter_a_total = posterior_matrix_filter["a_total"];
  NumericVector filter_b_total = posterior_matrix_filter["b_total"];
  NumericVector filter_sample_size = posterior_matrix_filter["sample_size"];
  NumericVector filter_posterior = posterior_matrix_filter["posterior"];
  
  for (int i = 0; i < num_groups; ++i) {
    double cur_posterior = filter_posterior[i];
    for (int j = 0; j < nrows; ++j) {
      if (group[j] == (i + 1) && a_total[j] == filter_a_total[i] && b_total[j] == filter_b_total[i] && n == filter_sample_size[i]) {
        posterior[j] = cur_posterior;
      }
    }
  }
  
  return DataFrame::create(_["a1"] = indicator["a1"],
                           _["a2"] = indicator["a2"],
                           _["b1"] = indicator["b1"],
                           _["b2"] = indicator["b2"],
                           _["a_total"] = indicator["a_total"],
                           _["b_total"] = indicator["b_total"],
                           _["posterior"] = posterior,    
                           _["group"] = indicator["group"]);
}

')
  indicator <- indicator %>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total) %>%
    mutate(sample_size = n)  %>%
    mutate(group = cur_group_id()) 
  
  posterior_matrix_filter <- posterior_matrix %>%
    filter(sample_size==n &a_total>r & b_total>r &a_total<b_total)%>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total)
  
  indicator <- add_posterior_values(indicator, posterior_matrix_filter, n)
  
  
  
  indicator_p <- indicator %>% 
    mutate(p.ra1rb1=dbinom(a1, n1, p=pa1)*dbinom(a2, n2, p=pa1)*
             dbinom(b1, n1, p=pb1)*dbinom(b2, n2, p=pb1),
           p.ra0rb0=dbinom(a1, n1, p=pa0)*dbinom(a2, n2, p=pa0)*
             dbinom(b1, n1, p=pb0)*dbinom(b2, n2, p=pb0))
  
  
  bW.aPbP.smry <- apply(indicator_p[indicator_p$posterior> delta,
                                    c("p.ra1rb1","p.ra0rb0")],2,sum)
  
  
  ## expected total sample size under p0
  ess <- n1*2*table0[1,1]+              
    (n1*2+n2*2)*(table0[2,2] +table0[2,3]+table0[3,2] +table0[3,3])+    
    (n1*2+n2)*(table0[1,2] +table0[1,3]+table0[2,1] +table0[3,1])  
  
  
  
  ## power and type I error of pick the true winner
  bpower_1 <- bW.aPbP.smry["p.ra1rb1"]
  bpower_0 <- bW.aPbP.smry["p.ra0rb0"]
  
  power <-  table1[1,3]+table1[2,3]+bpower_1 ## pa_1 pb_1 
  alpha <-  table0[1,3]+table0[2,3]+bpower_0 ## pa_0 pb_0 
  
  table1<- as.data.frame(table1)
  table1 <- round(table1,3)
  
  table1[3,3] <- paste0(table1[3,3], " (B claims as winner: ", round(bpower_1,3), ")")
  return(table1)
}


## return table 0
octable0 <- function(n, n1, r, r1, pa1, pb1, pa0, pb0, delta, posterior_matrix){
  
  if(pa1>pb1){
    stop("Error: pa1 must be less than pb1.")
  }
  
  n2 <- n-n1
  x12 <- data.frame(x1=rep(0:n1, each=n2+1), 
                    x2=rep(0:n2, n1+1)) 
  x12 <- x12 %>% filter(!(x1<=r1&x2>0)) %>%
    mutate(result=case_when(x1<=r1 ~ "F1", 
                            x1>r1&(x1+x2)<=r ~ "F2", 
                            x1>r1 & (x1+x2)>r ~ "P"),
           pa1=ifelse(x1<=r1, dbinom(x1, n1, p=pa1), 
                      dbinom(x1,n1, p=pa1)*dbinom(x2,n2, p=pa1)), 
           pb1=ifelse(x1<=r1, dbinom(x1, n1, p=pb1),
                      dbinom(x1,n1, p=pb1)*dbinom(x2,n2, p=pb1)),
           pa0=ifelse(x1<=r1, dbinom(x1, n1, p=pa0),
                      dbinom(x1,n1, p=pa0)*dbinom(x2,n2, p=pa0)),
           pb0=ifelse(x1<=r1, dbinom(x1, n1, p=pb0),
                      dbinom(x1,n1, p=pb0)*dbinom(x2,n2, p=pb0))
    )
  
  x12.smry <- x12[,-(1:2)] %>%
    group_by(result) %>%
    summarize(p.pa1=sum(pa1),
              p.pb1=sum(pb1),
              p.pa0=sum(pa0),
              p.pb0=sum(pb0)) %>% data.frame()
  
  row.names(x12.smry) <-x12.smry$result
  
  
  table1 <- outer(x12.smry$p.pa1,x12.smry$p.pb1) 
  
  table0 <- outer(x12.smry$p.pa0,x12.smry$p.pb0)
  
  
  #colnames(table0) <- c("B.fail.stage1", "B.fail.stage2", "B.pass")
  #rownames(table0) <- c("A.fail.stage1", "A.fail.stage2", "A.pass")
  ## when both 2 doses pase stage 1 
  nr.x12 <- dim(x12)[1] 
  ab <- data.frame(a1=rep(x12$x1,nr.x12),
                   a2=rep(x12$x2,nr.x12),
                   b1=rep(x12$x1,each=nr.x12),
                   b2=rep(x12$x2,each=nr.x12))
  
  
  ab.pf <- ab %>% 
    mutate(a.out=case_when(a1<=r1 ~ "aF1",
                           a1>r1&(a1+a2)<=r ~ "aF2",
                           a1>r1 & (a1+a2)>r ~ "aP"),
           b.out=case_when(b1<=r1 ~ "bF1",
                           b1>r1&(b1+b2)<=r ~ "bF2",
                           b1>r1 & (b1+b2)>r ~ "bP")) %>% 
    mutate(bH.aPbP=ifelse(a.out=="aP"& b.out=="bP",1,0))
  
  ## filter only (a1+a2)<(b1+b2)
  bW.aPbP <- ab.pf %>% 
    filter(bH.aPbP==1)  %>% 
    filter((a1+a2)<(b1+b2))
  
  ## Indicator matrix: 
  indicator <- data.frame(a1 = bW.aPbP$a1,
                          a2 = bW.aPbP$a2,
                          b1 = bW.aPbP$b1,
                          b2 = bW.aPbP$b2,
                          a_total = bW.aPbP$a1+bW.aPbP$a2,
                          b_total = bW.aPbP$b1+bW.aPbP$b2,
                          posterior = 0)
  
  ## Sort the data frame by ascending order of a_total and b_total, 
  ## create groups for unique a_total, b_total
  indicator <- indicator %>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total) %>%
    mutate(sample_size = n)  %>%
    mutate(group = cur_group_id()) 
  
  
  ## load posterior matrix
  ## changed simulation in p.rbH to 10000 times
  posterior_matrix_filter <- posterior_matrix %>%
    filter(sample_size==n &a_total>r & b_total>r &a_total<b_total)%>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total)
  
  ## add the "posterior" values from the "posterior_matrix_filter" to 
  ## the "indicator" based on matching values of
  ##"a_total", "b_total", and "sample_size"
  cppFunction('
  DataFrame add_posterior_values(DataFrame indicator, DataFrame posterior_matrix_filter, int n) {
  int nrows = indicator.nrows();
  NumericVector group = indicator["group"];
  NumericVector a_total = indicator["a_total"];
  NumericVector b_total = indicator["b_total"];
  NumericVector posterior = indicator["posterior"];
  
  int num_groups = posterior_matrix_filter.nrows();
  NumericVector filter_a_total = posterior_matrix_filter["a_total"];
  NumericVector filter_b_total = posterior_matrix_filter["b_total"];
  NumericVector filter_sample_size = posterior_matrix_filter["sample_size"];
  NumericVector filter_posterior = posterior_matrix_filter["posterior"];
  
  for (int i = 0; i < num_groups; ++i) {
    double cur_posterior = filter_posterior[i];
    for (int j = 0; j < nrows; ++j) {
      if (group[j] == (i + 1) && a_total[j] == filter_a_total[i] && b_total[j] == filter_b_total[i] && n == filter_sample_size[i]) {
        posterior[j] = cur_posterior;
      }
    }
  }
  
  return DataFrame::create(_["a1"] = indicator["a1"],
                           _["a2"] = indicator["a2"],
                           _["b1"] = indicator["b1"],
                           _["b2"] = indicator["b2"],
                           _["a_total"] = indicator["a_total"],
                           _["b_total"] = indicator["b_total"],
                           _["posterior"] = posterior,    
                           _["group"] = indicator["group"]);
}

')
  indicator <- indicator %>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total) %>%
    mutate(sample_size = n)  %>%
    mutate(group = cur_group_id()) 
  
  posterior_matrix_filter <- posterior_matrix %>%
    filter(sample_size==n &a_total>r & b_total>r &a_total<b_total)%>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total)
  
  indicator <- add_posterior_values(indicator, posterior_matrix_filter, n)
  
  
  
  indicator_p <- indicator %>% 
    mutate(p.ra1rb1=dbinom(a1, n1, p=pa1)*dbinom(a2, n2, p=pa1)*
             dbinom(b1, n1, p=pb1)*dbinom(b2, n2, p=pb1),
           p.ra0rb0=dbinom(a1, n1, p=pa0)*dbinom(a2, n2, p=pa0)*
             dbinom(b1, n1, p=pb0)*dbinom(b2, n2, p=pb0))
  
  
  bW.aPbP.smry <- apply(indicator_p[indicator_p$posterior> delta,
                                    c("p.ra1rb1","p.ra0rb0")],2,sum)
  
  
  ## expected total sample size under p0
  ess <- n1*2*table0[1,1]+              
    (n1*2+n2*2)*(table0[2,2] +table0[2,3]+table0[3,2] +table0[3,3])+    
    (n1*2+n2)*(table0[1,2] +table0[1,3]+table0[2,1] +table0[3,1])  
  
  
  
  ## power and type I error of pick the true winner
  bpower_1 <- bW.aPbP.smry["p.ra1rb1"]
  bpower_0 <- bW.aPbP.smry["p.ra0rb0"]
  
  power <-  table1[1,3]+table1[2,3]+bpower_1 ## pa_1 pb_1 
  alpha <-  table0[1,3]+table0[2,3]+bpower_0 ## pa_0 pb_0 
  
  table0<- as.data.frame(table0)
  table0 <- round(table0,3)
  
  table0[3,3] <- paste0(table0[3,3], " (B claims as winner: ", round(bpower_0,3), ")")
  return(table0)
  
}
















