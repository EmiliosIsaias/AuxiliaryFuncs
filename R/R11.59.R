library(rethinking)

N <- 500 
#simulate family incomes for each individual  
family_income <- runif(N)  
#assign a unique coefficient for each type of event
b <- c(-2,0,2)  
career <- rep(NA,N)  #empty vector of choices for each individual  
for (i in 1:N) {
  score <- 0.5*(1:3) + b*family_income[i]  
  p <- softmax(score[1],score[2],score[3])  
  career[i] <- sample( 1:3, size=1, prob=p) 
  
}  

code_m11.14 <- "
data{
  int N; // number of observations 
  int K; // number of outcome values  
  array[N] int career; // outcome
  array[N] real family_income;  
}  

parameters{  
  vector[K-1] a; // intercepts
  vector[K-1] b; // coefficients on family income
}

model{  
  vector[K] p;  
  vector[K] s;  
  a ~ normal(0,1.5);  
  b ~ normal(0,1);  
  for (i in 1:N) {  
    for (j in 1:(K-1)) s[j] = a[j] + b[j]*family_income[i];  
    s[K] = 0; // the pivot  
    p = softmax( s);  
    career[i] ~ categorical( p);  
  }
}
"  
dat_list <- list(N=N, K=3, career=career, family_income=family_income)
m11.14 <- stan( model_code = code_m11.14, data = dat_list, chains = 4)
precis(m11.14, 2)
