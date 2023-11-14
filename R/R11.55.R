library(rethinking)
#simulate career choices among 500 individuals
N <- 500 #number of individuals
income <- c(1,2,5) #expected income of each career
score <- 0.5*income #scores for each career, based on income  

#next line converts scores to probabilities
p <- softmax(score[1],score[2],score[3])  #now simulate choice
#outcome career holds event type values, not counts
career <- rep(NA,N) #empty vector of choices for each individual
#sample chosen career for each individual  
set.seed(34302)  
for (i in 1:N) career[i] <- sample( 1:3, size=1, prob=p) 

code_m11.13 <- "  
data{
  int N; // number of individuals
  int K; // number of possible careers
  array[N] int career; // outcome
  vector[K] career_income;  
}  

parameters{  
  vector[K-1] a; // intercepts  
  real<lower=0> b; // association of income with choice
} 

model{  
  vector[K] p;
  vector[K] s;  
  a ~ normal( 0, 1);  
  b ~ normal( 0, 0.5);  
  s[1] = a[1] + b*career_income[1];  
  s[2] = a[2] + b*career_income[2];  
  s[3] = 0; // pivot
  p = softmax( s);
  career ~ categorical( p);  
}  
"

dat_list <- list( N=N, K=3, career=career, career_income=income)
m11.13 <- stan( model_code=code_m11.13, data=dat_list, chains=4, cores=4)  
precis( m11.13, 2) 

post <- extract.samples( m11.13)
#set up logit scores
s1 <- with( post, a[,1] + b*income[1])
s2_orig <- with( post, a[,2] + b*income[2])
s2_new <- with( post, a[,2] + b*income[2]*2)  

#compute probabilities for original and counterfactual
p_orig <- sapply( 1:length(post$b), function(i)  
  softmax( c(s1[i],s2_orig[i],0)))  
p_new <- sapply( 1:length(post$b), function(i)  
  softmax( c(s1[i],s2_new[i],0)))  
#summarize  
p_diff <- p_new[2,] - p_orig[2,]  
precis( p_diff)