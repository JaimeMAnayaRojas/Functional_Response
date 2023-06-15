// Functional response model
// We want to estiate the effects of parasites on it. 
// Expose and Infected = E_If
// Exposed and not infected = E_nIf
data{
  int N; 
  int Nid;
  vector[N] Eaten;
  array[N] int ID;
  array[N] int I;
  array[N] int nI;
  array[N] int R;
}
parameters{
  real a;
  real b;
  real h;
  
  real a_I;
  real h_I;
  real a_i;
  real h_i;
 
  
  
  vector[Nid] a_Intercept;
  vector[Nid] b_Intercept;
  vector[Nid] h_Intercept;
  real<lower=0> sigma_ID;
  real sigma;
}

model{
  vector[N] FR;
  vector[N] A;
  vector[N] B;
  vector[N] H;
  
  sigma ~ cauchy( 0 , 1 );
  sigma_ID ~ cauchy( 0 , 1 );
  a_Intercept ~ normal( 0 , sigma_ID );
  b_Intercept ~ normal( 0 , sigma_ID );
  h_Intercept ~ normal( 0 , sigma_ID );

  a ~ normal( 0 , 5 );
  h ~ normal( 0 , 5 );

  a_I ~ normal( 0 , 1 );
  h_I ~ normal( 0 , 1 );
  
  a_i ~ normal( 0 , 1 );
  h_i ~ normal( 0 , 1 );

 
  for ( i in 1:N ) {
    A[i] = exp(a + a_I * I[i] + a_i * nI[i] + a_Intercept[ID[i]]);
    H[i] = exp(h + h_I * I[i] + h_i * nI[i] + h_Intercept[ID[i]]);
    B[i] = exp(b + b_Intercept[ID[i]]);
    

    // FR[i] = (A[i]*(R[i])) / (1 + A[i]*H[i]*(R[i]) );
    FR[i] = (A[i]*(R[i]^B[i])) / (1 + A[i]*H[i]*(R[i]^B[i]) );
  }
  
  Eaten ~ normal( FR , sigma );
}
