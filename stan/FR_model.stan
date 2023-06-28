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
  array[N] int i;
  array[N] int R;
  // array[N] int pD;
  // vector[N] cSize;

}
parameters{
  real a;
  real a_I; // effect of infection
  real a_i; // effect of exposure but not infection
  // real a_pD; // effect of previous diet

  real h;
  real h_i;
  real h_I;
 
  real b;
  // real b_I;
  // real b_i;
  
  
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

 
  for ( j in 1:N ) {
    A[j] = exp(a + a_I * I[j] + a_i *i[j] + a_Intercept[ID[j]]); // + a_pD *pD[j] 
    H[j] = exp(h + h_I * I[j] + h_i * i[j] + h_Intercept[ID[j]]);
    B[j] = exp(b + b_Intercept[ID[j]]);
    

    // FR[i] = (A[i]*(R[i])) / (1 + A[i]*H[i]*(R[i]) );
    FR[j] = (A[j]*(R[j]^B[j])) / (1 + A[j]*H[j]*(R[j]^B[j]) );
  }
  
  Eaten ~ normal( FR , sigma );
}
