// Functional response model
// We want to estiate the effects of parasites on it. 
// Expose and Infected = E_If
// Exposed and not infected = E_nIf
data{
  int N; 
  int Nid;
  vector[N] Eaten;
  array[N] int ID;
  array[N] int Block;
  array[N] int I;
  array[N] int i;
  array[N] int R;
  vector[N] pD; 
  vector[N] Egg;
  // vector[N] cSize;
  // vector[N] ISize;
}
parameters{
  real a;
  real a_I; // effect of infection
  real a_i; // effect of exposure but not infection
  real a_pD; // effect of previous diet
  // real a_cSize;
  // real a_ISize;
  real a_Egg;
  
  real h;
  real h_i;
  real h_I;
  real h_pD; // effect of previous diet
  // real h_cSize;
  // real h_ISize;
  real h_Egg;
  
  

  // real b_pD; // effect of previous diet
  // real b_cSize;
  // real b_ISize;
  // real b_Egg; 
  
  vector[Nid] a_Intercept;
  
  vector[Nid] h_Intercept;
  
  vector[Nid] ab_Intercept;
  
  vector[Nid] hb_Intercept;
  
  real<lower=0> sigma_ID;
  // real<lower=0> sigma_B;
  real sigma;
}

model{
  vector[N] FR;
  vector[N] A;
  // vector[N] B;
  vector[N] H;
  
  sigma ~ cauchy( 0 , 1 );
  sigma_ID ~ cauchy( 0 , 1 );
  // sigma_B ~ cauchy( 0 , 1 );
  a_Intercept ~ normal( 0 , sigma_ID );
  h_Intercept ~ normal( 0 , sigma_ID );
  
  // ab_Intercept ~ normal( 0 , sigma_B );
  
  // hb_Intercept ~ normal( 0 , sigma_B );

  a ~ normal( 0 , 1 );
  h ~ normal( 0 , 1 );

  a_I ~ normal( 0 , 1 );
  h_I ~ normal( 0 , 1 );
  
  a_i ~ normal( 0 , 1 );
  h_i ~ normal( 0 , 1 );

  a_pD ~ normal( 0 , 1 );
  h_pD ~ normal( 0 , 1 );
  
  a_Egg ~ normal( 0 , 1 );
  h_Egg ~ normal( 0 , 1 );
 
  for ( j in 1:N ) {
    A[j] = exp(a + a_I * I[j] + a_i * i[j] + a_Intercept[ID[j]] + a_pD *pD[j] + a_Egg *Egg[j]);
    H[j] = exp(h + h_I * I[j] + h_i * i[j] + h_Intercept[ID[j]] + h_pD *pD[j] + h_Egg *Egg[j]);
    FR[j] = (A[j]*(R[j])) / (1 + A[j]*H[j]*(R[j]) );
  }
  
  Eaten ~ normal( FR , sigma );
}
