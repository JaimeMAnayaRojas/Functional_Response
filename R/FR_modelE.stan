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
  array[N] int E;
  array[N] int R;
  vector[N] pD; 
  vector[N] Egg;
  // vector[N] cSize;
  // vector[N] ISize;
}
parameters{
  real a;
  real a_E; // effect of infection
  real a_pD; // effect of previous diet
  real a_Egg;
  
  real h;
  real h_E;
  real h_pD; // effect of previous diet
  real h_Egg;
  
  real b;
  real b_E;

  vector[Nid] a_Intercept;
  vector[Nid] b_Intercept;
  vector[Nid] h_Intercept;
  
  vector[Nid] ab_Intercept;
  vector[Nid] bb_Intercept;
  vector[Nid] hb_Intercept;
  
  real<lower=0> sigma_ID;
  real<lower=0> sigma_B;
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
  
  ab_Intercept ~ normal( 0 , sigma_B );
  bb_Intercept ~ normal( 0 , sigma_B );
  hb_Intercept ~ normal( 0 , sigma_B );

  a ~ normal( 0 , 5 );
  h ~ normal( 0 , 5 );
  b ~ normal( 0 , 5 );

  a_E ~ normal( 0 , 1 );
  a_E ~ normal( 0 , 1 );
  h_E ~ normal( 0 , 1 );


  a_pD ~ normal( 0 , 1 );
  h_pD ~ normal( 0 , 1 );

  a_Egg ~ normal( 0 , 1 );
  h_Egg ~ normal( 0 , 1 );
  // b_Egg ~ normal(0, 1 );
 
  for ( j in 1:N ) {
    A[j] = exp(a + a_E * E[j] + a_Intercept[ID[j]] + ab_Intercept[Block[j]] + a_pD *pD[j] + a_Egg *Egg[j]);// + a_pD *pD[j] + a_cSize *cSize[j] + a_ISize *ISize[j] + a_Egg *Egg[j]
    H[j] = exp(h + h_E * E[j] + h_Intercept[ID[j]] +  hb_Intercept[Block[j]] + h_pD *pD[j] + h_Egg *Egg[j]);//+ h_pD *pD[j] + h_cSize *cSize[j] + h_ISize *ISize[j] + h_Egg *Egg[j]
    B[j] = exp(b + bb_Intercept[ID[j]] + bb_Intercept[Block[j]]);
    
    FR[j] = (A[j]*(R[j]^B[j])) / (1 + A[j]*H[j]*(R[j]^B[j]) );
  }
  
  Eaten ~ normal( FR , sigma );
}
