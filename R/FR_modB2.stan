data{
  int N;
  int Nid;
  vector[N] Eaten;
  int ID[N];
  int E_If[N];
  int E_nIf[N];
  int R[N];
}
parameters{
  real a;
  // real b;
  real h;
  
  real a_E_If;
  real h_E_If;
  real a_E_nIf;
  real h_E_nIf;
 
  
  
  vector[Nid] a_Intercept;
// vector[Nid] b_Intercept;
  vector[Nid] h_Intercept;
  real<lower=0> sigma_ID;
  real sigma;
}

model{
  vector[N] FR;
  vector[N] A;
// vector[N] B;
  vector[N] H;
  
  sigma ~ cauchy( 0 , 2 );
  sigma_ID ~ cauchy( 0 , 2 );
  a_Intercept ~ normal( 0 , sigma_ID );
// b_Intercept ~ normal( 0 , sigma_ID );
  h_Intercept ~ normal( 0 , sigma_ID );

  a_E_If ~ normal( 0 , 10 );
  h_E_If ~ normal( 0 , 10 );
  
  a_E_nIf ~ normal( 0 , 10 );
  h_E_nIf ~ normal( 0 , 10 );

 
  for ( i in 1:N ) {
    A[i] = 10^(a + a_E_If * E_If[i] + a_E_nIf * E_nIf[i] + a_Intercept[ID[i]]);
    H[i] = 10^(h + h_E_If * E_If[i] + h_E_nIf * E_nIf[i] + h_Intercept[ID[i]]);

    FR[i] = A[i]*(R[i]) / (1 + A[i]*H[i]*(R[i]) );
  }
  Eaten ~ normal( FR , sigma );
}
