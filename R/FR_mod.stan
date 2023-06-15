data{
  int N;
  int Nid;
  int Eaten[N];
  int ID[N];
  int Exp[N];
  int R[N];
}
parameters{
  real a;
  // real b;
  real h;
  
  real a_Exp;
 // real b_Exp;
  real h_Exp;
  
  
  vector[Nid] a_Intercept;
// vector[Nid] b_Intercept;
  vector[Nid] h_Intercept;
  real<lower=0> sigma_ID;
//   real sigma;
}

model{
  vector[N] FR;
  vector[N] A;
// vector[N] B;
  vector[N] H;
  
//   sigma ~ cauchy( 0 , 2 );
  sigma_ID ~ cauchy( 0 , 2 );
  a_Intercept ~ normal( 0 , sigma_ID );
// b_Intercept ~ normal( 0 , sigma_ID );
  h_Intercept ~ normal( 0 , sigma_ID );

  a_Exp ~ normal( 0 , 10 );
// b_Exp ~ normal( 0 , 10 );
  h_Exp ~ normal( 0 , 10 );
 
  for ( i in 1:N ) {
    A[i] = 10^(a + a_Exp * Exp[i] + a_Intercept[ID[i]]);
// B[i] = 10^(b + b_Exp * Exp[i] + b_Intercept[ID[i]]);
    H[i] = 10^(h + h_Exp * Exp[i] + h_Intercept[ID[i]]);


    FR[i] = A[i]*(R[i]) / (1 + A[i]*H[i]*(R[i]) );
  }
//   Eaten ~ normal( FR , sigma );

  Eaten ~ binomial( R , FR );
}
