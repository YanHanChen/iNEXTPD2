#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double Dq0(double n, double f1, double f2, double g1,  double g2) {
  double ans;
  if(f1 == 0) {ans = 0;
  } else if ( ((g1*f2)/(2*f1) < g2) | (f1==0)){
    ans = ((n-1)/n)*(pow(g1,2)/(2*g2));
  } else{
    ans = ((n-1)/n)*(g1*(f1-1)/(2*(f2+1)));
  }
  return(ans);
}
// [[Rcpp::export]]
double Dq1_2(double n, double g1, double A) {
  double q1 = 0;
  double h2 = 0;
  if(A==1||g1==0){
    h2 = 0;
  }else{
    for(int r = 1; r < n; r++){
      q1 = q1 + pow((1-A),r)/r;
    }
    h2 = (g1/n)*(pow(1-A,(-n+1)))*(-log(A)-q1);
  }
  return(h2);
}
// [[Rcpp::export]]
double Dq2(NumericMatrix tmpaL, double n, double t_bar) {
  double ans = 0;
  for(int i = 0; i < tmpaL.nrow(); i++){
    ans = ans + (tmpaL(i,2)*tmpaL(i,1)*tmpaL(i,0)*(tmpaL(i,0)-1)/n/(n-1));
  }
  //Rcout << "ans: " << ans;
  ans = pow(t_bar,2)/ans;
  return(ans);
}

// [[Rcpp::export]]
double Dq_2nd(double n, double g1, double A, double q) {
  double qq = 0;
  double ans = 0;
  if(A==1||g1==0){
    ans = 0;
  }else{
    for(int r = 0; r < n; r++){
      qq = qq + Rf_choose(q-1,r)*pow((A-1),r);
      //Rcpp::Rcout << "qq: " << qq << std::endl;
    }
    ans = (g1/n)*(pow(1-A,(-n+1)))*(pow(A,q-1)-qq);
  }
  return(ans);
}
// [[Rcpp::export]]
NumericVector Dq(NumericMatrix tmpaL, int n,NumericVector qs,double g1, double A, double t_bar){
  int nrows = tmpaL.nrow(), z = 0 , qlength = qs.length();
  double delta = 0.0;
  //NumericMatrix deltas(nrows,n);
  NumericMatrix ans_i(nrows,qlength);
  for(int i = 0; i<nrows; i++){
    z = tmpaL(i,0);
    for(int k = 0; k<=n-z; k++){
      delta = tmpaL(i,1)*Rf_dhyper( 1, z, n - z, k+1, false )/(k+1);
      //deltas(i,k) = Rf_dhyper( 1, z, n - z, k+1, false )/(k+1);
      for(int i_q = 0; i_q<qlength; i_q++){
        ans_i(i,i_q) = ans_i(i,i_q) + tmpaL(i,2) * Rf_choose(k-qs[i_q],k)*delta;
      }
    }
  }
  NumericVector ans_1(qlength);
  for(int i_q = 0; i_q<qlength; i_q++){
    for(int i = 0; i<nrows; i++){
      ans_1(i_q) = ans_1(i_q) + ans_i(i,i_q);
    }
    ans_1(i_q) = (ans_1(i_q) + Dq_2nd(n,g1,A,qs[i_q]))/pow(t_bar,qs[i_q]);
    ans_1(i_q) = pow(ans_1(i_q),1/(1-qs[i_q]));
  }
  return(ans_1);
}

// [[Rcpp::export]]
double delta(NumericMatrix del_tmpaL, double k, double n){
  double ans = 0;
  for(int i = 0;i < del_tmpaL.nrow(); i++){
    ans = ans +
      del_tmpaL(i,2)*del_tmpaL(i,1)*exp(Rf_lchoose(n-k-1,del_tmpaL(i,0)-1)-Rf_lchoose(n, del_tmpaL(i,0)));
  }
  return(ans);
}

// [[Rcpp::export]]
NumericVector delta_part2(NumericVector ai, double k, double n){
  int nrow = ai.length();
  NumericVector ans(nrow);
  for(int i = 0;i < nrow; i++){
    ans[i] = exp(Rf_lchoose(n-k-1,ai[i]-1)-Rf_lchoose(n, ai[i]));
  }
  return(ans);
}

// [[Rcpp::export]]
NumericVector RPD_old(NumericMatrix x , int n  , int m , NumericVector q) {
  int nrow = x.nrow();
  double tbar=0;
  NumericVector ghat(m);
  int qlength = q.length();
  NumericVector out(qlength);
  for (int i = 0; i < nrow; i++) {
    tbar += x(i, 0)*x(i, 1)/n;
  }
  //Rcpp::Rcout << "tbar in cpp: " << tbar << std::endl;
  for (int k = 0; k < m ; k++) {
    for (int i = 0; i < nrow; i++) {
      if ( x(i,0) >= k+1 && x(i,0) <= n-m+k+1 )
      {
        ghat[k] +=  x(i,1)*exp(Rf_lchoose(x(i,0), k+1)+Rf_lchoose(n-x(i,0), m-k-1)-Rf_lchoose(n, m)) ;
      }
      else
      {
        ghat[k] += 0 ;
      }
    }
  }
  //Rcpp::Rcout << "ghat: " << ghat << std::endl;
  for (int j = 0; j < qlength; j++ ){
    for(int k = 0; k < m; k++){
      if(q[j] == 0){
        out[j] = ghat[k] + out[j];
      }else if(q[j] == 1){
        //Rcout << "q1 in cpp: " <<log ( (k+1) )<< std::endl;
        out[j] = -( (k+1) / (m*tbar) ) * log ( (k+1) / (m*tbar) ) * ghat[k] + out[j];
      }else if(q[j] == 2){
        out[j] = pow( ( (k+1) / (m*tbar) ),2) * ghat[k] + out[j];
      }else{
        out[j] = pow( ( (k+1) / (m*tbar) ),q[j]) * ghat[k] + out[j];
      }
    }
  }
  for(int j = 0; j < qlength; j++ ){
    if(q[j] == 0){
      out[j] = out[j] ;
    }else if(q[j] == 1){
      out[j] = exp(out[j]);
    }else if(q[j] == 2){
      out[j] = 1 / out[j];
    }else{
      out[j] = pow( (out[j]) , 1/(1-q[j]) );
    }
  }
  return out ;
}

// [[Rcpp::export]]
NumericMatrix RPD(NumericVector ai ,NumericMatrix Lis, int n , int m , NumericVector q) {
  int tau_l = Lis.ncol();
  int S = ai.length();
  int qlength = q.length();
  NumericMatrix ghat_pt2(m,S);
  NumericMatrix out(qlength,tau_l);
  for (int k = 0; k < m ; k++) {
    for (int i = 0; i < S; i++) {
      if ( ai[i] >= k+1 && ai[i] <= n-m+k+1 )
      {
        ghat_pt2(k,i) = exp(Rf_lchoose(ai[i], k+1)+Rf_lchoose(n-ai[i], m-k-1)-Rf_lchoose(n, m)) ;
      }
      else
      {
        ghat_pt2(k,i) = 0 ;
      }
    }
  }
  for(int l = 0; l < tau_l;l++){
    double tbar=0;
    NumericVector ghat(m);
    for (int i = 0; i < S; i++) {
      tbar += ai[i]*Lis(i, l)/n;
    }
    for(int k = 0; k < m ; k++) {
      for (int i = 0; i < S; i++) {
        ghat[k] +=  Lis(i, l)*ghat_pt2(k,i) ;
      }
    }
    //Rcpp::Rcout << "ghat: " << ghat << std::endl;
    for(int j = 0; j < qlength; j++ ){
      for(int k = 0; k < m; k++){
        if(q[j] == 0){
          out(j,l) += ghat[k];
        }else if(q[j] == 1){
          //Rcout << "q1 in cpp: " <<log ( (k+1) )<< std::endl;
          out(j,l) += -( (k+1) / (m*tbar) ) * log ( (k+1) / (m*tbar) ) * ghat[k];
        }else if(q[j] == 2){
          out(j,l) += pow( ( (k+1) / (m*tbar) ),2) * ghat[k];
        }else{
          out(j,l) += pow( ( (k+1) / (m*tbar) ),q[j]) * ghat[k];
        }
      }
    }
    for(int j = 0; j < qlength; j++ ){
      if(q[j] == 0){
        out(j,l) = out(j,l) ;
      }else if(q[j] == 1){
        out(j,l) = exp(out(j,l));
      }else if(q[j] == 2){
        out(j,l) = 1 / out(j,l);
      }else{
        out(j,l) = pow( out(j,l) , 1/(1-q[j]) );
      }
    }
    
  }
  
  return out ;
}

// [[Rcpp::export]]
NumericMatrix ghat_pt2(NumericVector ai , int n , int mmax ) {
  //int tau_l = Lis.ncol();
  int S = ai.length();
  //int qlength = q.length();
  //NumericMatrix ghat_pt2(m,S);
  NumericMatrix out(mmax,S);
  for (int k = 0; k < mmax ; k++) {
    for (int i = 0; i < S; i++) {
      if ( ai[i] >= k+1 && ai[i] <= n-mmax+k+1 )
      {
        out(k,i) = exp(Rf_lchoose(ai[i], k+1)+Rf_lchoose(n-ai[i], mmax-k-1)-Rf_lchoose(n, mmax)) ;
      }
      else
      {
        out(k,i) = 0 ;
      }
    }
  }
  return out ;
}

// [[Rcpp::export]]
NumericVector qD_MLE(NumericVector q,NumericVector ai){
  const int length = q.size();
  const int S = ai.size();
  NumericVector Q(length);
  NumericVector temp(S);
  for(int j = 0; j<length;j++){
    for(int i = 0 ; i<S;i++){
      temp[i] = pow(ai[i],q[j]);
    }
    Q[j] = pow(sum(temp),1/(1-q[j]));
  }
  return Q;
}