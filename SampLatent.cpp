#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
// Z = A matrix of testing responses whose jth row is of the form 
// (col1=Z_j, col2=cj, col3=assay used, col4:col(4+cj-1)=indices of the individuals assigned to the jth pool) 
// Y = matrix whose ith row is of the form 
// (col1=0, col2=number of pools ith individual involved in, col3:col(3+col2-1)=pools ith individual was assigned to)

// [[Rcpp::export]]
NumericVector SampLatent(int N,NumericVector p,NumericMatrix Y,NumericMatrix Z,
                         NumericVector U,NumericVector se,NumericVector sp,int na) {
  int i, np, l, j, Zj, cj, ybar, tid, id, t;
  float pi1, pi2, pistar, sej, spj;
  NumericVector WW(N);
  
  for(i=0;i<N;i++){
    pi1=p(i);
    pi2=1-p(i);
    np=Y(i,2-1);
    for(l=0;l<np;l++){
      j=Y(i,(l+2));
      Zj=Z(j-1,0);
      cj=Z(j-1,1);
      tid=Z(j-1,2);
      sej=se(tid-1);
      spj=sp(tid-1);
      ybar=0;
      Y(i,1-1)=0;
      for(t=0;t<cj;t++){
        id=Z(j-1,(t+3));
        ybar=ybar+Y(id-1,1-1);
      }
      pi1=pi1*(sej*Zj + (1-sej)*(1-Zj));
      if(ybar > 0){
        pi2=pi2*(sej*Zj + (1-sej)*(1-Zj));
      }else{
        pi2=pi2*((1-spj)*Zj + spj*(1-Zj));
      }
    }
    pistar=(pi1/(pi1+pi2));
    p(i)=pistar	;
    if(U(i)<pistar){
      Y(i,1-1)=1;
    }else{Y(i,1-1)=0;}
    WW(i)=Y(i,1-1);
  }  

  return WW;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


