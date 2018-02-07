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

// [[Rcpp::export]]
double maxd(NumericMatrix x, int d) {
  int p = x.nrow();
  double M=0;
  for(int i=0;i<(p-d);i++){
    for(int j=std::min(i+d,p-1);j<p;j++){
      M=std::max(std::abs(x(i,j)),M);//((abs(x(i,j))>M)?x(i,j):M);
    }
  }
  return M;
}

// [[Rcpp::export]]
NumericMatrix mdiag(int n, NumericVector d){
  NumericMatrix M(n,n);
  //std::fill(M.begin(), M.end(), 0);
  int m=d.size();
  for(int i=0;i<n;i++){
    for(int j=std::max(0,i-m+1);j<=std::min(n,i+m-1);j++){
      M(i,j)=d[std::abs(i-j)];
    }
  }
  return M;
}

// [[Rcpp::export]]
NumericMatrix mnorm(NumericMatrix M, int h, int type=1){
  int p=M.nrow();
  NumericMatrix N(p,p);//NumericMatrix N(clone(M));
  int il,iu,jl,ju;
  double s=0;
  for(int i=0;i<p;i++){
    for(int j=0;j<p;j++){
      il=std::max(i-h,0);
      iu=std::min(i+h,p-1);
      jl=std::max(j-h,0);
      ju=std::min(j+h,p-1);
      if(type==1){
        s=0;
        for(int i1=il;i1<=iu;i1++){
          for(int j1=jl;j1<=ju;j1++){
            s+=std::abs(M(i1,j1));
          }
        }
        N(i,j)=s/((iu-il+1)*(ju-jl+1));
      }else{
        s=0;
        for(int i1=il;i1<=iu;i1++){
          for(int j1=jl;j1<=ju;j1++){
            s=std::max(s,std::abs(M(i1,j1)));
          }
        }
        N(i,j)=s;
      }

      
    }
  }
  
  return N;
}


// [[Rcpp::export]]
NumericMatrix banding(NumericMatrix x, int k) {
  int p = x.nrow();
  NumericMatrix M(p,p);

  for(int i=0;i<p;i++){
    for(int j=0;j<p;j++){
      if(std::abs(i-j)<=k){
        M(i,j)=x(i,j);
      }
    }
  }
  return M;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

# /*** R
# timesTwo(42)
# */
