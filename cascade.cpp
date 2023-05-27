// [[Rcpp::plugins("cpp11")]]
#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <RcppArmadillo.h>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]

////////////////////////////////////////////////

// Object conversion
Rcpp::NumericMatrix arma_to_R(arma::mat m){
  
  Rcpp::NumericMatrix r_m(m.n_rows,m.n_cols);
  std::copy(m.begin(),m.end(),r_m.begin());
  
  return(r_m);
}
 
Rcpp::NumericVector arma_to_R(arma::vec v){
  
  Rcpp::NumericVector r_v(v.size());
  std::copy(v.begin(),v.end(),r_v.begin());
  
  return(r_v);
}
 
Rcpp::NumericVector arma_to_R(arma::uvec v){
  
  Rcpp::NumericVector r_v(v.size());
  std::copy(v.begin(),v.end(),r_v.begin());
  
  return(r_v);
}
 
////////////////////////////////////////////////
 
 
// out neighbors
arma::uvec out_neighbors(arma::mat Y,
                               int i){
  if(i<0 || i>=Y.n_cols){
    std::cout << "i must be a number in 0,1,2,...," << Y.n_cols - 1 << "." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  return arma::find(Y.row(i));
}

// graph statistics
int n_state(arma::vec X, int state){
  return sum(X == state);
}

int graph_edges_from_i_to_j(arma::mat Y,arma::vec X, int i, int j){
  
  arma::uvec in_i = (X == i);
  arma::uvec in_j = (X == j);
  
  return arma::as_scalar( in_i.t()*Y*in_j );
}


double graph_border_size(arma::mat Y,arma::vec X){
  
  int from_S_to_R = graph_edges_from_i_to_j(Y,X,3,2);
  int from_R_to_S = graph_edges_from_i_to_j(Y,X,2,3);
  
  return (from_S_to_R + from_R_to_S)/arma::accu(Y);
}

double graph_assortativity(arma::mat Y,arma::vec X){
  
  arma::mat f(4,4,arma::fill::zeros);
  for(int i = 0; i<4; i++)
    for(int j = 0; j<4; j++){
      f(i,j) = graph_edges_from_i_to_j(Y,X,i,j);
    }
  f = f/arma::accu(Y); // Convertir totales a proporciones
  
  double expected_vals = arma::as_scalar( 
    arma::sum(f,0)*arma::sum(f,1) 
  );
  
  return (arma::sum(f.diag()) - expected_vals)/(1 - expected_vals);
}

////////////////////////////////////////////////


// Search minimum
class minimum_jump{
  
  int u, v;
  double min_t = std::numeric_limits<double>::infinity();
  std::vector<double> jump_times;
  std::vector<double> jump_type;
    
  
  public:
    
    void update_min(double t,int new_u, int new_v){
      if(t < min_t){
        min_t = t;
        u = new_u;
        v = new_v;
      }
    }
    
    bool make_jump(arma::vec& X,std::ofstream& file){
      
      if(min_t == std::numeric_limits<double>::infinity()){
        
        return 0;
        
      }else{
        
        switch((int)X(v)){
        case 0:
          X(v) = 1;
          file << v << "," << 1 << "," << 0 << ",";
          break;
        case 1:
          X(v) = X(u);
          file << v << "," << (int)X(u) << "," << 1 << ",";
          break;
        default:
          X(v) = 1;
          file << v << "," << 1 << "," << 1 << ",";
          break;
        }
        
        file << min_t << ","; 
        min_t = std::numeric_limits<double>::infinity();
        
        return 1;
      }
      
    }
    
};

void update_stats(arma::mat Y, arma::vec X,
                  std::ofstream& file){
  file << n_state(X,0) << ",";
  file << n_state(X,1) << ",";
  file << n_state(X,2) << ",";
  file << n_state(X,3) << ",";
  file << graph_border_size(Y,X) << ",";
  file << graph_assortativity(Y,X) << "\n";
}
    

// Progress report
class progress_report{
  public:
  private:
    int bar_size;
    int N;
    int last_bar_fill;
  public:
    progress_report(int new_N){
      Rcpp::Function opt("options");
      Rcpp::List WD = opt("width"); 
      bar_size      = std::max(Rcpp::as<int>(WD["width"]) - 16,5);
      N             = new_N;
      last_bar_fill = 0;
      
      std::cout << "|S v R|: " << std::string(bar_size,'|') << "  0% \r" << std::flush;
    }
    void display(arma::vec X){
      int size_RvS = n_state(X,2) + n_state(X,3);
      int bar_fill = (size_RvS*bar_size)/N;
      
      if(bar_fill != last_bar_fill){
        last_bar_fill = bar_fill;
        int perc      = (size_RvS*100)/N;
        
        std::cout << "|S v R|: " << std::string(bar_fill,'/') <<
          std::string(bar_size - bar_fill,'|') << "  " <<
            perc << "% \r" << std::flush;
      }
      
    };
    void close(){
      std::cout << std::endl;
    }
};

class check_stationarity{
  private:
    int I,U,R,S;
    int tol;
  public:
    int streak;
    check_stationarity(arma::vec X,int new_tol){
      I = n_state(X,0);
      U = n_state(X,1);
      R = n_state(X,2);
      S = n_state(X,3);
      streak = 0;
      tol = new_tol;
    }
    void monitor(arma::vec X){
      if( (std::fabs(I - n_state(X,0))<tol) &&
          (std::fabs(U - n_state(X,1))<tol) &&
          (std::fabs(R - n_state(X,2))<tol) &&
          (std::fabs(S - n_state(X,3))<tol)
      ){
        streak++;
      }else{
        streak = 0;
        I = n_state(X,0);
        U = n_state(X,1);
        R = n_state(X,2);
        S = n_state(X,3);
      }
    }
};

////////////////////////////////////////////////


// [[Rcpp::export]]
Rcpp::NumericVector cpp_cascade(Rcpp::NumericMatrix rY,
                       Rcpp::NumericVector rX,
                       Rcpp::NumericVector rO,
                       Rcpp::NumericVector rI,
                       Rcpp::NumericMatrix rtau,
                       int seed,
                       std::string statsFile){
  
  // Initialize file
  std::ofstream outStats(statsFile);
  outStats << "node,destiny,jump_type,jump_len,I,U,R,S,borderSize,Assortativity\n";
  
  // Conversion to C++
  arma::mat Y    = Rcpp::as<arma::mat>(rY);
  arma::vec X    = Rcpp::as<arma::vec>(rX);
  arma::vec O    = Rcpp::as<arma::vec>(rO);
  arma::vec I    = Rcpp::as<arma::vec>(rI);
  arma::mat tau  = Rcpp::as<arma::mat>(rtau);
  
  // Initialize variables
  int N = Y.n_cols;
  double t;
  int u, v;
  minimum_jump min_jump; // argmin search
  int n_jumps = 0;
  progress_report RvS(N); // Active links monitor
  check_stationarity alarm(X,N/20);
  
  // Call R base rexp function
  Rcpp::Function rexp("rexp");
  Rcpp::Function set_seed("set.seed");
  
    
  // Iterate over jumps
  set_seed(seed);
  do{
    
    // Simulate jump times
    for(u = 0; u < N; u++){
      
      if(X(u) == 0) continue;
      
      arma::uvec N_out = out_neighbors(Y,u);
      
      for(int v_i = 0;v_i < N_out.size(); v_i++){
        v = N_out(v_i);
        
        if(X(v) == 0){ // => Information jump
          t = Rcpp::as<double>(rexp(Rcpp::Named("n",1),
                                    Rcpp::Named("rate", std::exp(O(u)))
                                    ));
          min_jump.update_min(t,u,v);
          continue;
        }
        
        if(X(u) == 1) continue;
        
        if(X(u) != X(v)){ // -> Influence jump
          t = Rcpp::as<double>(rexp(Rcpp::Named("n",1),
                                    Rcpp::Named("rate", std::exp(O(u) + tau(u,v)*I(v)))
                                    ));
          min_jump.update_min(t,u,v);
        }
      }
    }
    
    RvS.display(X);
    
    
    // Make jump and save time
    if(!min_jump.make_jump(X,outStats)){
      break;
    };
		// n_jumps++;
		
		alarm.monitor(X);
    
    // Compute statistics
    // stats_monitor.update_stats(Y,X,outStats);
    update_stats(Y,X,outStats);
    
    // Stop loops
    // if((n_state(X,2) + n_state(X,3))>0.8*N){
    //   SR_loop++;
    // }else{
    //   SR_loop = 0;
    // }
    
  // } while (SR_loop<5*N || n_state(X,0)>0.05*N); // Choosing N or any factor of N is an arbitrary decition that was taken based in observed results.
  } while (alarm.streak<3*N); // Choosing N or any factor of N is an arbitrary decition that was taken based in observed results.
  
  RvS.close();
  outStats.close();
  
  return arma_to_R(X);
}

// Aparentemente este algoritmo es O( N^2(N-1)*ln(N)^4 )
// El término log es producto de los algoritmos ya
// implementados en Armadillo y C++. Mi algoritmo como
// tal es O(N²(N-1)).

/////////////////////////////////////////////////


// [[Rcpp::export]]
Rcpp::NumericVector cpp_update_X(Rcpp::NumericVector rX,
                             Rcpp::NumericMatrix rJumps,
                             int until){
  // Convert Rcpp to arma
  arma::vec X = Rcpp::as<arma::vec>(rX);
  arma::mat Jumps = Rcpp::as<arma::mat>(rJumps);
  
  for(int i = 0; i<until;i++){
    X( Jumps(i,0) ) = Jumps(i,1);
  }
  return arma_to_R(X);
}

/*** R
cascade_simulator <- function(Y,X,O,I,
                              tau = (Y+t(Y))/2,
                              seed = 2023,
                              stats_file = ".tmp_results.txt"){
  X <- cpp_cascade(Y,X,O,I,tau,seed,stats_file)
  return(X)
}
update_X_until <- function(X,cascade,until = nrow(cascade)){
  Jumps <- as.matrix(cascade[,1:2])
  X <- cpp_update_X(X,Jumps,until-1)
  return(X)
}
*/