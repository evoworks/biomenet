#include "InferenceFunctions.h"

//unsigned int myRmultinom( NumericVector probs ){                                                                                                                                  
//  unsigned int rmultinomDraw = 0;                                                                                                                                                 
//                                                                                                                                                                                  
//  return rmultinomDraw ;                                                                                                                                                          
//}                                                                                                                                                                                 


#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/distributions/normal.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <vector>
#include <algorithm>
#include <numeric>

#include <iostream>
#include <fstream>

#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_errno.h>

#include <math.h>


using namespace boost::numeric::ublas;
using namespace boost::numeric;

boost::mt19937 gen(static_cast<unsigned int>(std::time(0)));


/*`
  Now define a function that simulates rolling a loaded die.
  Note that the C++0x library contains a `discrete_distribution`
  class which would be a better way to do this.
*/
int roll_weighted_die( const std::vector<double> probabilities) {
  std::vector<double> cumulative;
  std::partial_sum(probabilities.begin(), probabilities.end(),
		   std::back_inserter(cumulative));
  boost::uniform_real<> dist(0, cumulative.back());
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > die(gen, dist);

  return (std::lower_bound(cumulative.begin(), cumulative.end(), die()) - cumulative.begin()) ;
}



int Cwhich(IntegerVector binVec){
  int vecsize = binVec.size();
  for ( int i =0; i < vecsize; i++)
    if (binVec(i) == 1)
      return i;
  cout << "soemthing went wrong, no value equal to one was found !" << endl;
  return -1;
}


double Alpha_phi_Metropolis_Hastings_update(const double alpha_phi_val, const Rcpp::IntegerMatrix C_k_l, Rcpp::IntegerVector sum_C_k_l_row){

  //Rcpp::IntegerMatrix C_k_l(K, L); Rcpp::IntegerVector sum_C_k_l_row = rep(0,K); const matrix<unsigned int> Count_G_V, const std::vector<unsigned int> sum_Count_G_V

  //RNGScope scope;
  boost::normal_distribution<> nd(0.0, 1.0);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor(gen, nd);
  //var_nor.engine().seed(static_cast<unsigned int>(std::time(NULL) + getpid()));
  //var_nor.distribution().reset();
  boost::uniform_real<> dist(0, 1);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > unif_rand(gen, dist);
  using boost::math::normal;
  normal s;

  double SearchWidth = 0.1;
  unsigned G = C_k_l.ncol(); // Count_G_V.size1();
  unsigned V = C_k_l.nrow(); // Count_G_V.size2();  
  double pAccept = 0;
  bool Accept = 0;
  double AcceptanceRate = 0;
  double RET_alpha_phi_val = alpha_phi_val;
  double NEW_alpha_phi_val = RET_alpha_phi_val + var_nor() * RET_alpha_phi_val*SearchWidth;

  double RET_sum_alpha_phi = RET_alpha_phi_val * G;
  double NEW_sum_alpha_phi = NEW_alpha_phi_val * G;

  pAccept += V*( lgamma(NEW_sum_alpha_phi) - G*lgamma(NEW_alpha_phi_val) );
  pAccept -= V*( lgamma(RET_sum_alpha_phi) - G*lgamma(RET_alpha_phi_val) );
  for (unsigned v=0; v < V; v++){
    pAccept -= lgamma(sum_C_k_l_row[v] + NEW_sum_alpha_phi);
    pAccept += lgamma(sum_C_k_l_row[v] + RET_sum_alpha_phi);
    for (unsigned g=0; g < G; g++){
      pAccept += ( lgamma(NEW_alpha_phi_val + C_k_l(v, g)) );
      pAccept -= ( lgamma(RET_alpha_phi_val + C_k_l(v, g)) );
    }
  }
  pAccept = exp(pAccept);
  pAccept *= ( pdf(s, (RET_alpha_phi_val - NEW_alpha_phi_val) /(SearchWidth*NEW_alpha_phi_val) ) ) / ( pdf(s, (NEW_alpha_phi_val - RET_alpha_phi_val) /(SearchWidth*RET_alpha_phi_val) ) );
  if (pAccept >= 1)
    Accept = 1;
  else
    Accept = (unif_rand() < pAccept);

  if (Accept)
    return NEW_alpha_phi_val;
  else
    return RET_alpha_phi_val;
  
  
}
// alpha_phi must be replaced with alpha_delta, another similar function for alpha_gamma
double Alpha_delta_Metropolis_Hastings_update(const double alpha_delta, const Rcpp::IntegerMatrix N_S_c_l, const Rcpp::IntegerVector N_S_l){


  //Rcpp::IntegerMatrix N_S_c_l(N, L); Rcpp::IntegerVector N_S_l = rep(0, L); const matrix<unsigned int> Count_T_G, const std::vector<unsigned int> sum_Count_T_G
  //Rcpp::IntegerMatrix N_R_c_l(N, L); Rcpp::IntegerVector N_R_l = rep(0, L); 

  boost::normal_distribution<> nd(0.0, 1.0);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor(gen, nd);
  //var_nor.engine().seed(static_cast<unsigned int>(std::time(NULL) + getpid()));
  //var_nor.distribution().reset();
  boost::uniform_real<> dist(0, 1);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > unif_rand(gen, dist);
  using boost::math::normal;
  normal s;

  double SearchWidth = 0.1;
  double RET_alpha_delta = alpha_delta;
  unsigned T = N_S_c_l.nrow(); //Count_T_G.size1();
  unsigned G = N_S_c_l.ncol(); //Count_T_G.size2();
  double RET_sum_alpha_delta = alpha_delta * T;

  double pAccept = 0;
  bool Accept = 0;
  double NEW_alpha_delta = RET_alpha_delta + var_nor() * RET_alpha_delta*SearchWidth;
  double NEW_sum_alpha_delta = NEW_alpha_delta * T;

  pAccept += G*( lgamma(NEW_sum_alpha_delta) - T*lgamma(NEW_alpha_delta) );
  pAccept -= G*( lgamma(RET_sum_alpha_delta) - T*lgamma(RET_alpha_delta) );
  for (unsigned g=0; g < G; g++){
    pAccept -= lgamma(N_S_l[g] + NEW_sum_alpha_delta);
    pAccept += lgamma(N_S_l[g] + RET_sum_alpha_delta);
    for (unsigned t=0; t < T; t++){
      pAccept += ( lgamma(NEW_alpha_delta + N_S_c_l(t, g)) );
      pAccept -= ( lgamma(RET_alpha_delta + N_S_c_l(t, g)) );
    }
  }

  pAccept = exp(pAccept);
  pAccept *= ( pdf(s, (RET_alpha_delta - NEW_alpha_delta) /(SearchWidth*NEW_alpha_delta) ) ) / ( pdf(s, (NEW_alpha_delta - RET_alpha_delta) /(SearchWidth*RET_alpha_delta) ) );
  if (pAccept >= 1)
    Accept = 1;
  else
    Accept = (unif_rand() < pAccept);
  
  if (Accept)
    RET_alpha_delta = NEW_alpha_delta;
  
  return RET_alpha_delta;
}

double Alpha_gamma_Metropolis_Hastings_update(const double alpha_gamma, const Rcpp::IntegerMatrix N_R_c_l, const Rcpp::IntegerVector N_R_l){


  //Rcpp::IntegerMatrix N_S_c_l(N, L); Rcpp::IntegerVector N_S_l = rep(0, L); const matrix<unsigned int> Count_T_G, const std::vector<unsigned int> sum_Count_T_G
  //Rcpp::IntegerMatrix N_R_c_l(N, L); Rcpp::IntegerVector N_R_l = rep(0, L); 

  boost::normal_distribution<> nd(0.0, 1.0);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor(gen, nd);
  //var_nor.engine().seed(static_cast<unsigned int>(std::time(NULL) + getpid()));
  //var_nor.distribution().reset();
  boost::uniform_real<> dist(0, 1);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > unif_rand(gen, dist);
  using boost::math::normal;
  normal s;

  double SearchWidth = 0.1;
  double RET_alpha_gamma = alpha_gamma;
  unsigned T = N_R_c_l.nrow(); //Count_T_G.size1();
  unsigned G = N_R_c_l.ncol(); //Count_T_G.size2();
  double RET_sum_alpha_gamma = alpha_gamma * T;

  double pAccept = 0;
  bool Accept = 0;
  double NEW_alpha_gamma = RET_alpha_gamma + var_nor() * RET_alpha_gamma*SearchWidth;
  double NEW_sum_alpha_gamma = NEW_alpha_gamma * T;

  pAccept += G*( lgamma(NEW_sum_alpha_gamma) - T*lgamma(NEW_alpha_gamma) );
  pAccept -= G*( lgamma(RET_sum_alpha_gamma) - T*lgamma(RET_alpha_gamma) );
  for (unsigned g=0; g < G; g++){
    pAccept -= lgamma(N_R_l[g] + NEW_sum_alpha_gamma);
    pAccept += lgamma(N_R_l[g] + RET_sum_alpha_gamma);
    for (unsigned t=0; t < T; t++){
      pAccept += ( lgamma(NEW_alpha_gamma + N_R_c_l(t, g)) );
      pAccept -= ( lgamma(RET_alpha_gamma + N_R_c_l(t, g)) );
    }
  }

  pAccept = exp(pAccept);
  pAccept *= ( pdf(s, (RET_alpha_gamma - NEW_alpha_gamma) /(SearchWidth*NEW_alpha_gamma) ) ) / ( pdf(s, (NEW_alpha_gamma - RET_alpha_gamma) /(SearchWidth*RET_alpha_gamma) ) );
  if (pAccept >= 1)
    Accept = 1;
  else
    Accept = (unif_rand() < pAccept);
  
  if (Accept)
    RET_alpha_gamma = NEW_alpha_gamma;
  
  return RET_alpha_gamma;
}


//alpha_pi will be replaced with alpha_theta
double Alpha_theta_Metropolis_Hastings_update(const double alpha_theta, const Rcpp::IntegerMatrix C_g_k, const Rcpp::IntegerVector sum_C_g_k_row){
  
  //Rcpp::IntegerMatrix C_g_k(G, K); Rcpp::IntegerVector sum_C_g_k_row = rep(0,G); const matrix<unsigned int> Count_S_V, const std::vector<unsigned int> sum_Count_S_V
  

  boost::normal_distribution<> nd(0.0, 1.0);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor(gen, nd);
  //var_nor.engine().seed(static_cast<unsigned int>(std::time(NULL) + getpid()));
  //var_nor.distribution().reset();
  boost::uniform_real<> dist(0, 1);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > unif_rand(gen, dist);
  using boost::math::normal;
  normal s;

  double SearchWidth = 0.1;
  double RET_alpha_theta = alpha_theta;
  unsigned N = C_g_k.nrow(); 
  unsigned V = C_g_k.ncol(); 

  double pAccept = 0;
  bool Accept = 0;
  double NEW_alpha_theta = RET_alpha_theta + var_nor() * RET_alpha_theta*SearchWidth;
  
  for (unsigned n=0; n<N; n++){
    pAccept += ( lgamma( V * NEW_alpha_theta) - V * lgamma( NEW_alpha_theta) );
    pAccept -= ( lgamma( V * RET_alpha_theta) - V * lgamma( RET_alpha_theta) );

    pAccept -= lgamma(sum_C_g_k_row[n] + V * NEW_alpha_theta);
    pAccept += lgamma(sum_C_g_k_row[n] + V * RET_alpha_theta);

    for (unsigned x=0; x<V; x++){
      pAccept += lgamma( C_g_k(n, x) + NEW_alpha_theta);
      pAccept -= lgamma( C_g_k(n, x) + RET_alpha_theta);
    }
  }

  pAccept = exp(pAccept);
  pAccept *= ( pdf(s, (RET_alpha_theta - NEW_alpha_theta) /(SearchWidth*NEW_alpha_theta) ) ) / ( pdf(s, (NEW_alpha_theta - RET_alpha_theta) /(SearchWidth*RET_alpha_theta) ) );
  if (pAccept >= 1)
    Accept = 1;
  else
    Accept = (unif_rand() < pAccept);
  
  if (Accept)
    RET_alpha_theta = NEW_alpha_theta;
  
  return RET_alpha_theta;
  
}


// Metropolis-Hastings within Gibbs Sampling
RcppExport SEXP doMetropolisWithinGibbs(SEXP Rsenders, SEXP Rrecipients, SEXP RNodes, SEXP RK, SEXP RL, SEXP RmaxItr, SEXP RY, SEXP RZ, SEXP RBurnin, SEXP RLag){
  Environment stats("package:stats");
  Function rmultinom = stats["rmultinom"];
  Rcpp::RNGScope scope;

  Rcpp::List S(Rsenders);
  Rcpp::List R(Rrecipients);
  Rcpp::IntegerVector Nodes(RNodes);
  //vector<string> Nodes = as< vector<string> >(RNodes);

  Rcpp::List Y(RY);
  Rcpp::List Z(RZ);
  unsigned int K = Rcpp::as<int>(RK);
  unsigned int L = Rcpp::as<int>(RL);
  unsigned int max_iter = Rcpp::as<int>(RmaxItr);
  unsigned int burnin = Rcpp::as<int>(RBurnin);
  unsigned int iter_lag = Rcpp::as<int>(RLag);

  cout << "Initializing \n";                                                                                                                           

  //Nodes = as.character( union( unique(senders), unique(recipients) ));                                                                      
                                                                                                                               
  unsigned int  G = S.size(); //  G = length(S);
  cout << G << " networks\n";
  unsigned int N = Nodes.size(); //  N = length(Nodes);
  cout << N << " nodes\n";

  
  Rcpp::NumericVector alpha_theta = rep(0.1, K); //  alpha.theta = rep(0.1, K);
  double sum_alpha_theta = sum(alpha_theta);
  Rcpp::NumericVector alpha_phi = rep(0.1, L); //  alpha.phi = rep(0.1, L);
  double sum_alpha_phi = sum(alpha_phi);
  Rcpp::NumericVector alpha_delta = rep(0.1, N); //  alpha.delta = rep(0.1, N);
  double sum_alpha_delta = sum(alpha_delta);
  Rcpp::NumericVector alpha_gamma = rep(0.1, N); //  alpha.gamma = rep(0.1, N);
  double sum_alpha_gamma = sum(alpha_gamma);
 
//  # random starting values for Y and Z if not given
  Rcpp::IntegerVector reaction_count = rep(0, G); //  reaction.count = list();
  std::vector< std::vector<int> > interaction_count; interaction_count.resize(G);//Rcpp::List interaction_count; //  interaction.count = list();
  if ( Y.size() == 0){ //Y.isNULL()){ //  if (is.null(Y)){
    cout << "Starting a new Gibbs run, initializing group variables";
    cout << " for " << max_iter << " iterations, burnin iterations = " << burnin << " , saving every " << iter_lag << " iterations" << endl;

    std::vector<double> probZ(K, 1.0/K);
    std::vector<double> probY(L, 1.0/L);
    for ( unsigned cidx=0; cidx < G; cidx++) { //    for (c.idx in 1:G){
      reaction_count(cidx) = as<Rcpp::List>(S(cidx)).size(); //      reaction.count[[c.idx]] = length(S[[c.idx]]);
      cout << reaction_count(cidx) << " reactions for network " << cidx << endl;
      interaction_count[cidx].resize( reaction_count(cidx) );//interaction_count(cidx) = rep(0, reaction_count(cidx) ); //      interaction.count[[c.idx]] = rep(0, reaction.count[[c.idx]] );
      Z.insert(cidx, rep(0, reaction_count(cidx) )); //      Z[[c.idx]] = rep(0, reaction.count[[c.idx]] );
      Y.insert(cidx, rep(0, reaction_count(cidx) ));  //      Y[[c.idx]] = rep(0, reaction.count[[c.idx]] );  
      for ( unsigned r=0; r< reaction_count(cidx); r++ ){ //      for ( r in 1:(reaction.count[[c.idx]]) ){
	as<Rcpp::IntegerVector>(Z(cidx))(r) = roll_weighted_die(probZ); //Cwhich(rmultinom(1,1, rep(1.0/K,K) )); //        Z[[c.idx]][r] = which(rmultinom(1,1, rep(1,K)/K )==1);  
	as<Rcpp::IntegerVector>(Y(cidx))(r) = roll_weighted_die(probY); //Cwhich(rmultinom(1,1, rep(1.0/L,L) )); //        Y[[c.idx]][r] = which(rmultinom(1,1, rep(1,L)/L )==1);  
      } //      }
    } //    }
  } else { //  }
    //burnin = 0;
    //Rcpp::List Y = clone( RY );
    //Rcpp::List Z = clone( RZ );
    cout << "Initializing from a previous run";
    cout << " for " << max_iter << " iterations, saving every " << iter_lag << " iterations" << endl;
    for ( unsigned cidx=0; cidx < G; cidx++) { 
      reaction_count(cidx) = as<Rcpp::List>(S(cidx)).size(); 
      interaction_count[cidx].resize( reaction_count(cidx) );
    }
  }

  
//  # initializing counter variable;
  cout << "Initializing and filling counter structures" << endl;
  Rcpp::IntegerMatrix C_g_k(G, K); for ( unsigned cidx=0; cidx < G; cidx++) for (unsigned k=0; k< K; k++) C_g_k(cidx, k) = 0; //  C.g.k = matrix(0, nrow=G, ncol=K);
  Rcpp::IntegerVector sum_C_g_k_row = rep(0,G); 
  Rcpp::IntegerMatrix C_k_l(K, L); for (unsigned k=0; k< K; k++) for (unsigned l=0; l< L; l++) C_k_l(k,l) = 0; //  C.k.l = matrix(0, nrow=K, ncol=L);
  Rcpp::IntegerVector sum_C_k_l_row = rep(0,K); 
  Rcpp::IntegerMatrix N_S_c_l(N, L); for(unsigned n=0; n< N; n++) for(unsigned l=0; l< L; l++) N_S_c_l(n,l) = 0; //  N.S.c.l = matrix(0, nrow=N, ncol=L);
  Rcpp::IntegerMatrix N_R_c_l(N, L); for(unsigned n=0; n< N; n++) for(unsigned l=0; l< L; l++) N_R_c_l(n,l) = 0; //  N.R.c.l = matrix(0, nrow=N, ncol=L);
  Rcpp::IntegerVector N_S_l = rep(0, L); //  N.S.l = rep(0, L);
  Rcpp::IntegerVector N_R_l = rep(0, L); //  N.R.l = rep(0, L);

  std::vector< std::vector< std::vector<int> > > S_count; 
  std::vector< std::vector< std::vector<int> > > R_count; 
  std::vector< std::vector< std::list<int> > > S_lst;
  std::vector< std::vector< std::list<int> > > R_lst;

  for ( unsigned cidx=0; cidx < G; cidx++) { //  for (c.idx in 1:G){

    std::vector< std::vector<int> >  cur_G_S_count; 
    std::vector< std::vector<int> >  cur_G_R_count; 
    std::vector< std::list<int> >  cur_G_S_lst;
    std::vector< std::list<int> >  cur_G_R_lst;
    for ( unsigned r=0; r< reaction_count(cidx); r++ ){ //    for ( r in 1:(reaction.count[[c.idx]]) ){
      int curZ = as<Rcpp::IntegerVector>(Z(cidx))(r); //      cur.Z = Z[[c.idx]][r] ;
      int curY = as<Rcpp::IntegerVector>(Y(cidx))(r); //      cur.Y = Y[[c.idx]][r] ;
      interaction_count[cidx][r] =as<Rcpp::List>( as<Rcpp::List>( S(cidx) ) (r) ).size();//      interaction.count[[c.idx]][r] = length(S[[c.idx]][[r]]);
      C_g_k(cidx, curZ) = C_g_k(cidx, curZ) +1; //      C.g.k[c.idx, cur.Z] = C.g.k[c.idx, cur.Z] +1;
      C_k_l(curZ, curY) = C_k_l(curZ, curY) +1; //      C.k.l[cur.Z, cur.Y] = C.k.l[cur.Z, cur.Y] +1;

      std::vector<int> cur_S_count(N,0) ;
      std::vector<int> cur_R_count(N,0) ;
      std::list<int> cur_S_lst;
      std::list<int> cur_R_lst;
      for ( unsigned i=0; i< interaction_count[cidx][r]; i++ ){ //      for ( i in 1:(interaction.count[[c.idx]][r]) ){
	unsigned curS = as<Rcpp::IntegerVector>( as<Rcpp::List>( S(cidx) ) (r) ) (i) ; //        cur.S = S[[c.idx]] [[r]] [i];
	cur_S_lst.push_back(curS);
	cur_S_count[ curS ] = cur_S_count[ curS ] +1;
	N_S_c_l( curS, curY ) = N_S_c_l( curS, curY ) +1; //        N.S.c.l[ cur.S, cur.Y ] = N.S.c.l[ cur.S, cur.Y ] +1;
	N_S_l(curY) = N_S_l(curY) + 1; //        N.S.l[cur.Y] = N.S.l[cur.Y] + 1;

	unsigned curR = as<Rcpp::IntegerVector>( as<Rcpp::List>( R(cidx) ) (r) ) (i) ; //        cur.R = R[[c.idx]] [[r]] [i];
	cur_R_lst.push_back(curR);
	cur_R_count[ curR ] = cur_R_count[ curR ] +1;
	N_R_c_l( curR, curY ) = N_R_c_l( curR, curY ) +1; //        N.R.c.l[ cur.R, cur.Y ] = N.R.c.l[ cur.R, cur.Y ] +1;
	N_R_l(curY) = N_R_l(curY) +1; //        N.R.l[cur.Y] = N.R.l[cur.Y] +1;
      } //      }
      cur_S_lst.sort(); cur_S_lst.unique(); 
      cur_R_lst.sort(); cur_R_lst.unique();

      std::vector<int> cur_S_count_lst;
      std::vector<int> cur_R_count_lst;
      for (list<int>::iterator n_it=cur_S_lst.begin(); n_it != cur_S_lst.end(); n_it++){
	unsigned int n = *n_it;
	cur_S_count_lst.push_back( cur_S_count[n] );
      }
      for (list<int>::iterator n_it=cur_R_lst.begin(); n_it != cur_R_lst.end(); n_it++){
	unsigned int n = *n_it;
	cur_R_count_lst.push_back( cur_R_count[n] );
      }
      cur_G_S_count.push_back(cur_S_count_lst);
      cur_G_R_count.push_back(cur_R_count_lst);
      cur_G_S_lst.push_back(cur_S_lst);
      cur_G_R_lst.push_back(cur_R_lst);

    } //    }
    S_count.push_back(cur_G_S_count);
    R_count.push_back(cur_G_R_count);
    S_lst.push_back(cur_G_S_lst);
    R_lst.push_back(cur_G_R_lst);

    sum_C_g_k_row(cidx)  = sum ( C_g_k.row(cidx) );
    
  } //  }
  for (unsigned kidx=0; kidx < K; kidx++)
    sum_C_k_l_row(kidx) = sum( C_k_l.row(kidx) );
  

  cout << "Starting Gibbs Iterations" << endl; //  # start Gibbs
  unsigned ndraws = (max_iter - burnin) / iter_lag + 1;
  Rcpp::List allZ(ndraws); //Rcpp::List allZ;
  Rcpp::List allY(ndraws); //Rcpp::List allY;
  for ( unsigned itr=1; itr <= max_iter; itr++ ){ //  for ( itr in 1:max.iter){
    //if (itr % 10 == 0){ //    if (itr %% 5 == 0){
      cout << "iteration " << itr << endl; //      cat(paste("iteration #", itr, "\n"));
      //} //    }

    for ( unsigned cidx=0; cidx < G; cidx++) { //    for (c.idx in 1:G){
      //cout << "sampling for network " << cidx << endl;
      for ( unsigned r=0; r< reaction_count(cidx); r++ ){ //      for ( r in 1:(reaction.count[[c.idx]]) ){
	int curZ = as<Rcpp::IntegerVector>(Z(cidx))(r) ; //        cur.Z = Z[[c.idx]][r] ;
	int curY = as<Rcpp::IntegerVector>(Y(cidx))(r) ; //        cur.Y = Y[[c.idx]][r] ;
	
	C_g_k(cidx, curZ) = C_g_k(cidx, curZ) -1; //        C.g.k[c.idx, cur.Z] = C.g.k[c.idx, cur.Z] -1;
	sum_C_g_k_row(cidx) = sum_C_g_k_row(cidx) -1;

	C_k_l(curZ, curY) = C_k_l(curZ, curY) -1; //        C.k.l[cur.Z, cur.Y] = C.k.l[cur.Z, cur.Y] -1;
	sum_C_k_l_row(curZ) = sum_C_k_l_row(curZ) -1;

	unsigned cur_interaction_count = interaction_count[cidx][r]; //        cur.interaction.count = interaction.count[[c.idx]][r];
	for ( unsigned i=0; i<cur_interaction_count; i++ ){ //        for ( i in 1:cur.interaction.count ){
	  unsigned curS = as<Rcpp::IntegerVector>( as<Rcpp::List>( S(cidx) ) (r) ) (i) ; //          cur.S = S[[c.idx]] [[r]] [i];
	  N_S_c_l( curS, curY ) = N_S_c_l( curS, curY ) -1; //          N.S.c.l[ cur.S, cur.Y ] = N.S.c.l[ cur.S, cur.Y ] -1;
	  N_S_l(curY) = N_S_l(curY) - 1; //          N.S.l[cur.Y] = N.S.l[cur.Y] - 1;

	  unsigned curR = as<Rcpp::IntegerVector>( as<Rcpp::List>( R(cidx) ) (r) ) (i) ; //          cur.R = R[[c.idx]] [[r]] [i];
	  N_R_c_l( curR, curY ) = N_R_c_l( curR, curY ) -1; //          N.R.c.l[ cur.R, cur.Y ] = N.R.c.l[ cur.R, cur.Y ] -1;
	  N_R_l(curY) = N_R_l(curY) - 1; //          N.R.l[cur.Y] = N.R.l[cur.Y] - 1;
	} //        }
		
	std::vector<double> prob_vec(K*L, 1.0/(K*L)); //Rcpp::NumericVector prob_vec = rep(1.0/(K*L), K*L); //        prob.vec = rep(0, K*L);
	
	for (unsigned k=0; k<K; k++){ //        for (k in 1:K){
	  for (unsigned l=0; l<L; l++){ //          for(l in 1:L){
	    //prob_vec(k*L + l) = ( alpha_theta(k) + C_g_k(cidx, k) ) / ( sum( alpha_theta ) + sum ( C_g_k.row(cidx) ) ) * ( alpha_phi(l) + C_k_l(k,l) ) / ( sum( alpha_phi ) + sum( C_k_l.row(k) ) ); //            prob.vec[(k-1)*L + l] = ( alpha.theta[k] + C.g.k[c.idx, k] ) / ( sum( alpha.theta + C.g.k[c.idx, ] ) ) * ( alpha.phi[l] + C.k.l[k,l] ) / ( sum( alpha.phi + C.k.l[k,] ) );
	    //prob_vec(k*L + l) = ( alpha_theta(k) + C_g_k(cidx, k) ) / ( sum_alpha_theta  + sum_C_g_k_row(cidx)  ) * ( alpha_phi(l) + C_k_l(k,l) ) / ( sum_alpha_phi  + sum_C_k_l_row(k)  ); 
	    prob_vec[k*L + l] = ( alpha_theta(k) + C_g_k(cidx, k) ) / ( sum_alpha_theta  + sum_C_g_k_row(cidx)  ) * ( alpha_phi(l) + C_k_l(k,l) ) / ( sum_alpha_phi  + sum_C_k_l_row(k)  ); 

	    unsigned iiii=0;
	    for (list<int>::iterator n_it=S_lst[cidx][r].begin(); n_it != S_lst[cidx][r].end(); n_it++){ //            for ( n in cur.S.lst){
	      unsigned n = *n_it;
	      for (unsigned c=0; c < S_count[cidx][r][iiii]; c++){ //for (unsigned c=0; c < cur_S_count(n); c++){ //              for ( c in 1:cur.S.count[ n ] ){
		//prob_vec(k*L + l) = prob_vec(k*L + l) *  ( alpha_delta(n) + N_S_c_l(n, l) + c ); //                prob.vec[(k-1)*L + l] = prob.vec[(k-1)*L + l] * ( alpha.delta[n] + N.S.c.l[n, l] + c -1);
		prob_vec[k*L + l] = prob_vec[k*L + l] *  ( alpha_delta(n) + N_S_c_l(n, l) + c ); 
	      } //              }
	      iiii++;
	    } //            }
	    for (unsigned c=0; c< cur_interaction_count; c++){ //            for ( n in 1:cur.interaction.count ){ # check this n or c as index
	      //prob_vec(k*L + l) = prob_vec(k*L + l) / ( sum_alpha_delta + N_S_l(l) +c ); //              prob.vec[(k-1)*L + l] = prob.vec[(k-1)*L + l] / ( sum(alpha.delta) + N.S.l[l] +c -1 )
	      prob_vec[k*L + l] = prob_vec[k*L + l] / ( sum_alpha_delta + N_S_l(l) +c );
	    } //            }

	    iiii=0;
	    for (list<int>::iterator n_it=R_lst[cidx][r].begin(); n_it != R_lst[cidx][r].end(); n_it++){ //            for ( n in cur.R.lst){
	      unsigned n = *n_it;
	      for (unsigned c=0; c < R_count[cidx][r][iiii]; c++){ //for (unsigned c=0; c < cur_R_count(n); c++){ //              for ( c in 1:cur.R.count[ n ] ){
		//prob_vec(k*L + l) = prob_vec(k*L + l) * ( alpha_gamma(n) + N_R_c_l(n, l) +c ); //                prob.vec[(k-1)*L + l] = prob.vec[(k-1)*L + l] * ( alpha.gamma[n] + N.R.c.l[n, l] +c -1);
		prob_vec[k*L + l] = prob_vec[k*L + l] * ( alpha_gamma(n) + N_R_c_l(n, l) +c );
	      } //              }
	      iiii++;
	    } //            }
	    for (unsigned c=0; c< cur_interaction_count; c++){ //            for ( n in 1:cur.interaction.count ){
	      //prob_vec(k*L + l) = prob_vec(k*L + l) / ( sum_alpha_gamma + N_R_l(l) +c ); //              prob.vec[(k-1)*L + l] = prob.vec[(k-1)*L + l] / ( sum(alpha.gamma) + N.R.l[l] +c -1)
	      prob_vec[k*L + l] = prob_vec[k*L + l] / ( sum_alpha_gamma + N_R_l(l) +c );
	    } //            }
	    //          
	  } //          }
	} //        }
	//
	unsigned ZY = roll_weighted_die( prob_vec ); //unsigned ZY = Cwhich(rmultinom(1,1, prob_vec)); //        ZY = which(rmultinom(1,1, prob.vec)==1);
	curY = ZY % L; //        cur.Y = ZY %% L;
	curZ = (ZY - curY)/L; //        cur.Z = (ZY - cur.Y)/L;
        
	as<Rcpp::IntegerVector>(Z(cidx))(r) = curZ; //        Z[[c.idx]][r] = cur.Z;
	as<Rcpp::IntegerVector>(Y(cidx))(r) = curY; //        Y[[c.idx]][r] = cur.Y;

	//savedGibbsIterations(drawcount-1) =  Rcpp::List::create(Rcpp::Named("Z")=clone(Z), Rcpp::Named("X")=clone(X) ) ;
	
	//        
	for ( unsigned i=0; i< cur_interaction_count; i++){//        for ( i in 1:cur.interaction.count ){
	  unsigned curS = Rcpp::as<int>( as<Rcpp::List>( as<Rcpp::List>( S(cidx) ) (r) ) (i) ); //          cur.S = S[[c.idx]] [[r]] [i];
	  N_S_c_l( curS, curY ) = N_S_c_l( curS, curY ) +1; //          N.S.c.l[ cur.S, cur.Y ] = N.S.c.l[ cur.S, cur.Y ] +1;
	  N_S_l(curY) = N_S_l(curY) + 1; //          N.S.l[cur.Y] = N.S.l[cur.Y] + 1;
	  //
	  unsigned curR = Rcpp::as<int>( as<Rcpp::List>( as<Rcpp::List>( R(cidx) ) (r) ) (i) ); //          cur.R = R[[c.idx]] [[r]] [i];
	  N_R_c_l( curR, curY ) = N_R_c_l( curR, curY ) +1; //          N.R.c.l[ cur.R, cur.Y ] = N.R.c.l[ cur.R, cur.Y ] +1;
	  N_R_l(curY) = N_R_l(curY) + 1; //          N.R.l[cur.Y] = N.R.l[cur.Y] + 1;
	} //        }
	
	C_g_k(cidx, curZ) = C_g_k(cidx, curZ) +1; //        C.g.k[c.idx, cur.Z] = C.g.k[c.idx, cur.Z] +1;
	sum_C_g_k_row(cidx) = sum_C_g_k_row(cidx) +1;
	C_k_l(curZ, curY) = C_k_l(curZ, curY) +1; //        C.k.l[cur.Z, cur.Y] = C.k.l[cur.Z, cur.Y] +1;
	sum_C_k_l_row(curZ) = sum_C_k_l_row(curZ) +1;
	//      
	//      
      } //      }
    } //    }
    
    // SAMPLING HYPERPARAMETERS WITH METROPOLIS-HASTINGS UPDATES
    double new_alpha_gamma_val = Alpha_gamma_Metropolis_Hastings_update(alpha_gamma[0], N_R_c_l, N_R_l);
    alpha_gamma = rep(new_alpha_gamma_val, N);
    sum_alpha_gamma = sum(alpha_gamma);
    double new_alpha_delta_val = Alpha_delta_Metropolis_Hastings_update(alpha_delta[0], N_S_c_l, N_S_l);
    alpha_delta = rep(new_alpha_delta_val, N);
    sum_alpha_delta = sum(alpha_delta);
    double new_alpha_theta_val = Alpha_theta_Metropolis_Hastings_update(alpha_theta[0], C_g_k, sum_C_g_k_row);
    alpha_theta = rep(new_alpha_theta_val, K);
    sum_alpha_theta = sum(alpha_theta);
    double new_alpha_phi_val = Alpha_phi_Metropolis_Hastings_update(alpha_phi[0], C_k_l, sum_C_k_l_row);
    alpha_phi = rep(new_alpha_phi_val, L);
    sum_alpha_phi = sum(alpha_phi);

    //cout << "before inserting new sample in the list" << endl;
    if ( (itr >= burnin) && ( (itr - burnin) % iter_lag == 0)){
      allZ((itr-burnin)/iter_lag) = clone(Z); //allZ.insert( (itr-burnin)/iter_lag, Z);
      allY((itr-burnin)/iter_lag) = clone(Y); //allY.insert( (itr-burnin)/iter_lag, Y);
    }
    //cout << "after inserting new sample in the list" << endl;
    //
  } //  }
  return Rcpp::List::create(Rcpp::Named("Z") = allZ, Rcpp::Named("Y") = allY); //  list(Zs = Z, Ys = Y);

} 








//# Gibbs sampling implementation

RcppExport SEXP doGibbs(SEXP Rsenders, SEXP Rrecipients, SEXP RNodes, SEXP RK, SEXP RL, SEXP RmaxItr, SEXP RY, SEXP RZ, SEXP RBurnin, SEXP RLag){ //doGibbs <- function(S, R, Nodes, K, L, max.iter, Y=NULL, Z=NULL){
  Environment stats("package:stats");
  Function rmultinom = stats["rmultinom"];
  Rcpp::RNGScope scope;

  Rcpp::List S(Rsenders);
  Rcpp::List R(Rrecipients);
  Rcpp::IntegerVector Nodes(RNodes);
  //vector<string> Nodes = as< vector<string> >(RNodes);

  Rcpp::List Y(RY);
  Rcpp::List Z(RZ);
  unsigned int K = Rcpp::as<int>(RK);
  unsigned int L = Rcpp::as<int>(RL);
  unsigned int max_iter = Rcpp::as<int>(RmaxItr);
  unsigned int burnin = Rcpp::as<int>(RBurnin);
  unsigned int iter_lag = Rcpp::as<int>(RLag);

  cout << "Initializing \n";                                                                                                                           

  //Nodes = as.character( union( unique(senders), unique(recipients) ));                                                                      
                                                                                                                               
  unsigned int  G = S.size(); //  G = length(S);
  cout << G << " networks\n";
  unsigned int N = Nodes.size(); //  N = length(Nodes);
  cout << N << " nodes\n";

  
  Rcpp::NumericVector alpha_theta = rep(0.01, K); //  alpha.theta = rep(0.1, K);
  double sum_alpha_theta = sum(alpha_theta);
  Rcpp::NumericVector alpha_phi = rep(0.01, L); //  alpha.phi = rep(0.1, L);
  double sum_alpha_phi = sum(alpha_phi);
  Rcpp::NumericVector alpha_delta = rep(0.01, N); //  alpha.delta = rep(0.1, N);
  double sum_alpha_delta = sum(alpha_delta);
  Rcpp::NumericVector alpha_gamma = rep(0.01, N); //  alpha.gamma = rep(0.1, N);
  double sum_alpha_gamma = sum(alpha_gamma);
 
//  # random starting values for Y and Z if not given
  Rcpp::IntegerVector reaction_count = rep(0, G); //  reaction.count = list();
  std::vector< std::vector<int> > interaction_count; interaction_count.resize(G);//Rcpp::List interaction_count; //  interaction.count = list();
  if ( Y.size() == 0){ //Y.isNULL()){ //  if (is.null(Y)){
    cout << "Starting a new Gibbs run, initializing group variables";
    cout << " for " << max_iter << " iterations, burnin iterations = " << burnin << " , saving every " << iter_lag << " iterations" << endl;

    std::vector<double> probZ(K, 1.0/K);
    std::vector<double> probY(L, 1.0/L);
    for ( unsigned cidx=0; cidx < G; cidx++) { //    for (c.idx in 1:G){
      reaction_count(cidx) = as<Rcpp::List>(S(cidx)).size(); //      reaction.count[[c.idx]] = length(S[[c.idx]]);
      cout << reaction_count(cidx) << " reactions for network " << cidx << endl;
      interaction_count[cidx].resize( reaction_count(cidx) );//interaction_count(cidx) = rep(0, reaction_count(cidx) ); //      interaction.count[[c.idx]] = rep(0, reaction.count[[c.idx]] );
      Z.insert(cidx, rep(0, reaction_count(cidx) )); //      Z[[c.idx]] = rep(0, reaction.count[[c.idx]] );
      Y.insert(cidx, rep(0, reaction_count(cidx) ));  //      Y[[c.idx]] = rep(0, reaction.count[[c.idx]] );  
      for ( unsigned r=0; r< reaction_count(cidx); r++ ){ //      for ( r in 1:(reaction.count[[c.idx]]) ){
	as<Rcpp::IntegerVector>(Z(cidx))(r) = roll_weighted_die(probZ); //Cwhich(rmultinom(1,1, rep(1.0/K,K) )); //        Z[[c.idx]][r] = which(rmultinom(1,1, rep(1,K)/K )==1);  
	as<Rcpp::IntegerVector>(Y(cidx))(r) = roll_weighted_die(probY); //Cwhich(rmultinom(1,1, rep(1.0/L,L) )); //        Y[[c.idx]][r] = which(rmultinom(1,1, rep(1,L)/L )==1);  
      } //      }
    } //    }
  } else { //  }
    //burnin = 0;
    //Rcpp::List Y = clone( RY );
    //Rcpp::List Z = clone( RZ );
    cout << "Initializing from a previous run";
    cout << " for " << max_iter << " iterations, saving every " << iter_lag << " iterations" << endl;
    for ( unsigned cidx=0; cidx < G; cidx++) { 
      reaction_count(cidx) = as<Rcpp::List>(S(cidx)).size(); 
      interaction_count[cidx].resize( reaction_count(cidx) );
    }
  }

  
//  # initializing counter variable;
  cout << "Initializing and filling counter structures" << endl;
  Rcpp::IntegerMatrix C_g_k(G, K); for ( unsigned cidx=0; cidx < G; cidx++) for (unsigned k=0; k< K; k++) C_g_k(cidx, k) = 0; //  C.g.k = matrix(0, nrow=G, ncol=K);
  Rcpp::IntegerVector sum_C_g_k_row = rep(0,G); 
  Rcpp::IntegerMatrix C_k_l(K, L); for (unsigned k=0; k< K; k++) for (unsigned l=0; l< L; l++) C_k_l(k,l) = 0; //  C.k.l = matrix(0, nrow=K, ncol=L);
  Rcpp::IntegerVector sum_C_k_l_row = rep(0,K); 
  Rcpp::IntegerMatrix N_S_c_l(N, L); for(unsigned n=0; n< N; n++) for(unsigned l=0; l< L; l++) N_S_c_l(n,l) = 0; //  N.S.c.l = matrix(0, nrow=N, ncol=L);
  Rcpp::IntegerMatrix N_R_c_l(N, L); for(unsigned n=0; n< N; n++) for(unsigned l=0; l< L; l++) N_R_c_l(n,l) = 0; //  N.R.c.l = matrix(0, nrow=N, ncol=L);
  Rcpp::IntegerVector N_S_l = rep(0, L); //  N.S.l = rep(0, L);
  Rcpp::IntegerVector N_R_l = rep(0, L); //  N.R.l = rep(0, L);

  std::vector< std::vector< std::vector<int> > > S_count; 
  std::vector< std::vector< std::vector<int> > > R_count; 
  std::vector< std::vector< std::list<int> > > S_lst;
  std::vector< std::vector< std::list<int> > > R_lst;

  for ( unsigned cidx=0; cidx < G; cidx++) { //  for (c.idx in 1:G){

    std::vector< std::vector<int> >  cur_G_S_count; 
    std::vector< std::vector<int> >  cur_G_R_count; 
    std::vector< std::list<int> >  cur_G_S_lst;
    std::vector< std::list<int> >  cur_G_R_lst;
    for ( unsigned r=0; r< reaction_count(cidx); r++ ){ //    for ( r in 1:(reaction.count[[c.idx]]) ){
      int curZ = as<Rcpp::IntegerVector>(Z(cidx))(r); //      cur.Z = Z[[c.idx]][r] ;
      int curY = as<Rcpp::IntegerVector>(Y(cidx))(r); //      cur.Y = Y[[c.idx]][r] ;
      interaction_count[cidx][r] =as<Rcpp::List>( as<Rcpp::List>( S(cidx) ) (r) ).size();//      interaction.count[[c.idx]][r] = length(S[[c.idx]][[r]]);
      C_g_k(cidx, curZ) = C_g_k(cidx, curZ) +1; //      C.g.k[c.idx, cur.Z] = C.g.k[c.idx, cur.Z] +1;
      C_k_l(curZ, curY) = C_k_l(curZ, curY) +1; //      C.k.l[cur.Z, cur.Y] = C.k.l[cur.Z, cur.Y] +1;

      std::vector<int> cur_S_count(N,0) ;
      std::vector<int> cur_R_count(N,0) ;
      std::list<int> cur_S_lst;
      std::list<int> cur_R_lst;
      for ( unsigned i=0; i< interaction_count[cidx][r]; i++ ){ //      for ( i in 1:(interaction.count[[c.idx]][r]) ){
	unsigned curS = as<Rcpp::IntegerVector>( as<Rcpp::List>( S(cidx) ) (r) ) (i) ; //        cur.S = S[[c.idx]] [[r]] [i];
	cur_S_lst.push_back(curS);
	cur_S_count[ curS ] = cur_S_count[ curS ] +1;
	N_S_c_l( curS, curY ) = N_S_c_l( curS, curY ) +1; //        N.S.c.l[ cur.S, cur.Y ] = N.S.c.l[ cur.S, cur.Y ] +1;
	N_S_l(curY) = N_S_l(curY) + 1; //        N.S.l[cur.Y] = N.S.l[cur.Y] + 1;

	unsigned curR = as<Rcpp::IntegerVector>( as<Rcpp::List>( R(cidx) ) (r) ) (i) ; //        cur.R = R[[c.idx]] [[r]] [i];
	cur_R_lst.push_back(curR);
	cur_R_count[ curR ] = cur_R_count[ curR ] +1;
	N_R_c_l( curR, curY ) = N_R_c_l( curR, curY ) +1; //        N.R.c.l[ cur.R, cur.Y ] = N.R.c.l[ cur.R, cur.Y ] +1;
	N_R_l(curY) = N_R_l(curY) +1; //        N.R.l[cur.Y] = N.R.l[cur.Y] +1;
      } //      }
      cur_S_lst.sort(); cur_S_lst.unique(); 
      cur_R_lst.sort(); cur_R_lst.unique();

      std::vector<int> cur_S_count_lst;
      std::vector<int> cur_R_count_lst;
      for (list<int>::iterator n_it=cur_S_lst.begin(); n_it != cur_S_lst.end(); n_it++){
	unsigned int n = *n_it;
	cur_S_count_lst.push_back( cur_S_count[n] );
      }
      for (list<int>::iterator n_it=cur_R_lst.begin(); n_it != cur_R_lst.end(); n_it++){
	unsigned int n = *n_it;
	cur_R_count_lst.push_back( cur_R_count[n] );
      }
      cur_G_S_count.push_back(cur_S_count_lst);
      cur_G_R_count.push_back(cur_R_count_lst);
      cur_G_S_lst.push_back(cur_S_lst);
      cur_G_R_lst.push_back(cur_R_lst);

    } //    }
    S_count.push_back(cur_G_S_count);
    R_count.push_back(cur_G_R_count);
    S_lst.push_back(cur_G_S_lst);
    R_lst.push_back(cur_G_R_lst);

    sum_C_g_k_row(cidx)  = sum ( C_g_k.row(cidx) );
    
  } //  }
  for (unsigned kidx=0; kidx < K; kidx++)
    sum_C_k_l_row(kidx) = sum( C_k_l.row(kidx) );
  

  cout << "Starting Gibbs Iterations" << endl; //  # start Gibbs
  unsigned ndraws = (max_iter - burnin) / iter_lag + 1;
  Rcpp::List allZ(ndraws); //Rcpp::List allZ;
  Rcpp::List allY(ndraws); //Rcpp::List allY;
  for ( unsigned itr=1; itr <= max_iter; itr++ ){ //  for ( itr in 1:max.iter){
    //if (itr % 10 == 0){ //    if (itr %% 5 == 0){
      cout << "iteration " << itr << endl; //      cat(paste("iteration #", itr, "\n"));
      //} //    }

    for ( unsigned cidx=0; cidx < G; cidx++) { //    for (c.idx in 1:G){
      //cout << "sampling for network " << cidx << endl;
      for ( unsigned r=0; r< reaction_count(cidx); r++ ){ //      for ( r in 1:(reaction.count[[c.idx]]) ){
	int curZ = as<Rcpp::IntegerVector>(Z(cidx))(r) ; //        cur.Z = Z[[c.idx]][r] ;
	int curY = as<Rcpp::IntegerVector>(Y(cidx))(r) ; //        cur.Y = Y[[c.idx]][r] ;
	
	C_g_k(cidx, curZ) = C_g_k(cidx, curZ) -1; //        C.g.k[c.idx, cur.Z] = C.g.k[c.idx, cur.Z] -1;
	sum_C_g_k_row(cidx) = sum_C_g_k_row(cidx) -1;

	C_k_l(curZ, curY) = C_k_l(curZ, curY) -1; //        C.k.l[cur.Z, cur.Y] = C.k.l[cur.Z, cur.Y] -1;
	sum_C_k_l_row(curZ) = sum_C_k_l_row(curZ) -1;

	unsigned cur_interaction_count = interaction_count[cidx][r]; //        cur.interaction.count = interaction.count[[c.idx]][r];
	for ( unsigned i=0; i<cur_interaction_count; i++ ){ //        for ( i in 1:cur.interaction.count ){
	  unsigned curS = as<Rcpp::IntegerVector>( as<Rcpp::List>( S(cidx) ) (r) ) (i) ; //          cur.S = S[[c.idx]] [[r]] [i];
	  N_S_c_l( curS, curY ) = N_S_c_l( curS, curY ) -1; //          N.S.c.l[ cur.S, cur.Y ] = N.S.c.l[ cur.S, cur.Y ] -1;
	  N_S_l(curY) = N_S_l(curY) - 1; //          N.S.l[cur.Y] = N.S.l[cur.Y] - 1;

	  unsigned curR = as<Rcpp::IntegerVector>( as<Rcpp::List>( R(cidx) ) (r) ) (i) ; //          cur.R = R[[c.idx]] [[r]] [i];
	  N_R_c_l( curR, curY ) = N_R_c_l( curR, curY ) -1; //          N.R.c.l[ cur.R, cur.Y ] = N.R.c.l[ cur.R, cur.Y ] -1;
	  N_R_l(curY) = N_R_l(curY) - 1; //          N.R.l[cur.Y] = N.R.l[cur.Y] - 1;
	} //        }
		
	std::vector<double> prob_vec(K*L, 1.0/(K*L)); //Rcpp::NumericVector prob_vec = rep(1.0/(K*L), K*L); //        prob.vec = rep(0, K*L);
	
	for (unsigned k=0; k<K; k++){ //        for (k in 1:K){
	  for (unsigned l=0; l<L; l++){ //          for(l in 1:L){
	    //prob_vec(k*L + l) = ( alpha_theta(k) + C_g_k(cidx, k) ) / ( sum( alpha_theta ) + sum ( C_g_k.row(cidx) ) ) * ( alpha_phi(l) + C_k_l(k,l) ) / ( sum( alpha_phi ) + sum( C_k_l.row(k) ) ); //            prob.vec[(k-1)*L + l] = ( alpha.theta[k] + C.g.k[c.idx, k] ) / ( sum( alpha.theta + C.g.k[c.idx, ] ) ) * ( alpha.phi[l] + C.k.l[k,l] ) / ( sum( alpha.phi + C.k.l[k,] ) );
	    //prob_vec(k*L + l) = ( alpha_theta(k) + C_g_k(cidx, k) ) / ( sum_alpha_theta  + sum_C_g_k_row(cidx)  ) * ( alpha_phi(l) + C_k_l(k,l) ) / ( sum_alpha_phi  + sum_C_k_l_row(k)  ); 
	    prob_vec[k*L + l] = ( alpha_theta(k) + C_g_k(cidx, k) ) / ( sum_alpha_theta  + sum_C_g_k_row(cidx)  ) * ( alpha_phi(l) + C_k_l(k,l) ) / ( sum_alpha_phi  + sum_C_k_l_row(k)  ); 

	    unsigned iiii=0;
	    for (list<int>::iterator n_it=S_lst[cidx][r].begin(); n_it != S_lst[cidx][r].end(); n_it++){ //            for ( n in cur.S.lst){
	      unsigned n = *n_it;
	      for (unsigned c=0; c < S_count[cidx][r][iiii]; c++){ //for (unsigned c=0; c < cur_S_count(n); c++){ //              for ( c in 1:cur.S.count[ n ] ){
		//prob_vec(k*L + l) = prob_vec(k*L + l) *  ( alpha_delta(n) + N_S_c_l(n, l) + c ); //                prob.vec[(k-1)*L + l] = prob.vec[(k-1)*L + l] * ( alpha.delta[n] + N.S.c.l[n, l] + c -1);
		prob_vec[k*L + l] = prob_vec[k*L + l] *  ( alpha_delta(n) + N_S_c_l(n, l) + c ); 
	      } //              }
	      iiii++;
	    } //            }
	    for (unsigned c=0; c< cur_interaction_count; c++){ //            for ( n in 1:cur.interaction.count ){ # check this n or c as index
	      //prob_vec(k*L + l) = prob_vec(k*L + l) / ( sum_alpha_delta + N_S_l(l) +c ); //              prob.vec[(k-1)*L + l] = prob.vec[(k-1)*L + l] / ( sum(alpha.delta) + N.S.l[l] +c -1 )
	      prob_vec[k*L + l] = prob_vec[k*L + l] / ( sum_alpha_delta + N_S_l(l) +c );
	    } //            }

	    iiii=0;
	    for (list<int>::iterator n_it=R_lst[cidx][r].begin(); n_it != R_lst[cidx][r].end(); n_it++){ //            for ( n in cur.R.lst){
	      unsigned n = *n_it;
	      for (unsigned c=0; c < R_count[cidx][r][iiii]; c++){ //for (unsigned c=0; c < cur_R_count(n); c++){ //              for ( c in 1:cur.R.count[ n ] ){
		//prob_vec(k*L + l) = prob_vec(k*L + l) * ( alpha_gamma(n) + N_R_c_l(n, l) +c ); //                prob.vec[(k-1)*L + l] = prob.vec[(k-1)*L + l] * ( alpha.gamma[n] + N.R.c.l[n, l] +c -1);
		prob_vec[k*L + l] = prob_vec[k*L + l] * ( alpha_gamma(n) + N_R_c_l(n, l) +c );
	      } //              }
	      iiii++;
	    } //            }
	    for (unsigned c=0; c< cur_interaction_count; c++){ //            for ( n in 1:cur.interaction.count ){
	      //prob_vec(k*L + l) = prob_vec(k*L + l) / ( sum_alpha_gamma + N_R_l(l) +c ); //              prob.vec[(k-1)*L + l] = prob.vec[(k-1)*L + l] / ( sum(alpha.gamma) + N.R.l[l] +c -1)
	      prob_vec[k*L + l] = prob_vec[k*L + l] / ( sum_alpha_gamma + N_R_l(l) +c );
	    } //            }
	    //          
	  } //          }
	} //        }
	//
	unsigned ZY = roll_weighted_die( prob_vec ); //unsigned ZY = Cwhich(rmultinom(1,1, prob_vec)); //        ZY = which(rmultinom(1,1, prob.vec)==1);
	curY = ZY % L; //        cur.Y = ZY %% L;
	curZ = (ZY - curY)/L; //        cur.Z = (ZY - cur.Y)/L;
        
	as<Rcpp::IntegerVector>(Z(cidx))(r) = curZ; //        Z[[c.idx]][r] = cur.Z;
	as<Rcpp::IntegerVector>(Y(cidx))(r) = curY; //        Y[[c.idx]][r] = cur.Y;

	//savedGibbsIterations(drawcount-1) =  Rcpp::List::create(Rcpp::Named("Z")=clone(Z), Rcpp::Named("X")=clone(X) ) ;
	
	//        
	for ( unsigned i=0; i< cur_interaction_count; i++){//        for ( i in 1:cur.interaction.count ){
	  unsigned curS = Rcpp::as<int>( as<Rcpp::List>( as<Rcpp::List>( S(cidx) ) (r) ) (i) ); //          cur.S = S[[c.idx]] [[r]] [i];
	  N_S_c_l( curS, curY ) = N_S_c_l( curS, curY ) +1; //          N.S.c.l[ cur.S, cur.Y ] = N.S.c.l[ cur.S, cur.Y ] +1;
	  N_S_l(curY) = N_S_l(curY) + 1; //          N.S.l[cur.Y] = N.S.l[cur.Y] + 1;
	  //
	  unsigned curR = Rcpp::as<int>( as<Rcpp::List>( as<Rcpp::List>( R(cidx) ) (r) ) (i) ); //          cur.R = R[[c.idx]] [[r]] [i];
	  N_R_c_l( curR, curY ) = N_R_c_l( curR, curY ) +1; //          N.R.c.l[ cur.R, cur.Y ] = N.R.c.l[ cur.R, cur.Y ] +1;
	  N_R_l(curY) = N_R_l(curY) + 1; //          N.R.l[cur.Y] = N.R.l[cur.Y] + 1;
	} //        }
	
	C_g_k(cidx, curZ) = C_g_k(cidx, curZ) +1; //        C.g.k[c.idx, cur.Z] = C.g.k[c.idx, cur.Z] +1;
	sum_C_g_k_row(cidx) = sum_C_g_k_row(cidx) +1;
	C_k_l(curZ, curY) = C_k_l(curZ, curY) +1; //        C.k.l[cur.Z, cur.Y] = C.k.l[cur.Z, cur.Y] +1;
	sum_C_k_l_row(curZ) = sum_C_k_l_row(curZ) +1;
	//      
	//      
      } //      }
    } //    }
    //cout << "before inserting new sample in the list" << endl;
    if ( (itr >= burnin) && ( (itr - burnin) % iter_lag == 0)){
      allZ((itr-burnin)/iter_lag) = clone(Z); //allZ.insert( (itr-burnin)/iter_lag, Z);
      allY((itr-burnin)/iter_lag) = clone(Y); //allY.insert( (itr-burnin)/iter_lag, Y);
    }
    //cout << "after inserting new sample in the list" << endl;
    //
  } //  }
  return Rcpp::List::create(Rcpp::Named("Z") = allZ, Rcpp::Named("Y") = allY); //  list(Zs = Z, Ys = Y);
} //}
