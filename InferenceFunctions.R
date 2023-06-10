# Gibbs sampling implementation

doGibbs.R <- function(S, R, Nodes, K, L, max.iter, Y=NULL, Z=NULL){

  G = length(S);
  N = length(Nodes);
  
  alpha.theta = rep(0.1, K);
  alpha.phi = rep(0.1, L);
  alpha.delta = rep(0.1, N);
  alpha.gamma = rep(0.1, N);
 
  # random starting values for Y and Z if not given
  reaction.count = list();
  interaction.count = list();
  if (is.null(Y)){
    Z = list();
    Y = list();
    for (c.idx in 1:G){
      reaction.count[[c.idx]] = length(S[[c.idx]]);
      interaction.count[[c.idx]] = rep(0, reaction.count[[c.idx]] );
      Z[[c.idx]] = rep(0, reaction.count[[c.idx]] );
      Y[[c.idx]] = rep(0, reaction.count[[c.idx]] );  
      for ( r in 1:(reaction.count[[c.idx]]) ){
        Z[[c.idx]][r] = which(rmultinom(1,1, rep(1,K)/K )==1);  #which(rmultinom(1,1, theta[[c.idx]] )==1);
        Y[[c.idx]][r] = which(rmultinom(1,1, rep(1,L)/L )==1);  #which(rmultinom(1,1, phi.mat[ Z[[c.idx]][r], ] )==1);
      }
    }
  }


  
  # initializing counter variable;
  C.g.k = matrix(0, nrow=G, ncol=K);
  C.k.l = matrix(0, nrow=K, ncol=L);
  N.S.c.l = matrix(0, nrow=N, ncol=L);
  N.R.c.l = matrix(0, nrow=N, ncol=L);
  N.S.l = rep(0, L);
  N.R.l = rep(0, L);

  for (c.idx in 1:G){
    reaction.count[[c.idx]] = length(S[[c.idx]]);
    interaction.count[[c.idx]] = rep(0, reaction.count[[c.idx]] );  
    for ( r in 1:(reaction.count[[c.idx]]) ){
      cur.Z = Z[[c.idx]][r] ;
      cur.Y = Y[[c.idx]][r] ;
      interaction.count[[c.idx]][r] = length(S[[c.idx]][[r]]);
      C.g.k[c.idx, cur.Z] = C.g.k[c.idx, cur.Z] +1;
      C.k.l[cur.Z, cur.Y] = C.k.l[cur.Z, cur.Y] +1;
      for ( i in 1:(interaction.count[[c.idx]][r]) ){
        cur.S = S[[c.idx]] [[r]] [i];
        N.S.c.l[ cur.S, cur.Y ] = N.S.c.l[ cur.S, cur.Y ] +1;
        N.S.l[cur.Y] = N.S.l[cur.Y] + 1;
        cur.R = R[[c.idx]] [[r]] [i];
        N.R.c.l[ cur.R, cur.Y ] = N.R.c.l[ cur.R, cur.Y ] +1;
        N.R.l[cur.Y] = N.R.l[cur.Y] +1;
      }
    }
  }

  # start Gibbs
  for ( itr in 1:max.iter){
    if (itr %% 5 == 0){
      cat(paste("iteration #", itr, "\n"));
    }

    for (c.idx in 1:G){

      reaction.count[[c.idx]] = length(S[[c.idx]]);
      interaction.count[[c.idx]] = rep(0, reaction.count[[c.idx]] );
      #Z[[c.idx]] = rep(0, reaction.count[[c.idx]] );
      #Y[[c.idx]] = rep(0, reaction.count[[c.idx]] );
      
      for ( r in 1:(reaction.count[[c.idx]]) ){
        cur.Z = Z[[c.idx]][r] ;
        cur.Y = Y[[c.idx]][r] ;
        interaction.count[[c.idx]][r] = length(S[[c.idx]][[r]]);


        C.g.k[c.idx, cur.Z] = C.g.k[c.idx, cur.Z] -1;
        C.k.l[cur.Z, cur.Y] = C.k.l[cur.Z, cur.Y] -1;

        cur.S.count = rep(0, N);
        cur.R.count = rep(0, N);
        cur.interaction.count = interaction.count[[c.idx]][r];
        for ( i in 1:cur.interaction.count ){
          cur.S = S[[c.idx]] [[r]] [i];
          cur.S.count[ cur.S ] = cur.S.count[ cur.S ] +1;
          N.S.c.l[ cur.S, cur.Y ] = N.S.c.l[ cur.S, cur.Y ] -1;
          N.S.l[cur.Y] = N.S.l[cur.Y] - 1;

          cur.R = R[[c.idx]] [[r]] [i];
          cur.R.count[ cur.R ] = cur.R.count[ cur.R ] +1;
          N.R.c.l[ cur.R, cur.Y ] = N.R.c.l[ cur.R, cur.Y ] -1;
          N.R.l[cur.Y] = N.R.l[cur.Y] - 1;
        }
        cur.S.lst = unique(S[[c.idx]] [[r]]);
        cur.R.lst = unique(R[[c.idx]] [[r]]);
        prob.vec = rep(0, K*L);
        for (k in 1:K){
          for(l in 1:L){
            prob.vec[(k-1)*L + l] = ( alpha.theta[k] + C.g.k[c.idx, k] ) / ( sum( alpha.theta + C.g.k[c.idx, ] ) ) * ( alpha.phi[l] + C.k.l[k,l] ) / ( sum( alpha.phi + C.k.l[k,] ) );

            for ( n in cur.S.lst){
              for ( c in 1:cur.S.count[ n ] ){
                prob.vec[(k-1)*L + l] = prob.vec[(k-1)*L + l] * ( alpha.delta[n] + N.S.c.l[n, l] + c -1);
              }
            }
            for ( n in 1:cur.interaction.count ){
              prob.vec[(k-1)*L + l] = prob.vec[(k-1)*L + l] / ( sum(alpha.delta) + N.S.l[l] +c -1 )
            }
            for ( n in cur.R.lst){
              for ( c in 1:cur.R.count[ n ] ){
                prob.vec[(k-1)*L + l] = prob.vec[(k-1)*L + l] * ( alpha.gamma[n] + N.R.c.l[n, l] +c -1);
              }
            }
            for ( n in 1:cur.interaction.count ){
              prob.vec[(k-1)*L + l] = prob.vec[(k-1)*L + l] / ( sum(alpha.gamma) + N.R.l[l] +c -1)
            }
          
          }
        }

        ZY = which(rmultinom(1,1, prob.vec)==1);
        cur.Y = ZY %% L;
        cur.Z = (ZY - cur.Y)/L;
        if (cur.Y == 0){
          cur.Y = L;
        } else{
          cur.Z = cur.Z +1;
        }
        
        Z[[c.idx]][r] = cur.Z;
        Y[[c.idx]][r] = cur.Y;
        C.g.k[c.idx, cur.Z] = C.g.k[c.idx, cur.Z] +1;
        C.k.l[cur.Z, cur.Y] = C.k.l[cur.Z, cur.Y] +1;
        
        for ( i in 1:cur.interaction.count ){
          cur.S = S[[c.idx]] [[r]] [i];
          N.S.c.l[ cur.S, cur.Y ] = N.S.c.l[ cur.S, cur.Y ] +1;
          N.S.l[cur.Y] = N.S.l[cur.Y] + 1;

          cur.R = R[[c.idx]] [[r]] [i];
          N.R.c.l[ cur.R, cur.Y ] = N.R.c.l[ cur.R, cur.Y ] +1;
          N.R.l[cur.Y] = N.R.l[cur.Y] + 1;
        }
      
      
      }
    }

  }
  list(Zs = Z, Ys = Y);
}

estimate.theta <- function(Z.samples, K, alpha.theta, itr){
  G = length(Z.samples[[1]]);
  if (length(alpha.theta) == 1)
    alpha.theta = rep(alpha.theta, K);
  
  est.theta = matrix(0, nrow=G, ncol=K)
  for (i in 1:G) {
    z.samples.G = factor(Z.samples[[itr]][[i]], levels=seq(0, (K-1)) )
    est.theta[i,] = (table(z.samples.G) + alpha.theta) / ( length(Z.samples[[itr]][[i]]) + sum(alpha.theta) );
    #est.theta[i,] = (table(Z.samples[[itr]][[i]]) + alpha.theta) / ( length(Z.samples[[itr]][[i]]) + sum(alpha.theta) );
  }

  est.theta
}

estimate.phi <- function(Z.samples, Y.samples, K, L, alpha.phi, itr){
  G = length(Z.samples[[1]]);
  if (length(alpha.phi) == 1)
    alpha.phi = rep(alpha.phi, L);
  
  est.phi = matrix(0, nrow=K, ncol=L)
  for (i in 1:G) { 
    z.samples.G = factor(Z.samples[[itr]][[i]], levels=seq(0, (K-1)) )
    y.samples.G = factor(Y.samples[[itr]][[i]], levels=seq(0, (L-1)) )
    est.phi = est.phi + table( z.samples.G, y.samples.G ) ; 
    #est.phi = est.phi + table( Z.samples[[itr]][[i]], Y.samples[[itr]][[i]] ) ; 
  }

  for (k in 1:K) { est.phi[k,] = ( est.phi[k,] + alpha.phi) / (sum(est.phi[k,]) + sum(alpha.phi) ); }

  est.phi
}

estimate.delta <- function(Y.samples, S, Nodes, L, alpha.delta, itr){
  Y = Y.samples[[itr]];
  G = length(S);
  N = length(Nodes);
  if (length(alpha.delta) == 1)
    alpha.delta = rep(0.1, N);

  est.delta = matrix(0, nrow=L, ncol=N)
  colnames(est.delta) = names(Nodes)

  for (c.idx in 1:G){
    for ( r in 1:(length(S[[c.idx]])) ){
      cur.Y = Y[[c.idx]][r] +1;
      for ( i in 1:(length(S[[c.idx]][[r]])) ){
        cur.S = S[[c.idx]] [[r]] [i] + 1
        est.delta[cur.Y, cur.S] = est.delta[cur.Y, cur.S] +1;
      }
    }
  }

  for (l in 1:L){ est.delta[l,] = ( est.delta[l,] + alpha.delta ) / ( sum(est.delta[l,]) + sum(alpha.delta) ); }

  est.delta

  
}

estimate.gamma <- function(Y.samples, R, Nodes, L, alpha.gamma, itr){
  Y = Y.samples[[itr]];
  G = length(R);
  N = length(Nodes);
  if (length(alpha.gamma) == 1)
    alpha.gamma = rep(0.1, N);

  est.gamma = matrix(0, nrow=L, ncol=N)
  colnames(est.gamma) = names(Nodes)

  for (c.idx in 1:G){
    for ( r in 1:(length(R[[c.idx]])) ){
      cur.Y = Y[[c.idx]][r] + 1;
      for ( i in 1:(length(R[[c.idx]][[r]])) ){
        cur.R = R[[c.idx]] [[r]] [i] + 1
        est.gamma[cur.Y, cur.R] = est.gamma[cur.Y, cur.R] +1;
      }
    }
  }

  for (l in 1:L){ est.gamma[l,] = ( est.gamma[l,] + alpha.gamma ) / ( sum(est.gamma[l,]) + sum(alpha.gamma) ); }

  est.gamma

  
}

estimate.reaction.membership <- function(Y.samples, Reactions, L, itr){
  Y = Y.samples[[itr]];
  G = length(Y);

  Y.unlisted = unlist(Y) 
  Reactions.unlisted = unlist(Reactions)

  reaction.set = unique(Reactions.unlisted)

  reaction.set.membership = table(Reactions.unlisted, Y.unlisted ) 
  
  #reaction.set.membership = matrix(0, nrow=length(reaction.set), ncol=L)
  #rownames(reaction.set.membership) = reaction.set

  
  #for ( i in 1:length(Reactions.unlisted) ){
  #  reaction.set.membership[ Reactions.unlisted[i], Y.unlisted ] = reaction.set.membership[ Reactions.unlisted[i], Y.unlisted ] +1
  #}
  
  #for (c.idx in 1:G){
  #  for ( r in 1:(length(Reactions[[c.idx]])) ){
  #    cur.Y = Y[[c.idx]][r] + 1;
  #    reaction.set.membership[ Reactions[[c.idx]][[r]], cur.Y] = reaction.set.membership[ Reactions[[c.idx]][[r]], cur.Y] +1
  #  }
  #}
  
  reaction.set.membership
}

estimate.reaction.membership.sample.group <- function(reaction.membership, est.phi) {
  Reactions.count = nrow(reaction.membership)
  L = ncol(reaction.membership)
  K = nrow(est.phi)
  
  reaction.membership.sample.group = reaction.membership %*% t(est.phi)

  reaction.membership.sample.group
  
  
}
