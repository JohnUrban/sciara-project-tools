

puffRstring = '''
impute_zeros <- function(x, y, bw){
    k <- ksmooth(x=x, y=y, bandwidth=bw)
    y[y == 0] <- k$y[y == 0]
    return(y)
}
mednorm <- function(x){x/median(x)}
mednorm.ksmooth <-function(x,y,bw){mednorm(ksmooth(x=x,y=y,bandwidth = bw)$y)}
mednorm.ksmooth.norm <-function(x,y,bw,norm.y){mednorm.ksmooth(x,y,bw)/mednorm.ksmooth(x,norm.y,bw)}
inner.quant.mean.norm <- function(x, inner=c(0.4,0.6)){
  innerq <- quantile(x=x, probs = inner)
  x/mean(x[x >= innerq[1] & x <= innerq[2]])
}

slope <- function(x1,x2,y1,y2){
  (y2-y1)/(x2-x1)
}



prob.data <- function(Forward){
  num.emits <- dim(Forward)[2]
  return(sum(Forward[,num.emits]))
}




# HMM Functions
viterbi.puff <- function(emissions, transitions, initial, states, emitted.data, emodel="normal", logprobs=FALSE){
  ## logprobs = whether or not transition and initial matrices are logged yet, assumes false
  if(!(logprobs)){
    initial = log(initial)
    transitions = log(transitions)
  }
  num.states <- length(states)
  num.emits <- length(emitted.data) 
  ## need to add log_probs instead of multiply probs to prevent underflow
##  write.table(emissions, file="del1")
##  write.table(transitions, file="del2")
##  write.table(initial,file="del3")
##  write.table(num.states, file="del4")
##  write.table(num.emits, file="del5")
  if (emodel == "normal"){V=viterbi.normal(emissions, transitions, initial, emitted.data, num.states, num.emits)}
  else if (emodel == "poisson"){V=viterbi.poisson(emissions, transitions, initial, emitted.data, num.states, num.emits)}
  ## NOTE: for following two emission_prob_means that would otherwise be used for normal or poissoin are inversed (1/mu) for exponential and geometric
  else if (emodel == "exponential"){
      emissions[1,] <- 1/emissions[1,]
      V=viterbi.exponential(emissions, transitions, initial, emitted.data, num.states, num.emits)
  }
  else if (emodel == "geometric"){
      emissions[1,] <- 1/emissions[1,]
      V=viterbi.geometric(emissions, transitions, initial, emitted.data, num.states, num.emits)
  } else if (emodel == "gamma"){
      params <- emissions
      ## Estimate shape
      params[1,] <- emissions[1,]^2 / emissions[2,]^2
      ## Estimate scale
      params[2,] <- emissions[2,]^2 / emissions[1,]
      ## Stay calm and carry on
      V=viterbi.gamma(params, transitions, initial, emitted.data, num.states, num.emits)
  }
  
  viterbi_path <- matrix(data = rep(0, num.emits), nrow = 1)
  viterbi_path[num.emits] <- which.max(V$Viterbi[,num.emits]);  
  viterbi_prob <- V$Viterbi[viterbi_path[num.emits], num.emits]
  
  for (j in num.emits:2){
    viterbi_path[j-1] <- V$pointer[j,viterbi_path[j]]
  }

  return(list(viterbi_path=viterbi_path, viterbi_prob=viterbi_prob))
}

##NORMAL
viterbi.normal <- function(emissions, transitions, initial, emitted.data, num.states, num.emits){
  pointer <- matrix(rep(0, num.emits*num.states), nrow = num.emits)
  Viterbi <- matrix(rep(0, num.states*num.emits), nrow = num.states)  
  Viterbi[ ,1] <- initial + dnorm(emitted.data[1], mean = emissions[1, ], sd = emissions[2, ], log = TRUE)
  pointer[1, ] <- 1
  f <- function(x){i <- which.max(x); y <- x[i]; return(c(i,y))}
  for (j in 2:num.emits){
    selection <- Viterbi[,j-1] + transitions
    for (i in 1:num.states){
      maxstate <- which.max(selection[,i])
      Viterbi[i,j] <- dnorm(emitted.data[j], mean = emissions[1, i], sd = emissions[2, i], log = TRUE) + selection[maxstate,i]
      pointer[j,i] <- maxstate 
    }
  }   
  return(list(Viterbi=Viterbi, pointer=pointer))
}

##POISSON
viterbi.poisson <- function(emissions, transitions, initial, emitted.data, num.states, num.emits){
  ## rounds floats to integers
  pointer <- matrix(rep(0, num.emits*num.states), nrow = num.emits)
  Viterbi <- matrix(rep(0, num.states*num.emits), nrow = num.states)  
  Viterbi[ ,1] <- initial + dpois(round(emitted.data[1]), lambda = emissions[1, ], log=TRUE)
  pointer[1, ] <- 1
  f <- function(x){i <- which.max(x); y <- x[i]; return(c(i,y))}
  for (j in 2:num.emits){
    selection <- Viterbi[,j-1] + transitions
    for (i in 1:num.states){
      maxstate <- which.max(selection[,i])
      Viterbi[i,j] <- dpois(round(emitted.data[j]), lambda = emissions[1, i], log = TRUE) + selection[maxstate,i]
      pointer[j,i] <- maxstate 
    }
  }  
  return(list(Viterbi=Viterbi, pointer=pointer))
}

## EXPONENTIAL
viterbi.exponential <- function(emissions, transitions, initial, emitted.data, num.states, num.emits){
  pointer <- matrix(rep(0, num.emits*num.states), nrow = num.emits)
  Viterbi <- matrix(rep(0, num.states*num.emits), nrow = num.states)  
  Viterbi[ ,1] <- initial + dexp(emitted.data[1], rate = emissions[1, ], log=TRUE)
  pointer[1, ] <- 1
  f <- function(x){i <- which.max(x); y <- x[i]; return(c(i,y))}
  for (j in 2:num.emits){
    selection <- Viterbi[,j-1] + transitions
    for (i in 1:num.states){
      maxstate <- which.max(selection[,i])
      Viterbi[i,j] <- dexp(emitted.data[j], rate = emissions[1, i], log = TRUE) + selection[maxstate,i]
      pointer[j,i] <- maxstate 
    }
  }  
  return(list(Viterbi=Viterbi, pointer=pointer))
}


## GEOMETRIC
viterbi.geometric <- function(emissions, transitions, initial, emitted.data, num.states, num.emits){
  ## rounds floats to integers
  pointer <- matrix(rep(0, num.emits*num.states), nrow = num.emits)
  Viterbi <- matrix(rep(0, num.states*num.emits), nrow = num.states)  
  Viterbi[ ,1] <- initial + dgeom(round(emitted.data[1]), prob = emissions[1, ], log=TRUE)
  pointer[1, ] <- 1
  f <- function(x){i <- which.max(x); y <- x[i]; return(c(i,y))}
  for (j in 2:num.emits){
    selection <- Viterbi[,j-1] + transitions
    for (i in 1:num.states){
      maxstate <- which.max(selection[,i])
      Viterbi[i,j] <- dgeom(round(emitted.data[j]), prob = emissions[1, i], log = TRUE) + selection[maxstate,i]
      pointer[j,i] <- maxstate 
    }
  }  
  return(list(Viterbi=Viterbi, pointer=pointer))
}

## GAMMA
viterbi.gamma <- function(emissions, transitions, initial, emitted.data, num.states, num.emits){
  pointer <- matrix(rep(0, num.emits*num.states), nrow = num.emits)
  Viterbi <- matrix(rep(0, num.states*num.emits), nrow = num.states)  
  Viterbi[ ,1] <- initial + dgamma(emitted.data[1], shape = emissions[1, ], scale =emissions[2, ], log=TRUE)
  pointer[1, ] <- 1
  f <- function(x){i <- which.max(x); y <- x[i]; return(c(i,y))}
  for (j in 2:num.emits){
    selection <- Viterbi[,j-1] + transitions
    for (i in 1:num.states){
      maxstate <- which.max(selection[,i])
      Viterbi[i,j] <- dgamma(emitted.data[j], shape = emissions[1, i], scale = emissions[2, i], log = TRUE) + selection[maxstate,i]
      pointer[j,i] <- maxstate 
    }
  }  
  return(list(Viterbi=Viterbi, pointer=pointer))
}


#### VECTORIZED FORWARDS ########
forward.puff <- function(emissions, transitions, initial, states, emitted.data, emodel="normal"){
  num.states <- length(states)
  num.emits <- length(emitted.data)
  Forward <- matrix(data = rep(0, num.states*num.emits), nrow = num.states)
  scalefactors <- matrix(data = rep(0, num.emits*2), nrow = 2)
  
  #model
  if (emodel == "normal"){emodel.fxn <- puff.normal}
  else if (emodel == "exponential"){
      emissions[1,] <- 1/emissions[1,]
      emodel.fxn <- puff.exponential
  } else if (emodel == "poisson"){emodel.fxn <- puff.poisson}
  else if (emodel == "geometric"){
      emissions[1,] <- 1/emissions[1,]
      emodel.fxn <- puff.geometric
  } else if (emodel == "gamma") {
      params <- emissions
      ## Estimate shape
      params[1,] <- emissions[1,]^2 / emissions[2,]^2
      ## Estimate scale
      params[2,] <- emissions[2,]^2 / emissions[1,]
      ## Stay calm and carry on
      emissions <- params
      emodel.fxn <- puff.gamma
  }
  
  ## initial
  Forward[, 1] <- initial*emodel.fxn(emitted.data[1], emissions)
  ## scale to prevent underflow -- keep track of scaling
  scalefactors[1,1] <- sum(Forward[, 1])
  scalefactors[2,1] <- log(scalefactors[1,1])
  Forward[,1] <- Forward[,1]/scalefactors[1,1]
  
  ## iterate
  for(k in 2:num.emits){
    emit <- emodel.fxn(emitted.data[k], emissions)
    Forward[, k] <- emit* Forward[,k-1] %*% transitions ## same as emit* Forward[,k-1] * colSums(transitions)
    scalefactors[1,k] <- sum(Forward[, k])
    scalefactors[2,k] <- log(scalefactors[1,k]) + scalefactors[2,k-1]
    Forward[,k] <- Forward[,k]/scalefactors[1,k]
  }
  
  return(list(forward=Forward, scales=scalefactors))
  ## mutiply forward column by row2,samecol in scale factors OR by exp(row3,samecol in scalfactors) to get actual value for forward
  ## update: OR actually I think multiply fwd column by product of [row1,1:samecol] in scale factors
  ## I must have at one point taken out row2 of scale factors (when row3 was logs of row2)
}




#### VECTORIZED BACKWARDS ########
backward.puff <- function(emissions, transitions, initial, states, emitted.data, emodel="normal"){
  num.states <- length(states)
  num.emits <- length(emitted.data)
  Backward = matrix(data = rep(0, num.states*num.emits), nrow = num.states)
  scalefactors <- matrix(data = rep(0, num.emits*2), nrow = 2)
  
  #model
  if (emodel == "normal"){emodel.fxn <- puff.normal}
  else if (emodel == "exponential"){
      emissions[1,] <- 1/emissions[1,]
      emodel.fxn <- puff.exponential
  }
  else if (emodel == "poisson"){emodel.fxn <- puff.poisson}
  else if (emodel == "geometric"){
      emissions[1,] <- 1/emissions[1,]
      emodel.fxn <- puff.geometric
  }
  else if (emodel == "gamma"){
      params <- emissions
      ## Estimate shape
      params[1,] <- emissions[1,]^2 / emissions[2,]^2
      ## Estimate scale
      params[2,] <- emissions[2,]^2 / emissions[1,]
      ## Stay calm and carry on
      emissions <- params
      emodel.fxn <- puff.gamma
  }
  
  ## initial
  Backward[ , num.emits] <- 1
  
  ## scale to prevent underflow -- keep track of scaling
  scalefactors[1,num.emits] <- sum(Backward[, num.emits])
  scalefactors[2,num.emits] <- log(scalefactors[1,num.emits])
  Backward[,num.emits] <- Backward[,num.emits]/scalefactors[1,num.emits]
  
  ## iterate
  for(k in (num.emits-1):1){
    #     emit <- matrix(dnorm(emitted.data[k+1], mean = emissions[1, ], sd = emissions[2, ]))
    emit <- matrix(emodel.fxn(emitted.data[k+1], emissions))
    #     print(Backward[, k+1] * emit)
    Backward [, k] <- transitions %*% (Backward[, k+1] * emit)
    scalefactors[1,k] <- sum(Backward[, k])
    scalefactors[2,k] <- log(scalefactors[1,k]) + scalefactors[2,k+1]
    Backward[,k] <- Backward[,k]/scalefactors[1,k]
  }
  return(list(backward=Backward, scales=scalefactors))
}


puff.normal <- function(x, emissions){
  dnorm(x, mean = emissions[1, ], sd = emissions[2, ], log=FALSE)
}

puff.exponential <- function(x, emissions){
  dexp(x = x, rate = emissions[1, ], log = FALSE)
}

puff.poisson <- function(x, emissions){
  dpois(x = round(x), lambda = emissions[1, ], log = FALSE)
}

puff.geometric <- function(x, emissions){
  dgeom(x = round(x), prob = emissions[1, ], log = FALSE)
}


puff.gamma <- function(x, emissions){
  dgamma(x = x, shape = emissions[1, ], scale = emissions[2, ], log=FALSE)
}

###
posterior <- function(Forward, Backward, states){
  ## F and B matrices are from small sequence fxns -- not scaled
  num.states <- length(states)
  num.emits <- dim(Forward)[2]
  posterior.path <- matrix(data = rep(0, num.emits), nrow = 1)
  probs <- matrix(data = rep(0, num.emits), nrow = 1)
  pd <- prob.data(Forward = Forward)
  for (i in 1:num.emits){
    fb <- Forward[,i]*Backward[,i]
    max.state <- which.max(fb)
    posterior.path[i] <- max.state
    #     probs[i] <- max(fb)/pd ## should be divided by prob.data...?
  }
  return(list(posterior.path=posterior.path, probs=probs))  
}

compare.statepath <- function(sp1, sp2){
  # where sp is a numeric vector or char vector with each state its own element -- not seq format
  total <- length(sp1)
  ident <- sum(sp1 == sp2)
  edit.dist <- total - ident
  return(list(edit.dist=edit.dist, identical.count=ident, pct.id=100*ident/total))
}


baum.welch.puff <- function(emissions, transitions, initial, states, emitted.data, emodel="normal"){
  #model
  if (emodel == "normal"){emodel.fxn <- puff.normal}
  else if (emodel == "exponential"){emodel.fxn <- puff.exponential}
  else if (emodel == "poisson"){emodel.fxn <- puff.poisson}
  else if (emodel == "geometric"){emodel.fxn <- puff.geometric}
  # emissions, transitions, and initial probs given are "best guesses" (or randomly chosen)
  # calculate log-likelihood of model
  n.states <- length(states)
  n.emits <- length(emitted.data)
  c <- 0.00001
  new.L <- 0
  old.L <- 100 #arbitrary number producing difference > c
  while (abs(new.L - old.L > c)){
    old.L <- new.L
    # emissions, transitions, and initial probs given are "best guesses" (or randomly chosen)
    # calculate log-likelihood of model
    # Get fwd, backward, and prob(data)
    fwd <- forward.puff(emissions, transitions, initial, states, emitted.data, emodel)
    bck <- backward.puff(emissions, transitions, initial, states, emitted.data, emodel)
    p <- prob.data(Forward = fwd$forward)
    new.L <- log10(p)
    #update initial, transition, emissions
    # calc new log likelihood of model
    # calc difference between new and old log likelihood
    # if diff > cutoff, return to fwd/bck/p step; else done 
    TRANS <- update.transitions(n.states, n.emits, fwd, bck, transitions, emissions)
    EMIS <- update.emissions()
  }
}


# update.transitions <- function(n.states, n.emits, fwd, bck, transitions, emissions){
#   TRANS <- matrix(rep(0, n.states*n.states), nrow=n.states)
#   for (i in 1:n.states){
#     for (k in 1:n.states){
#       for (m in 1:(n.emits-1)){
#         TRANS[i,k] <- TRANS[i,k] + fwd[i,m] * transitions[i,k] * ##TODO#emissions[k,emitted.data[1,m+1]## * bck[k,m+1]
#       }
#     }
#   }
#   return(TRANS)
# }

# update.emissions <- function(n.states, n.emits){
#   EMIS <- matrix(rep(0, n.states*n.emits), nrow=n.states)
#   for (i in 1:n.states){
#     for (k in 1:)
#   }
# }
## if took N=10 backsamples -- ie. N state paths
## could get mean and std dev (emission dist) for each state by simple count stats
## could also get transition probs by simple count stats

backsample.puff <- function(Forward, transitions, n.states=NA, states=NA, n.emits=NA){
  ## TODO vectorize -- e.g. eliminate for loop
  if(is.na(n.states) || is.na(states) || is.na(n.emits)){
      dim.fwd <- dim(Forward)
      n.states <- dim.fwd[1]
      states <- 1:n.states
      n.emits <- dim.fwd[2] 
  }
  #initialization
  b.sample <- rep(0, n.emits)
  # Randomly sample a state for Sn according to: P(Sn=sn|X1:n) = P(Sn=sn,X1:n)/P(X1:n)
  ## p(data) not nec since it scales all by same number and all proportionally the same
  b.sample[n.emits] <- sample(x = states, size = 1, prob = Forward[,n.emits]) #/p)
  
  ## Iterate for k in n.emits-1 to 1
  ## Randomly sample a state for Sk according to: P(Sk=sk|Sk+1=sk+1, X1:n) 
  ## = P(Sk=sk, X1:k)*P(Sk+1=sk+1|Sk=sk)/ [sum(Sk) P(Sk, X1:k)*P(Sk+1=sk+1|Sk)]
  ## = fwd_k(Sk=sk) * trans(Sk+1=sk+1 | Sk=sk) / sum(all states using numerator terms)
  for (k in (n.emits-1):1){
    b.sample[k] <- sample(x = states, size = 1, prob = Forward[,k]*transitions[,b.sample[k+1]]) # no need to re-scale fwd values since they would still be proportionally same
  }
  
 return(b.sample)
}

n.backsamples <- function(Forward, transitions, N=1000){
  #Forward is Forward matrix (not list object with F and scales)
  # transitions is trans prob matrix
  ## TODO vectorize -- e.g. can pick N states for each at once
  dim.fwd <- dim(Forward)
  n.states <- dim.fwd[1]
  states <- 1:n.states
  n.emits <- dim.fwd[2]
  b.samples <- matrix(rep(0, N*n.emits), nrow = N)
  for(i in 1:N){
    b.samples[i,] <- backsample.puff(Forward, transitions, n.states=n.states, states=states, n.emits=n.emits)
  }
  return(b.samples)
}

backsample.state.freq <- function(b.samples, n.states, states, N=NA, n.emits=NA){
  #b.samples is a n.bsamples x n.emits matrix output from n.backsamples
  if(is.na(N) || is.na(n.emits)){
    d <- dim(b.samples)
    N <- d[1]
    n.emits <- d[2]
  }
  freq <- matrix(rep(0, n.states*n.emits), nrow = n.states)
  for(i in 1:n.emits){
    freq[,i] <-  table(c(states,b.samples[,i]))-1 #adding all states in to ensure all levels are represented followed by subtracting 1 from all counts
  }
  return(freq)
}


backsample.max.freq.path <- function(freq){
  apply(X = freq, MARGIN = 2, FUN = which.max)
}


##params for 3state
##initial <- matrix(c(1,1,1)/3, nrow=1)
emissions <- matrix(rep(0, 6), nrow=2)
emissions[1,] <- c(-0.5,0,0.5)
emissions[2,] <- c(0.5,0.5,0.5)


##transitions <- matrix(rep(0,9),nrow=3)
##transitions[1,] <- c(0.99, 0.005, 0.005)
##transitions[2,] <- c(0.005,0.99,0.005)
##transitions[3,] <- c(0.005,0.005,0.99)
#
initial <- matrix(c(0.006,0.988,0.006), nrow=1)
transitions <- matrix(rep(0,9),nrow=3)
transitions[1,] <- c(0.99998, 0.00001, 0.00001)
transitions[2,] <- c(0.000000125,0.99999975,0.000000125)
transitions[3,] <- c(0.00001,0.00001,0.99998)


## params for 7state
##initial7 <- matrix(rep(1,7)/7, nrow=1)
initial7 <- matrix(c(0.997,rep(0.0005,6)), nrow=1)

s <- c(1,sqrt(2),2,sqrt(8),4,sqrt(32),8) 
m <- c(1,2,4,8,16,32,64) 

##For exponential and geometric
###m <- 1/m


emissions7 <- matrix(rep(0, 14), nrow=2)
emissions7[1,] <- m
emissions7[2,] <- s



transitions7 <- matrix(rep(0,49),nrow=7)
for(i in 1:7){
  transitions7[i,1:7] <- 0.001 #0.000001
  transitions7[i,i] <- 0.999 #0.999999
##  transitions7[i,1:7] <- 0.000001
##  transitions7[i,i] <- 0.999999
  #   if(i>1){transitions7[i,i-1] <- 0.000005}
  #   if(i<7){transitions7[i,i+1] <- 0.000005}
  transitions7[i,] <- transitions7[i,]/sum(transitions7[i,])
}

##transitions7 <- matrix(rep(0,49),nrow=7)
##for(i in 1:7){
##  transitions7[i,1:7] <- 0.0001 
##  transitions7[i,i] <- 0.9999 
##  transitions7[i,] <- transitions7[i,]/sum(transitions7[i,])
##}

##transitions7 <- matrix(rep(0,49),nrow=7)
##for(i in 1:7){
##  if(i>1){
##    for (j in 1:(i-1)){transitions7[i,(i-j)] <- 0.0001/j}
##  }
##  if (i<7){
##    for (j in 1:(7-i)){transitions7[i,(i+j)] <- 0.0001/j}
##  }
##  transitions7[i,i] <- 0.9999 
##  transitions7[i,] <- transitions7[i,]/sum(transitions7[i,])
##}

##transitions7 <- matrix(rep(0,49),nrow=7)
##for(i in 1:7){
##  if(i>1){
##    for (j in 1:(i-1)){transitions7[i,(i-j)] <- 0.0001/j^3}
##  }
##  if (i<7){
##    for (j in 1:(7-i)){transitions7[i,(i+j)] <- 0.0001/j^3}
##  }
##  transitions7[i,i] <- 0.9999 
##  transitions7[i,] <- transitions7[i,]/sum(transitions7[i,])
##}

##transitions7 <- matrix(rep(0,49),nrow=7)
##for(i in 1:7){
##  if(i>1){
##    for (j in 1:(i-1)){transitions7[i,(i-j)] <- 0.001/j^3}
##  }
##  if (i<7){
##    for (j in 1:(7-i)){transitions7[i,(i+j)] <- 0.001/j^3}
##  }
##  transitions7[i,i] <- 0.999
##  transitions7[i,] <- transitions7[i,]/sum(transitions7[i,])
##}


generate.normal <- function(n, mu, sig){
  rnorm(n, mean = mu, sd = sig)
}

generate.exponential <- function(n, mu, sig){
  ## assumes mu already 1/mu_given
  ## sig is just dummy var
  rexp(n, rate=mu)
}

generate.poisson <- function(n, mu, sig){
  ## mu is rounded
  rpois(n, lambda = round(mu))
}

generate.geometric <- function(n, mu, sig){
  ## assumes mu is 1/mu_given
  ## sig is just dummy var
  rgeom(n, prob = mu)
}

##generate_statepath <- function(transitions, initial, states, len=10){
##  statenums <- 1:length(states)
##  statepath <- vector(mode="integer", length=len)
##  # INITIAL
##  statepath[1] <- sample(statenums, size = 1, prob = initial)
##    ## TRANSITIONS
##  for (i in 2:len){
##    statepath[i] <- sample(statenums, size=1, prob = transitions[statepath[i-1], ])
##  }
##  return(statepath)
##}
##
##generate_emitted_data <- function(emissions, statepath, emodel = 'normal'){
##  #model
##  if (emodel == "normal"){emodel.fxn <- generate.normal}
##  else if (emodel == "exponential"){emodel.fxn <- generate.exponential}
##  else if (emodel == "poisson"){emodel.fxn <- generate.poisson}
##  else if (emodel == "geometric"){emodel.fxn <- generate.geometric}
##
##  statepathlen = length(statepath)
##  emitted_data <- vector(mode='numeric', length=statepathlen)
##  for (i in 1:statepathlen){
##    emitted_data[i] <- emodel.fxn(n=1, mu=emissions[1, statepath[i]], sig=emissions[2, statepath[i]])
##  }
##  return(emitted_data)
##}


generate <- function(emissions, transitions, initial, states, statepathlen=10, emodel="normal"){
  #model
  if (emodel == "normal"){emodel.fxn <- generate.normal}
  else if (emodel == "exponential"){emodel.fxn <- generate.exponential}
  else if (emodel == "poisson"){emodel.fxn <- generate.poisson}
  else if (emodel == "geometric"){emodel.fxn <- generate.geometric}

  ## Ensure states are indexes
  statenums <- 1:length(states)
  statepath <- vector(mode="integer", length=statepathlen)
  emitted_data <- vector(mode='numeric', length=statepathlen)
  
  # INITIAL
  statepath[1] <- sample(statenums, size = 1, prob = initial)
  emitted_data[1] <- emodel.fxn(n=1, mu=emissions[1, statepath[1]], sig=emissions[2, statepath[1]])
  
    ## TRANSITIONS
  for (i in 2:statepathlen){
    statepath[i] <- sample(statenums, size=1, prob = transitions[statepath[i-1], ])
    emitted_data[i] <- emodel.fxn(n=1, mu=emissions[1, statepath[i]], sig=emissions[2, statepath[i]])
  }
  return(list(statepath, emitted_data))
}
'''

