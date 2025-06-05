library(tidyverse)
library(matlib)

# this is a function that calculates the probability that each player
# is eliminated in a game with n players
lightningelim <- function(SS_make, LS_make, n) {
  p <- LS_make
  q <- SS_make
  Tr <- diag(6*n)
  for (i in 1:(5*n)) Tr[i,i] = 0
  for (i in 1:(n-2)) Tr[i,i+2] = p^2
  Tr[n-1,1] = p^2
  Tr[n,2] = p^2
  for (i in 1:(n-2)) Tr[i+n, i+2] = q^2
  Tr[2*n-1,1] = q^2
  Tr[2*n, 2] = q^2
  for (i in 1:(n-2)) Tr[i+3*n,i+2] = q*p
  Tr[4*n-1,1] = q*p
  Tr[4*n,2] = q*p
  for (i in 1:n) Tr[i,i+n] = (1.0-p)^2
  for (i in 1:n) Tr[i+n,i+n] = (1.0-q)^2
  for (i in 1:n) Tr[i+3*n, i+n] = (1.0-p)*(1.0-q)
  for (i in 1:n) Tr[i+2*n, i+2*n] = (1.0-q)^2
  for (i in 1:n) Tr[i+4*n, i+2*n]= (1.0-p)*(1.0-q)
  for (i in 1:(n-1)) Tr[i+2*n, i+3*n+1] = (1.0-q)*q
  Tr[3*n,3*n+1] = (1.0-q)*q
  for (i in 1:(n-1)) Tr[i+4*n,i+3*n+1] = (1.0-p)*q
  Tr[5*n,3*n+1] = (1.0-p)*q
  for (i in 1:(n-1)) Tr[i, i+4*n+1] = p*(1.0-p)
  Tr[n,4*n+1] = p*(1.0-p)
  for (i in 1:(n-1)) Tr[i+n, i+4*n+1] = q*(1.0-q)
  Tr[2*n,4*n+1] = q*(1.0-q)
  for (i in 1:(n-1)) Tr[i+3*n, i+4*n+1] = q*(1.0-p)
  Tr[4*n,4*n+1] = q*(1.0 - p)
  for (i in (1:n)) Tr[i,i+5*n] = (1.0-p)*p
  for (i in 1:n) Tr[i+n, i+5*n] = (1.0-q)*q
  for (i in 1:n) Tr[i+2*n,i+5*n] = q
  for (i in 1:n) Tr[i+3*n, i+5*n] = (1.0-q)*p
  for (i in 1:n) Tr[i+4*n, i+5*n] = p
  
  Q <- Tr[1:(5*n), 1:(5*n)]
  R <- Tr[1:(5*n), (5*n+1):(6*n)]
  I <- diag(5*n)
#  N <- inv(I-Q)[1,]
#  return(N %*% R)
  x <- solve(I-Q,R)
  return(x[1,])
}

# Calculates each player's probability of winning
lightning <- function(SS_make, LS_make, n) {
 W <- matrix(c(1.0/(3.0-LS_make), 1-1.0/(3.0-LS_make)), nrow = 1, ncol = 2)
 count <- 2
 while (count < n) {
  oldW <- W
  count <- count + 1
  E <- lightningelim(SS_make,LS_make,count)
  W <- matrix(c(1:count), nrow = 1, ncol = count)
  for (j in 1:count) {
   prob <- 0
   for (i in 1:(count-1)) {
    prob <- prob + E[(((j+i-1)%%count) +1)]*oldW[1,(((count-1-1-i) %% (count-1)) +1)]
   }
   W[1,j] = prob
  } 
 }
 return(W)
}

# number of steps in directed graph before a player is eliminated
lightningtimetoelim <- function(SS_make, LS_make, n) {
  p <- LS_make
  q <- SS_make
  #need to make function for n=2
  if (n == 2) {
    rounds <- (-1)*(p+q-3)*(p^2-p*q-2*p+2*q+1) / (q*(p-3)*(p-1)*(p-q+1))
    return(rounds)
  }
  else {
    Tr <- diag(6*n)
    for (i in 1:(5*n)) Tr[i,i] = 0
    for (i in 1:(n-2)) Tr[i,i+2] = p^2
    Tr[n-1,1] = p^2
    Tr[n,2] = p^2
    for (i in 1:(n-2)) Tr[i+n, i+2] = q^2
    Tr[2*n-1,1] = q^2
    Tr[2*n, 2] = q^2
    for (i in 1:(n-2)) Tr[i+3*n,i+2] = q*p
    Tr[4*n-1,1] = q*p
    Tr[4*n,2] = q*p
    for (i in 1:n) Tr[i,i+n] = (1.0-p)^2
    for (i in 1:n) Tr[i+n,i+n] = (1.0-q)^2
    for (i in 1:n) Tr[i+3*n, i+n] = (1.0-p)*(1.0-q)
    for (i in 1:n) Tr[i+2*n, i+2*n] = (1.0-q)^2
    for (i in 1:n) Tr[i+4*n, i+2*n]= (1.0-p)*(1.0-q)
    for (i in 1:(n-1)) Tr[i+2*n, i+3*n+1] = (1.0-q)*q
    Tr[3*n,3*n+1] = (1.0-q)*q
    for (i in 1:(n-1)) Tr[i+4*n,i+3*n+1] = (1.0-p)*q
    Tr[5*n,3*n+1] = (1.0-p)*q
    for (i in 1:(n-1)) Tr[i, i+4*n+1] = p*(1.0-p)
    Tr[n,4*n+1] = p*(1.0-p)
    for (i in 1:(n-1)) Tr[i+n, i+4*n+1] = q*(1.0-q)
    Tr[2*n,4*n+1] = q*(1.0-q)
    for (i in 1:(n-1)) Tr[i+3*n, i+4*n+1] = q*(1.0-p)
    Tr[4*n,4*n+1] = q*(1.0 - p)
    for (i in (1:n)) Tr[i,i+5*n] = (1.0-p)*p
    for (i in 1:n) Tr[i+n, i+5*n] = (1.0-q)*q
    for (i in 1:n) Tr[i+2*n,i+5*n] = q
    for (i in 1:n) Tr[i+3*n, i+5*n] = (1.0-q)*p
    for (i in 1:n) Tr[i+4*n, i+5*n] = p
  
    Q <- Tr[1:(5*n), 1:(5*n)]
    I <- diag(5*n)
    t <- I-Q
    ones <- matrix(1,5*n,1)
    x <- solve(I-Q,ones)
    return(x[1,])
  }
}
# number of rounds for longest game
lightningtimetoelim(0.9,0.4,2)*700

## Graphs/data visualization
data10players <- data.frame(
  Position = c(1:10),
  Probability = c(lightning(0.9,0.5,10))
)
data10players
lightning(0.9,0.5,10)
ggplot(data10players,aes(Position,Probability)) +
  geom_point(size=3) + 
  geom_line()

data20players <- data.frame(
  Position = c(1:20),
  Probability = c(lightning(0.9,0.5,20))
)
lightning(0.9,0.5,6)
ggplot(data20players,aes(Position,Probability)) +
  geom_point(size=3) + 
  geom_line()

lightninggraph <- function(SS_Make, LS_Make, n) {
  playerdata <- data.frame(
    Position = c(1:n),
    Probability = c(lightning(SS_Make, LS_Make, n))
  )
  graph <- ggplot(playerdata,aes(Position,Probability)) +
    geom_point(size=2) + 
    geom_line() +
    scale_x_continuous(breaks = seq(0,n))
  print(graph)
  return(playerdata)
}

lightninggraph(0.3,0.2,10)

data10players9020 <- data.frame(
  Position = c(1:10),
  Probability = c(lightning(0.9,0.2,10))
)
data10players9050 <- data.frame(
  Position = c(1:10),
  Probability = c(lightning(0.9,0.5,10))
)
data10players9080 <- data.frame(
  Position = c(1:10),
  Probability = c(lightning(0.9,0.8,10))
)
ggplot(data10players9020,aes(Position,Probability)) +
  geom_point(size=2, shape = 16) + 
  geom_line() +
  geom_point(data=data10players9050,size=2,color='black',shape = 15) +
  geom_line(data=data10players9050,color='black') +
  geom_point(data=data10players9080,size=3,color='black',shape=17) +
  geom_line(data=data10players9080,color='black')+
  scale_x_continuous(breaks = seq(0,10))+
  theme_bw()+
  labs(title="10 Player Game",
       shape="Long Shot")
