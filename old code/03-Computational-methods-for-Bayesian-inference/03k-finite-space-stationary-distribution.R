### Example of Finite-space Discrete Markov chains
### Taken from Section 6.2 of Bayesian Computing with R by 
### Jim Albert

# Transition matrix for six-state discrete Markov process
P = matrix(c(.5,.5,0,0,0,0,.25,.5,.25,0,0,0,0,.25,.5,.25,0,0,
           0,0,.25,.5,.25,0,0,0,0,.25,.5,.25,0,0,0,0,.5,.5),
           nrow=6,ncol=6,byrow=TRUE)

# number of steps
B = 50000
theta = numeric(B) # store each sampled value

# initial state
theta[1] = 3

for (j in 2:B) {
   theta[j] = sample(1:6, size = 1, prob = P[theta[j - 1],])
}

# look at marginal distribution for different number of iterations
m = c(500,2000,8000,50000)
for (i in 1:4) print(table(theta[1:m[i]])/m[i])

# stationary distribution
w = matrix(c(.1,.2,.2,.2,.2,.1), nrow = 1, ncol = 6)
w %*% P

