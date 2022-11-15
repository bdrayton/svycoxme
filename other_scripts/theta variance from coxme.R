



# theta - the random effect variance is estimated by coxme, but the variance
# for this estimate is not estimated.

# I have an analytical expression for the variance of theta. 
# it depends on the derivative and inverse of the theta penalty
# matrix. 
#
# Also it depends on the second partial derviative wrt b
#
# Also the random effects.

# After fitting a coxme model, I think these things could be 
# computed.

# Assume a shared frailty model,
# D is a diagonal matrix theta * I.
# solve(D) is 1/theta * I
# dD/dtheta is I
# Kppl is part of the returned hessian
# b is returned.























