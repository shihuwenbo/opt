import math
import numpy as np

# create simulated data
mu = 1.0
sigma2 = 100.0
m = 10000
data = np.random.normal(mu, math.sqrt(sigma2), m)

# compute sample mean, sample variance
sample_mean = np.mean(data)
sample_var = np.var(data)
print sample_mean, sample_var

# iterate through different rho
all_rho = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
for rho in all_rho:
    
    # estimate mu and sigma2
    mu_est = 0.0
    sigma2_est = 0.0

    for i in range(0,10000):
    	sigma2_est = np.sum(np.square(data-mu_est))/m
	tmp = rho*sigma2_est/m
	if(sample_mean >= -tmp and sample_mean <= tmp):
	    mu_est = 0.0
	elif(sample_mean > tmp):
	    mu_est = sample_mean - tmp
	else:
	    mu_est = sample_mean + tmp
    
    # print out result
    print rho, mu_est, sigma2_est
