%% Binomial & Poisson Distributions
% Given n independent trials and a probability of success p,
% compute the probability that we observe k = 0,1,...,n successes
% in n independent trials.

% Random variable X_i ~ Bernoulli(p)
% => Random variable X = (X_1 + ... + X_n) ~ Bin(n,p)
% For large n and small p, X is approximately Poisson(n*p)-distributed

% Input                     n - The number of independent trials
%                           p - P(X_i = 1); i.e., the probability of 
%                           success for any given trial
%                           k - The number of desired successes
% Output                    pk - 

function [ mean, var, bin_pk, poiss_pk, bin_ck, poiss_ck, error ] = bin_poiss( n, p, k ) 

% Compute expectation E[X] = mu = lambda (Poisson parameter)
mean = n * p;

% Compute variance sigma^2 = n*p(1 - p)
var = mean * (1 - p);

% Compute probability assuming X ~ Bin(n,p)
bin_pk = factorial(n) / (factorial(k) * factorial(n - k)) * (p^k) * (1 - p)^(n - k);

% Compute probability assuming X ~ Poisson(np)
poiss_pk = (mean^k) * exp(-mean) / factorial(k);

% Compute P(X <= k)
% Initialize bin_ck, poiss_ck
bin_ck = 0; poiss_ck = 0;
for i = 0:k
   % Assume X ~ Bin(n,p)
   bin_ck = bin_ck + ...
       factorial(n) / (factorial(i) * factorial(n - i)) * (p^i) * (1 - p)^(n - i);
   
   % Assume X ~ Poisson(np)
   poiss_ck = poiss_ck + (mean^i) * exp(-mean) / factorial(i);
end

% Compute the error between bin_pk & poiss_pk
error = abs(bin_pk - poiss_pk);

end