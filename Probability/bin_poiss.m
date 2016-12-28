%% Binomial & Poisson (Approximation) Distributions
% Given n independent trials and a probability of success p,
% compute the probability that we observe k = 0,1,...,n successes
% in n independent trials.

% Random variable X_i ~ Bernoulli(p)
% => Random variable X = (X_1 + ... + X_n) ~ Bin(n,p)
% For large n and small p, X is approximately Poisson(n*p)-distributed

% Input:                    n - The integer number of independent trials
%                           p - P(X_i = 1); i.e., the probability of success for any given trial
%                           k - The integer number of desired successes
% Output:                   mean - Expected value (mean) of random variable X
%                           var - Variance of random variable variable X
%                           bin_pk - P(X = k) if X ~ Bin(n,p)
%                           poiss_pk - P(X = k) if X ~ Poiss(n*p)
%                           bin_ck - P(X <= k) if X ~ Bin(n,p)
%                           poiss_ck - P(X <= k) if X ~ Poiss(n*p)
%                           error_pk - | bin_pk - poiss_pk |
%                           error_ck - | bin_ck - poiss_ck |

function [ mean, var, bin_pk, poiss_pk, bin_ck, poiss_ck, error_pk, error_ck ] = bin_poiss( n, p, k ) 

% Compute expectation E[X] = mu = lambda (Poisson parameter)
mean = n * p;

% Compute variance sigma^2 = n*p*(1 - p)
var = mean * (1 - p);

%% PDF

% Compute probability of exactly k successes, assuming X ~ Bin(n,p)
bin_pk = factorial(n) / (factorial(k) * factorial(n - k)) * (p^k) * (1 - p)^(n - k);

% Compute probability of exactly k successes, assuming X ~ Poisson(np)
poiss_pk = (mean^k) * exp(-mean) / factorial(k);

%% CDF

% Compute P(X <= k); i.e., compute the probability of k or fewer successes
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
error_pk = abs( bin_pk - poiss_pk );
error_ck = abs( bin_ck - poiss_ck );

end