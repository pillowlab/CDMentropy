function [Hdir,Hvar] = computeHdirBaseEquiv(nn,base,alphas,unWtBase)
% Compute posterior mean of P(H|n,alpha), the expected entropy under a
% fixed Dirichlet prior with Dirichlet parameter alpha (see eq. 14 of 
%  Archer,Park,Pillow 2013).
%
% Accepts a vector "alphas" and returns posterior mean at each alpha
%
% Inputs: (Note: K is the number of equivalence classes.)
% 
%       nn [Kx1] = count data, ie, as in multinomial samples.
%     base [Kx1] = # total bins in distribution
%             * nn and base must have the same length and order.
%         alphas = scalar (or vector) of Dirichlet parameters
% unWtBase [Kx1] = 'base' represents the probability distribution
%                over equivalence classes. 'unWtBase' is the probability
%                of each element in the class. When this parameter is 
%                nonempty, each element in base is interpreted as:
%                     base (i) = (# elements in class i) * unWtBase(i)
%          Example: 
%               Take base = P, where P is some probability distribution.
%               Then, if unWtBase = P/2 this indicates that each partition
%               is an equivalence class with 2 bins.
%
% Ouptuts:
%          Hdir = mean entropy at each alpha
%          Hvar = variance of entropy at each alpha
%
% Copyright Pillow Lab 2011-2014. All rights reserved.

base = base(:);
unWtBase = unWtBase(:);
alphas = alphas(:);
K = length(base);

if isempty(nn)
    nn = zeros(K,1);
end
nn = nn(:);

assert(length(nn) == K)

N = sum(nn);
A = N+alphas.*(sum(base)); % number of effective samples (vector)

Hdir = zeros(size(alphas));
trm2 = zeros(size(alphas));

% Compute posterior mean over entropy
H0 = digamma(A+1);

invA = 1./A;
for idx = 1:length(Hdir)
    aa = nn + alphas(idx)*base;
    trm2(idx) = invA(idx)*((aa'*digamma(nn + alphas(idx)*unWtBase + 1)));
    Hdir(idx) = H0(idx) - trm2(idx);
end

end
