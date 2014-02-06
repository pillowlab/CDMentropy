function [prior wts uwts] = computeBinomialPrior(alphas, P, N)
% prior = computeBinomialPrior(alphas, P, N)
% Accepts a vector "alphas" and returns Binomial prior for fixed P. 
%
% Inputs:
%
% Ouptuts:
%   prior = prior (d/d\alpha) E[H|\alpha, p]
%   wts = probability of falling into each partition
%   uwts = wts without the normalizer which is the size of each partition
%     (such that wts ./ uwts results in the size of the partition)
% 
% $Id$

spkRng = (0:N)';
if(isscalar(P))
    uwts = (P.^(spkRng) .* (1-P).^(N-spkRng));
    wts = binopdf(spkRng, N, P);
    wts2 = wts.*uwts;
else
    % we assume P is the count histogram - so, to get uwts, we have to
    % divide it by the binomial coefficients.
    assert(length(P) == (N+1))
    P = P/sum(P);
    wts = P(:)';
    % compute base measure
    uwts = log(wts) - gammaln(N+1)  + gammaln((0:N)+1) + gammaln(N - (0:N) +1);
    wts2  = exp(log(wts) + uwts); 
    uwts = exp(uwts);
end
uwts = uwts(:); wts = wts(:); wts2 = wts2(:);
Z = trigamma(alphas+1);
prior = zeros(size(alphas));
% Hdir = zeros(size(alphas));
for ii = 1:length(alphas)
%     Hdir(ii) = digamma(alphas(ii)+1) - (wts./uwts)'*digamma(alphas(ii)*uwts+1);
    prior(ii) = Z(ii) - wts2'*trigamma( alphas(ii)*uwts + 1);
end
