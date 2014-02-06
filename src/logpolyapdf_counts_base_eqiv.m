function logp = logpolyapdf_counts_base_eqiv(nn, base, alphas, unWtBase)
% logp = logpolyapdf_counts(nn,K,alphas)
%
% Compute the log-polya probability of data nn (a histogram) given Dirichlet
% prior with parameter alpha.

% Make alphas a column vector, zcts and ncts row vectors
alphas = alphas(:);
nn     = nn(:);
base   = base(:);

N = sum(nn);
A = alphas.*sum(base); % number of effective samples (vector)

% trm1 = gammaln(N+1) - sum(gammaln(nn+1)) ...
%     + gammaln(A) - gammaln(N+A);
trm1 = gammaln(N+1) - sum(gammaln(nn+1)) - gammalndiff(A,N);

trm2 = zeros(size(alphas));
nzIdx = nn>0;

for idx = 1:length(alphas)    
    trm2(idx) = sum(gammalndiff(alphas(idx)*unWtBase(nzIdx), nn(nzIdx)));
end

logp = trm2 + trm1;
