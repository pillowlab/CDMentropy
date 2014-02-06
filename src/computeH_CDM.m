function  [Hbls, Hvar, CIhandle, internal, Hsamples, opts] = computeH_CDM(nn, ocnts, ncells, opts)
% [Hbls, Hvar, CIhandle, internal, Hsamples, opts] = computeH_CDM(words, [], [], opts)
% [Hbls, Hvar, CIhandle, internal, Hsamples, opts] = computeH_CDM(nn, ocnts, ncells, opts)
% 
% Centered Dirichlet Mixture (CDM) entropy estimator as discussed in:
%
%     Evan Archer, Il Memming Park, and Jonathan Pillow
%     Bayesian entropy estimation for binary spike train data using parametric 
%     prior knowledge
%     Advances in Neural Information Processing Systems (NIPS), 2013.
% 
% This estimator is only for binary observations. Two different base measures
% can be applied: (1) binomial (Dber), (2) synchrony distribution (Dsyn).
%
% If only the Hbls is requested, no sampling is done (faster).
%
% Sampling can be done with MATLAB parallel toolbox. (it uses parfor)
%
% INPUT:
%   words - [time x ncells] matrix of zeros and ones
%      nn - histogram (if nn is empty, computes prior expected entropy)
%   ocnts - number of spikes for each bin in the histogram
%  ncells - length of binary observations 
%	    (number of cells for simultaneously-recorded neural data)
%    opts - Option structure (optional)
%	    opts.verbose: true for more verbose processing
%	    opts.isDBer - (optional/false)
%               *  true: use binomial prior (Dber)
%               * false: use empirical synchrony distrbution  (Dsyn)
%                      (note: synchrony distribution is regularized
%                        by adding a pseudocount of 1/2 to each bin)
%
%	    == numerical integration options ==
%	    opts.nAlpha: (default/500) # of grid points for integration
%	    opts.alphas: [opts.nAlpha x 1] (default/"Chebyshev")
%		specify grid points manually. Note that Chebyshev grid is
%		dynamic; it depends on data.
%		e.g.) opts.alphas = ...
%			logspace(log10(1e-5), log10(1e3*sum(nn)), opts.nAlpha)';
%
%	    == sampling related options ==
%	    opts.nMC (default/999) number of independent samples to draw
%	    opts.hPrecision: (default/1e-2) precision of sampled entropy
%		smaller hPrecision makes the sampling slower.
%
% OUTPUT:
%     Hbls - point estimate of entropy in bits (Bayes' Least Squares)
%	    (multiply by log(2) to obtain nats)
%     Hvar - posterior variance of entropy (estimated via sampling)
% CIhandle - @(c) function handle that returns credible interval vector
%	    E.g. CIhandle(0.05) returns 95% credible interval around the
%	    median.
% internal - an internal data structure for debugging and peeking into
%	    the computational details of CDM
% Hsamples - actual samples from the posterior
%     opts - options strucuture filled with default values, just in case.
%
% Copyright 2013-2014 Pillow Lab. See LICENSE.txt license information.

if nargin < 2 || isempty(ocnts)
    assert(nargin < 2 || isempty(ncells), ...
	'if second argument is empty, so must the third');
    words = nn;
    ncells = size(nn, 2);
    [nn ocnts] = words2nnOcnts(words);
end

if nargin < 4
    opts = struct();
end

if ~isfield(opts, 'verbose');
    opts.verbose = true;
end

if ~isfield(opts, 'isDBer') || isempty(opts.isDBer)
    opts.isDBer = false;
end

ocnts = ocnts(:); % make sure it's a column vector
if ~isempty(ocnts)
    assert(max(ocnts) <= ncells, ...
	'Cannot have more spikes than the the number of bins');
end

%% Check whether there may be numerical problems due to potentially high entropy
P = dot(ocnts, nn)/sum(nn)/ncells;
Hmax = -(P*log(P) + (1-P)*log((1-P))) * ncells;

if Hmax > 23 % warning value of 23 chosen heuristically
    warning('Maximum possible entropy indicates potential numerical problems.');
end

%% Setup the grid points and weights for the integration over alpha
if ~isfield(opts, 'nAlpha')
    opts.nAlpha = 500; % number of grid points for alpha
end

if ~isfield(opts, 'alphas') || isempty(opts.alphas)
    % Use Chebyshev spectral integration method (Boyd 1986)
    maxval = cot(pi/(opts.nAlpha+1))^2;
    aw = zeros(1,opts.nAlpha);
    tt = pi*(1:opts.nAlpha)/(opts.nAlpha+1);
    jj = 1:opts.nAlpha;
    ax = cot(tt/2).^2; % the gridpoints at which we evaluate the integrand
    L = 2^min(23, ncells)/maxval;
    L = max(L, sum(nn));
    % Heuristic: max alpha must be orders of magnitude larger than the number of
    % observations, so that the base measure can override the histogram part
    % otherwise, its performance will be dominated by data 
    % (i.e. plugin estimator)
    alphas = L*ax; alphas = alphas(:);
    z0 = (1-cos(pi*jj))./jj;
    for idx = 1:opts.nAlpha
	aw(idx) =  sum(sin(idx*tt).*z0);
    end
    wts = L*( 2*sin(tt)./(1-cos(tt)).^2 ) .* aw .* ( 2/(opts.nAlpha+1) );
    wts = wts(:);
else
    % Use user specified grid for integration
    alphas = opts.alphas;
    wts = diff([0; alphas]);
end

%% Compute prior over the grid of alphas
spkRng = (0:ncells)';
if opts.isDBer
    if isempty(nn)
	warning('CDM:no_data', 'No valid data is provided. Assuming prior.');
	Phat = 0.5;
    else
	Phat = dot(ocnts, nn)/sum(nn)/ncells;
    end
    if opts.verbose
       	fprintf('DBer mode: %f\n', Phat);
    end

    if Phat == 0 || Phat == 1
	warning('CDM:constant_neuron', 'This is a crazy constant neuron');
	Hbls = 0; Hvar = 0; internal = []; Hsamples = [];
	CIhandle = @(a) [NaN, NaN];
	return;
    end
else
    % Compute empirical count histogram
    empCntHist = accumarray([ocnts;spkRng]+1, [nn;zeros(length(spkRng),1)])';
    % how many of each class of word (equal # spks) did we see?
    if opts.verbose
       	fprintf('DSyn mode\n');
    end
    Phat = empCntHist+1/length(empCntHist); Phat = Phat/sum(Phat);
    Phat(:) = Phat(:)';
end

[dHdir, H, unWgtProb] = computeBinomialPrior(alphas,Phat,ncells);
H = reshape(H,1,[]); unWgtProb = reshape(unWgtProb,1,[]);
% note that (H ./ unWgtProb) is the size of each partition

%%
binCnt = hist(ocnts, spkRng); binCnt = binCnt(:)';
% Compute a new base measure, removing the counts of those words we have
% observed. This may be questionable numerically.
dH = H - binCnt.*unWgtProb;
% We don't need to consider equivalence classes of which we've observed all
% elements.
postH = dH;
postH(abs(postH) < 2 * eps) = 0; % numerical error happens

nn_concat = [nn(:); zeros(length(postH(dH>eps)),1)];

% Here we create the base measure for the posterior which reflects our
% observations. There's some redundancy right now in that the prior
% probabilities of the observed data are represented twice. 
postBase = [unWgtProb(ocnts+1)'; postH(dH>eps)'];
postUnWgtBase = [unWgtProb(ocnts+1)'; unWgtProb(dH>eps)']; 
postHdir = computeHdirBaseEquiv(nn_concat, postBase, alphas, postUnWgtBase);

% Compute evidence p(n|a)
logpn_a = logpolyapdf_counts_base_eqiv(nn_concat, postBase, alphas, postUnWgtBase);
[maxlogp,imx] = max(logpn_a);
pn_a = exp(logpn_a - maxlogp);

% Compute posterior p(a|n)
pa_n = pn_a.*dHdir;

% compute posterior mean E[H]
Z = wts' * pa_n;
Hbls = wts' * (postHdir .* pa_n) / Z;
Hbls = Hbls / log(2);

%% Pack intermediate variables in a struct for debugging purposes
if nargout > 2
    internal.pa = dHdir;
    internal.alphas = alphas;
    internal.logpn_a = logpn_a;
    internal.pn_a = pn_a;
    internal.pa_n = pa_n;
    internal.ml_a = alphas(imx);
end

if nargout > 1
    %% Sample the posterior to get confidence interval and variance and what not
    % We divide the words into two cases: the observed, and the unobserved
    %
    % For each sample, we take the following steps
    % 1. sample an alpha from pa_n = pn_a * pa / Z
    % 2. get the base measure over the observed words and the partitions
    % 3. sample from Dirichlet
    % 4. Do stick breaking within each partition until enough probability mass
    %    is covered
    % 5. compute entropy of the sample distribution

    if ~isfield(opts, 'nMC') || isempty(opts.nMC)
	opts.nMC = 999;
    end

    if ~isfield(opts, 'hPrecision') || isempty(opts.hPrecision)
	opts.hPrecision = 1e-2;
    end

    Hsamples = nan(opts.nMC, 1);

    parfor kMC = 1:opts.nMC % make this a parfor if you wish to parallelize
	if opts.verbose; fprintf('[%d/%d]\r', kMC, opts.nMC); end
	% Step 1
	a = randsample(alphas, 1, true, pa_n.*wts);

	% Step 2
        if isempty(nn)
	    aGeq = a * postH';
	else
	    aGeq = [nn + a * unWgtProb(ocnts+1)'; a * postH'];
        end
        aGeq(aGeq < 0) = 0; % numerical problems

	% Step 3
	s = warning('query', 'PYM:dirchletrnd');
	warning('off', s.identifier);
	p = dirichletrnd(aGeq, numel(aGeq), 1);
	warning(s.state, s.identifier);

	% Step 4
	pCell = cell(numel(postH) + 1, 1); % pi's within each class
	pCell{1} = p(1:numel(nn))';
	HnotSampled = 0;
	for kPartition = 1:numel(postH)
	    % break a stick from Dir(a, a, ..., a)
	    % with binCnt(kPartition) number of a's

	    nA = round(dH(kPartition) / unWgtProb(kPartition));
	    aPartition = aGeq(numel(nn) + kPartition) / nA;
	    % if alpha is large, use the expected entropy instead of sampling
	    isAdded = false;
	    pThreshold = opts.hPrecision / ...
		abs(log(dH(kPartition)) - log(unWgtProb(kPartition)));

	    if p(numel(nn)+kPartition) > 10*eps ...
		    && pThreshold < p(numel(nn)+kPartition)
		[Hdir, Hvar] = computeHdir([], [], nA, aPartition);
		if p(numel(nn)+kPartition) * sqrt(Hvar) < opts.hPrecision
		    HnotSampled = HnotSampled + ...
			(Hdir - log(p(numel(nn)+kPartition))) ...
			    * p(numel(nn)+kPartition);
		    if opts.verbose; fprintf('.'); end
		    isAdded = true;
		end

		if ~isAdded && dH(kPartition) > 10 * eps
		    % number of alphas within this partition:
		    if nA <= 1
			continue;
		    end
		    if opts.verbose; fprintf('o'); end

		    %%%%%% DEBUG - no size biased sampling =======
		    % pCell{1 + kPartition} = ...
		    %     p(numel(nn) + kPartition) * dirichletrnd(aPartition, nA, 1)';
		    % continue;
		    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		    Apartition = aPartition * nA;
		    kaPartition = -aPartition;

		    % draw the first stick
		    x = betarnd(1 + aPartition, Apartition + kaPartition);
		    pPartition = zeros(1, min(100, nA));
		    pPartition(1) = x;
		    pRemaining = 1 - x;
		    for kSticks = 2:nA
			kaPartition = kaPartition - aPartition;
			x = betarnd(1+aPartition, max(Apartition + kaPartition, 0));

			pPartition(kSticks) = x * pRemaining;
			pRemaining = (1 - x) * pRemaining;

			% if enough probability mass is sampled, stop.
			if kSticks ~= nA && pRemaining < pThreshold
			    pPartition(kSticks + 1) = pRemaining;
			    kSticks = kSticks + 1;
			    break;
			end
		    end
		    if isempty(kSticks); kSticks = 1; end % this is when nA < 2

		    pPartition = pPartition(1:kSticks);
		    pCell{1+kPartition} = p(numel(nn)+kPartition) * pPartition;
		end
	    else
		pCell{1 + kPartition} = p(numel(nn)+kPartition);
		if opts.verbose; fprintf('+'); end % too small to care!
	    end
	end
	pp = [pCell{:}]; pp = pp(:); % collect sticks from each partition

	% Step 5
	vidx = (pp ~= 0) & ~isnan(pp);
        Hsamples(kMC) = -sum(pp(vidx).*log(pp(vidx))) + HnotSampled;
    end % kMC

    Hsamples = Hsamples / log(2); % convert from nats to bits
    CIhandle = @(a) quantile(Hsamples, [a/2, 1 - a/2]);
    HmeanSampled = mean(Hsamples);
    if opts.verbose; 
	ci = CIhandle(0.05);
	fprintf('BLS [%g] sample mean [%g], 95%% CI [%g][%g]\n', ...
		Hbls, HmeanSampled, ci(1), ci(2));
    end

    Hvar = var(Hsamples);
end

end
