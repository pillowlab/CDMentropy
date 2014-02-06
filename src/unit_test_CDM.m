%%%%%%%%%%%%%%%%%%%%%%%%%%
%% computeHdirBaseEquiv %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% First test computeHdirBaseEquiv which is the core of BDentropy code

%% SIMPLEST CASE POSSIBLE!
P = .1;
alpha = 3;
pvec = [P^2 P*(1-P) P*(1-P)  (1-P)^2];
counts = [4 45 1 0];

aa = alpha*pvec + counts;
A = sum(alpha*pvec + counts);

Hhat = digamma(A+1) - sum(aa.*digamma(aa+1))/A;

assert(abs(Hhat - computeHdirBaseEquiv(counts, pvec, alpha, pvec)) < 1e-10)
fprintf('.');

%% Check that equiv class for works
P = .1;
alpha = 3;
G    = [P^2 2*P*(1-P) (1-P)^2];
pvec = [P^2   P*(1-P) (1-P)^2];

a2 = alpha*pvec;
aw = alpha*G;
A = sum(alpha*G);

Hhat = digamma(A+1) - sum(aw.*digamma(a2+1))/A;

assert( Hhat == computeHdirBaseEquiv([0 0 0], G, alpha, pvec))
fprintf('.');

%% Check that equiv class works ... WITH ONE DATAPOINTseS!!!!
P = .1;
alpha = 3;
G    = [P^2 2*P*(1-P) (1-P)^2];
pvec = [P^2   P*(1-P) (1-P)^2];
counts = [1 0 0 0];

A =  alpha + sum(counts);

Hhat = digamma(A+1) - (1+alpha*pvec(2))*digamma(1 + alpha*pvec(2)+1)/A -  alpha*sum(pvec.*digamma(alpha*pvec+1))/A;

assert(Hhat == computeHdirBaseEquiv(counts, [pvec(2) pvec], alpha, [pvec(2) pvec]))
fprintf('.');

%%%%%%%%%%%%%%%%%%%%
%% computeH_CDM %%
%%%%%%%%%%%%%%%%%%%%
% Fun high dimensional convergence testy thing. 

opts = struct('verbose', false);
opts1 = opts; opts1.isDBer = true; % DBer
opts2 = opts; opts2.isDBer = false; % DSyn

%% Call with a matrix or histogram
words = zeros(20, 4);
words(3,1) = 1;
H1 = computeH_CDM(words, [], [], opts1);
[nn ocnts] = words2nnOcnts(words);
H2 = computeH_CDM(nn, ocnts, size(words, 2), opts1);
assert(H1 == H2);
fprintf('.');

H1 = computeH_CDM(words, [], [], opts2);
[nn ocnts] = words2nnOcnts(words);
H2 = computeH_CDM(nn, ocnts, size(words, 2), opts2);
assert(H1 == H2);
fprintf('.');

% The test vectors are created by fixing the seed and generating 
kTestCase = 0; testCase = [];
kTestCase = kTestCase + 1;
testCase(kTestCase).TEST_DATA = 'SAME_P';
testCase(kTestCase).Nt = 1e3;
testCase(kTestCase).P = .01;
testCase(kTestCase).ncells = 10;

kTestCase = kTestCase + 1;
testCase(kTestCase).TEST_DATA = 'RANDOM_VECTOR_P';
testCase(kTestCase).Nt = 1e3;
testCase(kTestCase).P = .01;
testCase(kTestCase).ncells = 10;

kTestCase = kTestCase + 1;
testCase(kTestCase).TEST_DATA = 'RANDOM_VECTOR_P';
testCase(kTestCase).Nt = 1e3;
testCase(kTestCase).P = .01;
testCase(kTestCase).ncells = 100;

kTestCase = kTestCase + 1;
testCase(kTestCase).TEST_DATA = 'RANDOM_VECTOR_P';
testCase(kTestCase).Nt = 1e4;
testCase(kTestCase).P = .01;
testCase(kTestCase).ncells = 200;

kTestCase = kTestCase + 1;
testCase(kTestCase).TEST_DATA = 'RANDOM_VECTOR_P';
testCase(kTestCase).Nt = 5e4;
testCase(kTestCase).P = .02;
testCase(kTestCase).ncells = 200;

nTestCase = numel(testCase);
Htrue = zeros(nTestCase, 1);
for kTestCase = 1:nTestCase
    % sampling messes up, so we need to reset the random seed in order to
    % generate the same data
    rand('seed', 19247918740); 
    TEST_DATA = testCase(kTestCase).TEST_DATA;
    Nt = testCase(kTestCase).Nt;
    P = testCase(kTestCase).P;
    ncells = testCase(kTestCase).ncells;

    switch TEST_DATA
	case 'SAME_P'
	    words = binornd(1, P, Nt, ncells);
	    Htrue(kTestCase) = -(P*log(P) + (1-P)*log((1-P))) * ncells;
	case 'RANDOM_VECTOR_P'
	    Pvec = P*rand(1,ncells);
	    words = binornd(1, repmat(Pvec, Nt,1));
	    Htrue(kTestCase) = -sum(Pvec.*log(Pvec) + (1-Pvec).*log((1-Pvec)));
    end

    %[mm icts] = multiplicitiesFromCounts(fastWords2Counts(uint16(words'), 2));
    [nn ocnts] = words2nnOcnts(words);

    Hbdp(kTestCase, 1) = computeH_CDM(nn, ocnts, ncells, opts1);
    Hbdp(kTestCase, 2) = computeH_CDM(nn, ocnts, ncells, opts2);

    %% Test sampling
    [HbdpS Vbdp, CIhandle, prior, Hsamples] = computeH_CDM(nn, ocnts, ncells, opts1);
    if abs(HbdpS - mean(Hsamples)) > sqrt(Vbdp)
	fprintf('Sampling error is too large? Numerical integral [%g], sample mean [%g], SD-sampled [%g]\n', HbdpS, mean(Hsamples), sqrt(Vbdp));
    else
	fprintf('.');
    end
end

% These results are obtained with 5000 grid points
Hbdp_expected = [...
    0.572835995741053   0.579728139078523; ...
    0.344433267887705   0.347457575135835; ...
    2.842692769576906   2.794169816933513; ...
    5.672440776777673   5.683504112074000; ...
    9.810961547855117   9.804525157169088; ...
] / log(2);

assert(max(abs(Hbdp(:) - Hbdp_expected(:))) < 1e-4, 'Test vectors');
fprintf('.\n');

fprintf('[[[ Unit tests passed! ]]]\n')
