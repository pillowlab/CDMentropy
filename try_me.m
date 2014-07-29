% This is a toy script to guide you to using CDMentropy
% with a toy binary word observations.

%% Add the library path
startup

%% Compile the MEX files if not done already
if exist('discreteTimeSeries2Words') ~= 3
    disp('MEX file not found. Attempting to compile it');
    try
	makeMex
    catch
	disp('Failed to compile MEX');
    end
end

%% Toy data to play with; binary vectors of length 5
words = [...
    [0 1 0 0 0]; ...
    [1 0 0 1 0]; ...
    [0 0 0 0 0]; ...
    [0 0 1 0 0]; ...
    [1 0 0 0 1]; ...
    [0 0 0 1 0]; ...
    [1 0 0 0 1]; ...
    [1 0 0 0 1]; ...
    [1 0 0 0 1]; ...
    [1 1 1 1 0]; ...
    [1 0 0 1 0]; ...
];
m = size(words, 2); % get the dimension

%% Basic usage for estimating entropy
% opts.verbose = false; % supress verbose output
H1 = computeH_CDM(words);
fprintf('CDM entropy estimate [%f bits] under synchrony dist prior\n', H1);

%% Using the Bernoulli prior
opts.isDBer = true;
H2 = computeH_CDM(words, [], [], opts);
fprintf('CDM entropy estimate [%f bits] under Bernoulli prior\n', H2);

%% Using the compressed representation
[nn ocnts] = words2nnOcnts(words); % compact representation
H3 = computeH_CDM(nn, ocnts, m, opts);
fprintf('CDM entropy estimate [%f bits] under synchrony dist prior\n', H3);

%% Visualizing the word count histogram
% Convert the words to unique identities
symSeq = discreteTimeSeries2Words(uint16(words'), 2);
[xrange, uidx] = unique(symSeq); % Unique words
counts = histc(symSeq, xrange);
wordstrs = reshape(sprintf('%d', words(uidx,:))', [], m);

figure(8490);
bh = bar(xrange, counts);
set(bh, 'FaceColor', 'none');
set(gca, 'XTickLabel', wordstrs, 'YTick', 0:max(counts)+1);
title('Binary word histogram');
ylim([0 max(counts)+1]); ylabel('# of observations');

% If you need the variance 
opts.nMC = 999; % control the # of samples from posterior
opts.verbose = true;
[H4, Hvar] = computeH_CDM(nn, ocnts, m, opts);
fprintf('CDM entropy estimate [%f +/- %f bits]\n', H4, sqrt(Hvar));

% Generating the multiplicities form (perhaps for using other estimators)
counts = fastWords2Counts(uint16(words'), 2); % equivalent to above counts
[mm, icts] = multiplicitiesFromCounts(counts);
% [mm, icts] = words2multiplicities(words, 2); % equivalently
H5 = computeH_PYM(mm, icts); % estimate entropy using PYM for fun

% If you need the samples from the posterior
opts.verbose = false;
[H6, Hvar, CIhandle, ~, Hsamples, optsOut] = computeH_CDM(nn, ocnts, m, opts);
figure(1470); cla; hold on;
ksdensity(Hsamples);
ph = zeros(3, 1);
ph(1) = plot(H1 * [1 1], [0 0.1], 'o-');
ph(2) = plot(H2 * [1 1], [0 0.2], 'x-');
ph(3) = plot(H5 * [1 1], [0 0.3], '*-');
xlabel('entropy (bits)'); ylabel('distribution');
legend(ph, 'CDM (Syn)', 'CDM (Ber)', 'PYM');

if optsOut.isDBer
    priorStr = 'Bernoulli';
else
    priorStr = 'Synchrony';
end
title(sprintf('Posterior entropy distribution [prior: %s] nMC = %d, precision = %g', priorStr, optsOut.nMC, optsOut.hPrecision));
