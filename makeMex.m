% Compiles the mex files for fast histogram generation

cd lib/PYMentropy/src/
try
    fprintf('Mexing fastWords2Counts\n');
    mex fastWords2Counts.c fastWordsTree.c
    fprintf('Done!\n');
    fprintf('Mexing discreteTimeSeries2Words\n');
    mex discreteTimeSeries2Words.c fastWordsTree.c
    fprintf('Done!\n');
catch
    fprintf('Error!\n');
end
cd ../../..
