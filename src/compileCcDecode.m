% Compile all C-file to mex binaries

inputPath = '../src/c/';
outputPath = '../bin/private/';

% compilation of all c-files
srcFiles = {'ccDecode.c'};

for i = 1:length(srcFiles)
    display(['Compiling ', srcFiles{i}]);
    mex([inputPath srcFiles{i}],'-outdir',outputPath);
end