% Compile all C-file to mex binaries

inputPath = '../src/c/';
outputPath = '../bin/private/';

% compilation of all c-files
srcFiles = dir(inputPath);

for i = 3:length(srcFiles)
    display(['Compiling ', srcFiles(i).name]);
    mex([inputPath srcFiles(i).name],'-outdir',outputPath);
end