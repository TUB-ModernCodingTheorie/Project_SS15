% Compile all C-file to mex binaries

inputPath = '../src/c/';
outputPath = '../bin/private/';

% compilation of all c-files
srcFiles = dir(inputPath);

for i = 1:length(srcFiles)
    if ~isempty(regexp(srcFiles(i).name,'\w*(\.c)$','once'))
        disp(['Compiling ', srcFiles(i).name]);
        mex([inputPath srcFiles(i).name],'-outdir',outputPath);
    end
end