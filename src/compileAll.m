inputPath = '..\src\c\';
outputPath = '..\bin\';
tmpPath = [outputPath,'tmp\'];

if ~exist(tmpPath,'dir')
    mkdir(tmpPath);
end

delete([tmpPath,'*']);

srcFiles = dir(inputPath);

%-------------------- compiler --------------------------

for i = 3:length(srcFiles)
    copyfile([inputPath,srcFiles(i).name],[tmpPath,srcFiles(i).name]);
    
    mex([tmpPath,srcFiles(i).name],'-outdir',outputPath);
end

%--------------------------------------------------------