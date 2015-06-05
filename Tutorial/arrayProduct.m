function c = arrayProduct(multiplier, matrix)
% ARRAYPRODUCT  Multiply a scalar and a row matrix
%
% SYNOPSIS  c = arrayProduct(multiplier,matrix)
%
% INPUTS    multiplier = numeric scalar
%           matrix = numeric row vector
%
% OUTPUTS   c = numeric row vector
%
% This function is a matlab tutorial function.
% It's a good way to learn how it works
%
% More detailed help is in the <a href="matlab: help helloWorld>extended_help">extended help</a>.

% Examples:
% c = arrayProduct(4,[1,3,13]) 
%

assert(isscalar(multiplier) && isnumeric(multiplier), ...
        '1st argument should be a scalar');
assert(isrow(matrix) && isnumeric(matrix), ...
        '2nd argument should be a row vector of numbers');
assert((exist('private/arrayProduct.mexa64', 'file') > 0) || ...                 
        exist('private/arrayProduct.mexa32', 'file') > 0, ...
      'Put or compile the executable file in the private folder. You can use: ''cd private'' then ''mex arrayProduct.c''');

multiplier = double(multiplier);
matrix = double(matrix);

c = arrayProduct(multiplier, matrix);

function extended_help
%EXTENDED_HELP Some additional technical details and examples
%
% Here is where you would put additional examples, technical discussions,
% documentation on obscure features and options, and so on.

error('This is a placeholder function just for helptext');