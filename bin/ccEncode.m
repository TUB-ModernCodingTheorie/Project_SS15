function [c, sN] = ccEncode(fwd, s0, seq)
% HELLOWORLD  Encode the sequence seq with the treliis fwd
%
% SYNOPSIS  [c, sN] = ccEncode(fwd, s0, seq)
%
% INPUTS    fwd = forward Trellis matrix
%           s0  = initial state
%           seq = input sequence
%
% OUTPUTS   c  = codewords
%           sN = final state
%
%
% More detailed help is in the <a href="matlab: help helloWorld>extended_help">extended help</a>.


assert(isa(fwd,'double'), '1st argument should be of type double');
assert(isa(s0,'double'), '2nd argument should be of type double');
assert(isa(seq,'double'), '3rd argument should be of type double');
assert((exist('private/ccEncode.mexa64', 'file') > 0) || ...    
        (exist('private/ccEncode.mexw64', 'file') > 0) || ... 
        (exist('private/ccEncode.mexw32', 'file') > 0) || ... 
        exist('private/ccEncode.mexa32', 'file') > 0, ...
      'Put or compile the executable file in the private folder. You can use: ''cd private'' then ''mex ccEncode.c''');

[c, sN] = ccEncode(fwd, s0, seq);

function extended_help
%EXTENDED_HELP Some additional technical details and examples
%
% Here is where you would put additional examples, technical discussions,
% documentation on obscure features and options, and so on.

error('This is a placeholder function just for helptext');