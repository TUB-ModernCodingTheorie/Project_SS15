function [c, sN] = ccDecode(encodedFrame, bwd, initState, finalState)
% CCDECODE  Decode the frame encodedFrame with the trellis bwd
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


% assert(isa(fwd,'double'), '1st argument should be of type double');
% assert(isa(s0,'double'), '2nd argument should be of type double');
% assert(isa(seq,'double'), '3rd argument should be of type double');
assert((exist('private/ccDecode.mexa64', 'file') > 0) || ...    
        (exist('private/ccDecode.mexw64', 'file') > 0) || ... 
        (exist('private/ccDecode.mexw32', 'file') > 0) || ... 
        exist('private/ccDecode.mexa32', 'file') > 0, ...
      'Put or compile the executable file in the private folder.');

[c, sN] = ccDecode(encodedFrame, bwd, initState, finalState);

function extended_help
%EXTENDED_HELP Some additional technical details and examples
%
% Here is where you would put additional examples, technical discussions,
% documentation on obscure features and options, and so on.

error('This is a placeholder function just for helptext');