function [forwardStruct, backwardStruct] = ccInitialize(A,B,C,D)
% CCINITIALIZE  
%
% SYNOPSIS  treillis = ccInitialize(A,B,C,D)
%
% INPUTS    A
%           B
%           C
%           D
%
% OUTPUTS   forwardStruct = struct:
%                               forward: array      
%                               ldStates: int scalar = m
%                               ldOutputs: int scalar = n
%                               ldInputs: int scalar = k
%           forwardStruct = struct:
%                               backward: array      
%                               ldStates: int scalar = m
%                               ldOutputs: int scalar = n
%                               ldInputs: int scalar = k
%
%
% More detailed help is in the <a href="matlab: help helloWorld>extended_help">extended help</a>.


% assert(isa(fwd,'double'), '1st argument should be of type double');
% assert(isa(s0,'double'), '2nd argument should be of type double');
% assert(isa(seq,'double'), '3rd argument should be of type double');
assert((exist('private/ccInitialize.mexa64', 'file') > 0) || ...    
        (exist('private/ccInitialize.mexw64', 'file') > 0) || ... 
        (exist('private/ccInitialize.mexw32', 'file') > 0) || ... 
        exist('private/ccInitialize.mexa32', 'file') > 0, ...
      'Put or compile the executable file in the private folder.');

[forwardStruct, backwardStruct] = ccInitialize(A,B,C,D);

function extended_help
%EXTENDED_HELP Some additional technical details and examples
%
% Here is where you would put additional examples, technical discussions,
% documentation on obscure features and options, and so on.

error('This is a placeholder function just for helptext');