function [forwardStruct, backwardStruct] = ccInitialize(A,B,C,D)
% CCINITIALIZE  
%
% SYNOPSIS  [forwardStruct, backwardStruct] = ccInitialize(A,B,C,D)
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

assert((exist('private/ccInitialize.mexa64', 'file') > 0) || ...    
        (exist('private/ccInitialize.mexw64', 'file') > 0) || ... 
        (exist('private/ccInitialize.mexw32', 'file') > 0) || ... 
        exist('private/ccInitialize.mexa32', 'file') > 0, ...
      'Put or compile the executable file in the private folder.');

[forwardStruct, backwardStruct] = ccInitialize(A,B,C,D);
