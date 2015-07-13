function [c, finalState] = ccEncode(fwdStruct, seq, initialState)
% CCENCODE  Encode the sequence seq with the trellis fwd
%
% SYNOPSIS  [c, finalState] = ccEncode(fwdStruct, seq, initialState)
%
% INPUTS    fwdStructs = (struct) forward Trellis structure
%           seq = (scalar) input sequence
%           initialState  = (scalar) initial state
%
% OUTPUTS   c  = (row-vector) codewords
%           finalState = (scalar) final state

assert((exist('private/ccEncode.mexa64', 'file') > 0) || ...    
        (exist('private/ccEncode.mexw64', 'file') > 0) || ... 
        (exist('private/ccEncode.mexw32', 'file') > 0) || ... 
        exist('private/ccEncode.mexa32', 'file') > 0, ...
      'Put or compile the executable file in the private folder.');

[c, finalState] = ccEncode(fwdStruct, seq, initialState);
