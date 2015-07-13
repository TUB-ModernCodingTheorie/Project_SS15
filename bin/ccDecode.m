function [m, c] = ccDecode(bwdStruct, encodedFrame, metric, initState, finalState)
% CCDECODE  Decode the frame encodedFrame with the trellis bwd
%
% SYNOPSIS  [m, c] = ccDecode(bwdTrellis, encodedFrame, metric, initState, finalState)
%
%
% INPUTS    fwdStructs = (struct) backward Trellis structure
%           encodedFrame = (row-vector) frame to decode
%           metric  = (matrix) codewords metric
%           initState = (scalar) initial state for decoding
%           finalState = (scalar) final state for decoding
%
% OUTPUTS   m = (row-vector) decoded sequence
%           c  = (row-vector) codewords

assert((exist('private/ccDecode.mexa64', 'file') > 0) || ...    
        (exist('private/ccDecode.mexw64', 'file') > 0) || ... 
        (exist('private/ccDecode.mexw32', 'file') > 0) || ... 
        exist('private/ccDecode.mexa32', 'file') > 0, ...
      'Put or compile the executable file in the private folder.');

[m, c] = ccDecode(bwdStruct, encodedFrame, metric, initState, finalState);
