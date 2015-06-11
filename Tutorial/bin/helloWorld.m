function helloWorld(name)
% HELLOWORLD  Print "hello" followed by the value of name
%
% SYNOPSIS  helloWorld(name)
%
% INPUTS    name = string
%
% OUTPUTS   none
%
% This was just a test function to check that everything works fine.
% It's a good way to learn how it works
%
% More detailed help is in the <a href="matlab: help helloWorld>extended_help">extended help</a>.

% Examples:
% helloWorld('Olivier')
%

assert(ischar(name), '1st argument should be of type char');
assert((exist('private/helloWorld.mexa64', 'file') > 0) || ...                 
        exist('private/helloWorld.mexa32', 'file') > 0, ...
      'Put or compile the executable file in the private folder. You can use: ''cd private'' then ''mex helloWorld.c''');

helloWorld(name);

function extended_help
%EXTENDED_HELP Some additional technical details and examples
%
% Here is where you would put additional examples, technical discussions,
% documentation on obscure features and options, and so on.

error('This is a placeholder function just for helptext');