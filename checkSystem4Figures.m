function [figOpts, ofgOpts] = checkSystem4Figures()
%CHECKSYSTEM verifies that the code is running in either a 64 bit windows
%based or Linux based machine and sets figure options to be invisible
%   When running MATLAB code in the HPC, it's not a good idea to display
%   figures or any graphic object.
figOpts = {'Visible', 'on'};
ofgOpts = {'visible'};
if ~strcmp( computer, 'PCWIN64')
    figOpts{2} = {'off'};
    ofgOpts{2} = {'invisible'};
end

end