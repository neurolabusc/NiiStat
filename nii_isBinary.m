function [isBinary, dat] = nii_isBinary (dat)
%returns true if array x has only two values
%Inputs
% x : array to be examined
%Outputs:
% isBinary : true of false
% x : array - if binary this will be 0 or 1, see examples
%Examples
% nii_isBinary([1 1 1 1]); %false - only one value
% nii_isBinary([1 0 1 1]); %true - two values
% nii_isBinary([10 0 0 0]); %true - two values
% nii_isBinary([1 0 2 1]); %false - three values
% nii_isBinary([nan nan 1 1]); %false - nans present
% [i, x] = nii_isBinary([inf -inf inf -inf]); %true x=[1 0 1 0]
% [i, x] = nii_isBinary([-1 -1 0 0]) %true, x=[0 0 1 1]

isBinary = false; 
mn = min(dat(:));
mx = max(dat(:));
if mn == mx
    return; %not binary - only one value
end
nMin = sum(dat(:)==mn);
nMax = sum(dat(:)==mx);
if (nMin+nMax) ~= numel(dat(:)) %more than two values
   return; 
end
isBinary = true; 
if nargout < 2, return; end;
if (mx == 0) %special case, requires more memory...  [i, x] = nii_isBinary([-1 -1 0 0])
    isMin = (dat == mn);
    dat(isMin) = 0;
    dat(~isMin) = 1;
    return;
end
dat(dat == mn) = 0;
dat(dat == mx) = 1;
    
