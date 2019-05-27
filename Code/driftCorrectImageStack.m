function [imStackCorr,maxCorr,imTform] = driftCorrectImageStack(imStack,varargin)
%DRIFTCORRECTIMAGESTACK performs drift correction / image stabilization / intra-stack registration on the input image stack
%
% imStackCorr = driftCorrectImageStack(imStack)
% [imStackCorr,xyTranslation,refineXform] = driftCorrectImageStack(imStack)
%
% This function attempts to correct for drift/misalignment in the input
% image stack in a two-step approach. First the images are aligned by
% maximum cross-correlation (pixel-level translation only), and then this
% intiial alignment is (optionally) refined using imregister (sub-pixel,
% supports translation in addition to rigid, affine, etc.)
%
% Input:
%
%   imStack - 3D stack of 2D images. (Though it wouldn't be too much work
%   to generalize this to 4D stacks of 3D... feel free to do so!!)
%
%   Also, a bunch of parameters that you'll have to open the function to
%   see descriptions of because I'm too lazy to re-document here.
%
%Hunter Elliott
%7/2013
%
% Copyright (C) 2019, Jaqaman & Danuser Labs - UTSouthwestern 
%
% This file is part of CondensateAnalysis.
% 
% CondensateAnalysis is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% CondensateAnalysis is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with CondensateAnalysis.  If not, see <http://www.gnu.org/licenses/>.
% 
% 


%% -------------- Input ------------ %%

ip = inputParser;
ip.addRequired('imStack')
ip.addParamValue('ReferenceImage',-1,@(x)(isscalar(x) && isequal(x,round(x)) && x ~= 0 ));%Image to correct to. Negative numbers are relative e.g. -1 = previos image
ip.addParamValue('DoRefinement',true,@islogical);%Whether to run imregister after the correlation-based translation to improve the registration
ip.addParamValue('Verbose',false,@islogical);%Whether to display text status
ip.addParamValue('TransformType','translation');%See imregister.m
ip.addParamValue('Modality','monomodal')%See imregister.m
ip.addParamValue('FillValue',0,@isscalar);%Value to fill in pixels which are empty due to correction
ip.parse(imStack,varargin{:});

p = ip.Results;

%% ----------- Init --------- %%

imSize = size(imStack);

imStackCorr = imStack;

%Get indices of reference images and corresponding images to correct
if p.ReferenceImage > 0
    %If a constant reference image was specified
    refInd = ones(imSize(end),1) * p.ReferenceImage;
    corrInd = [1:p.ReferenceImage-1 p.ReferenceImage+1:imSize(end)];
else
    refInd = 1:imSize(end)+p.ReferenceImage;
    corrInd = -p.ReferenceImage+1:imSize(end);    
end

nCorr = numel(corrInd);

[optimizer,metric] = imregconfig(p.Modality);

%We need to keep track of filled in values in case they're NaNs
trackFill = false(imSize);

maxCorr = nan(imSize(end),2);
imTform = cell(imSize(end),1);

%% ----------- Correction --------- %%

for j = 1:nCorr
       
    currRef = double(imStackCorr(:,:,refInd(j)));    
    currCorr = double(imStack(:,:,corrInd(j)));
    
    if p.Verbose;tic;end;
        
    %Get cross-orrelation between images to do an initial translation since imregister is slow (optimization
    %based)    
    imXcorr = convnfft(currRef - mean(currRef(:)),rot90(currCorr,2)-mean(currCorr(:)));
    %Find max correlation and shift accordingly
    [~,iMax] = max(imXcorr(:));    
    [maxCorr(corrInd(j),1), maxCorr(corrInd(j),2)] = ind2sub(size(imXcorr),iMax);        
    [Xi,Yi] = meshgrid((1:imSize(2))+imSize(2)-maxCorr(corrInd(j),2),(1:imSize(1))+imSize(1)-maxCorr(corrInd(j),1));
    currCorr = interp2(currCorr,Xi,Yi,'nearest');%Lazy way to do shifting + fill NaNs
    trackFill(:,:,corrInd(j)) = isnan(currCorr);%Log these so they don't mess up the correlations/registrations and can be put back alter
    currCorr(isnan(currCorr)) = 0;
       
    %Optionally improve the registration
    if p.DoRefinement
        %We can't call imregister directly because we want to specify fill
        %values
        imTform{corrInd(j)} = imregtform(currCorr,currRef,p.TransformType,optimizer,metric);        
        Rfixed = imref2d(size(currRef));
        currCorr = imwarp(currCorr,imTform{corrInd(j)},'OutputView',Rfixed,'FillValues',NaN);
        trackFill(:,:,corrInd(j)) = trackFill(:,:,corrInd(j)) | isnan(currCorr);        
        currCorr(isnan(currCorr)) = 0;
    end
    
    imStackCorr(:,:,corrInd(j)) = currCorr;
    
    if p.Verbose;disp(['Finished registering image ' num2str(corrInd(j)) ' to image ' num2str(refInd(j)) ',took ' num2str(toc) ' seconds']);end
    
end

imStackCorr(trackFill) = p.FillValue;


