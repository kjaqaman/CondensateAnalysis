function [lmax]  = mapIntensity(MD,channelNum) 
%MAPINTENSITY takes local maxima from detection process of one channel and applies to other channels to read intensity from them
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

%Input
%   MD: MovieData Object with detections acquired by SubResolution Process 
%
%   channelNum: Channel number to be analyzed
%
%Output:
%
%   lmax: a structure that contains the following fields for each local
%   maxima for every time frame
%        .amp : the mean intensity of a local maxima
%        .ampBG : the mean intensity of the local background
%        .ampCorr : the difference of .amp and .ampBG
%        .stdBG : the standard deviation of the local background
%   Note: lmax is automatically saved at the location [MD.outputDirectory_   '/localMaxima_' num2str(channelNum) '.mat']
%
% Tony Vega 2017




%Get pixelIndx
iM = MD.getProcessIndex('ThreshLocMaxProcess',1,0);
p = parseProcessParams(MD.processes_{iM},[]);
inDetectDir = MD.processes_{iM}.outFilePaths_(p.ChannelIndex);
load(inDetectDir{1});

%Get movieInfo
iM = MD.getProcessIndex('SubResolutionProcess',1,0);
p = parseProcessParams(MD.processes_{iM},[]);
inDetectDir = MD.processes_{iM}.outFilePaths_(p.ChannelIndex);
load(inDetectDir{1});

%First go through each frame
for m = 1:length(localMaxima)%pixelRef
    
    %Load image
    imageRaw = double(MD.channels_(channelNum).loadImage(m));
    
    %Get pixel coordinates of local maxima
    xyCoord =[round(movieInfo(m).yCoord(:,1)),round(movieInfo(m).xCoord(:,1))];
    
    %Convert to indices
    lmInd = sub2ind(MD.imSize_,xyCoord(:,1),xyCoord(:,2));
    
    %Get binary mask where 1 is cluster area and 0 is background
    segMask = pixelIndx(:,:,m);
    
    %Now get similar mask from detection
    maskLocalMax = zeros(MD.imSize_);
    maskLocalMax(lmInd(:)) = 1;
    SE = strel('disk',2);
    
    %Make composite to ensure total area is accounted for
    maskLocalMaxDil = imdilate(maskLocalMax,SE);
    segMaskCombo= maskLocalMaxDil+segMask;
    segMaskCombo= segMaskCombo>0;
    
    %Dilate again and get label identity
    segMaskComboDil=imdilate(segMaskCombo,SE);
    l = bwlabel(segMaskComboDil);
    
    %Initialize variables to collect intensity, background and background intensity for each
    %detection
    frameLM = zeros(length(movieInfo(m).xCoord),1);
    frameBG = zeros(length(movieInfo(m).xCoord),1);
    frameBGStd = zeros(length(movieInfo(m).xCoord),1);
    
    %In every frame, look at each local maxima and try to locate
    %nearest segmented object
    for k = 1:length(lmInd)
        %First get intensity from local maxima
        maskIm = zeros(MD.imSize_);
        maskIm(lmInd(k)) = 1;%pixelRef{m}{k}
        maskIm = logical(maskIm);
        SE = strel('disk',2);
        maskLm = imdilate(maskIm,SE);
        statsLM = regionprops(maskLm,imageRaw,'MeanIntensity','PixelValues');
        frameLM(k,1) = statsLM(1).MeanIntensity;%sum(statsLM(1).PixelValues);
        %Now check to see if local maxima exists in a segmented blob
        lmArea = l(lmInd(k));
        if lmArea > 0
            %If the local maxima is within a segmented blob, use a ring
            %around the segmented blob to get the background intensity,
            lmMask = (l == lmArea);
            maskBg = imdilate(lmMask,SE);
            maskBgF = logical(maskBg-lmMask);
            statsBG = regionprops(maskBgF,imageRaw,'MeanIntensity','PixelValues');
            frameBG(k,1) = statsBG(1).MeanIntensity;%sum(statsBG(1).PixelValues);
            frameBGStd(k,1) = std(statsBG(1).PixelValues);
        else
            %Otherwise, use a ring around maskLm object to get
            %background intensity
            maskBg = imdilate(maskLm,SE);
            maskBgF = logical(maskBg-maskLm);
            statsBG = regionprops(maskBgF,imageRaw,'MeanIntensity','PixelValues');
            frameBG(k,1) = statsBG(1).MeanIntensity;%sum(statsBG(1).PixelValues);
            frameBGStd(k,1) = std(statsBG(1).PixelValues);
        end
        
    end
    lmax(m,1).amp = frameLM;
    lmax(m,1).ampBG = frameBG;
    lmax(m,1).ampCorr = frameLM-frameBG;
    lmax(m,1).stdBG= frameBGStd;
end

save([MD.outputDirectory_ '/localMaxima_' num2str(channelNum) '.mat'],'lmax')

end

