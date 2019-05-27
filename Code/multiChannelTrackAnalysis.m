function [multiChannelStruct] =  multiChannelTrackAnalysis(tracksFinal,movieInfoL,movieInfoO,fileListEdge,fileD,fileListCenter)
%MULTICHANNELTRACKANALYSIS collects and stores track measurements from different channels in single structure
%SYNOPSIS: [multiChannelStruct] =  multiChannelTrackAnalysis(tracksFinal,movieInfoL,movieInfoO,fileListEdge,fileD,fileListCenter)
%
%INPUT:
%
%       tracksFinal:tracking data (ex.tracksFinal)
%
%       movieInfoL: saved output from function mapIntensity for LAT channel
%       (i.e channel that was used for detection)
%
%       movieInfoO:saved output from function mapIntensity for other
%       channel (Grb2 or Nck)
%
%       fileListEdge: cell structure containing list of names for synapse mask of
%       each frame from time-lapse
%
%       fileD:Name of MD file (ex. "J2 1")
%
%       fileListCenter:cell structure containing list of names for cSMAC mask of
%       each frame from time-lapse
%
%OUTPUT:
%
%       multiChannelStruct: structure output containing following fields
%       for each track-
%
%     .rawIntensityL = raw intensity values from LAT channel
%     .rawIntensityO = raw intensity values from other channel (Grb2 or Nck)
%     .stdBackground = standard deviation of background
%     .latIntensity = intensity normalized by the first time point for LAT
%     .otherIntensity = intensity normalized by the first time point for
%     other channel
%     .normDistance = normalized distances of detected positions in track
%     .realDistance = real distance of detected positions in track
%     .normCSMACDistance = normalized distances of detected positions from
%     cSMAC
%     .realCSMACDistance = real distance of detected positions in track
%     from cSMAC
%     .cellEdge = full distance from cell edge to center
%     .frames = frames that a track exists in
%     .initialInt = initial intensity
%
%
%Tony Vega 2018
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
warning('off')

%Change format of tracking data
T = TracksHandle(tracksFinal);

%Go through each track
for m = 1:size(T,1) %Track

    for k = 1:length(T(m).f) %track position
        if T(m).f(k) == 0 %If there is a gap, assume track has disappeared
                          %And thus gets no value
            testL(k)=0;
            testO(k)=0;
            testNoise(k)=0;
            nDist(k)=0;
            nCDist(k) = 0;
            rDist(k)=0;
            rCDist(k)=0;
            cellEdge(k) = 0;
            centEdge(k) = 0;
            
        else %If track exists, get normalized distance of position
        
        %load frame masks
            maskingFile = imread([fileD '/Cell Mask/' fileListEdge{T(m).f(k)}]);
            centerFile = imread([fileD '/Center Mask/' fileListCenter{T(m).f(k)}]);
            
        %Find detections in QP that lie inside maskList
            [row, col] = find(maskingFile); 
            maskList = [row, col];
            QP = [T(m).y(k), T(m).x(k)];
            lia1 = ismember(round(QP),maskList,'rows');
            
            if lia1 %If position within mask, get normalized distance and record
                    % intensities in each channel
                
                %Get LAT Intensity through movieInfo
                    testL(k) = movieInfoL(T(m).f(k)).ampCorr(T(m).tracksFeatIndxCG(k));
                %Get Other Intensity through movieInfo
                    testO(k) = movieInfoO(T(m).f(k)).ampCorr(T(m).tracksFeatIndxCG(k));
                % Get std background
                    testNoise(k) = movieInfoO(T(m).f(k)).stdBG(T(m).tracksFeatIndxCG(k));
                %Get Norm Distance
                    indexQP = sub2ind(size(maskingFile),round(QP(:,1)),round(QP(:,2)));
                    erMask = bwmorph(maskingFile,'erode');
                %Get cell edge pixels
                    borderMask = maskingFile - erMask;
                %Get center(cSMAC) edge pixels
                    erMaskC = bwmorph(centerFile,'erode');
                    centerMask = centerFile-erMaskC;
                %Get center of cell
                    stats = regionprops(centerFile,'Centroid'); %Make weighted perhaps
                %---------------------------------------------------
                %Make center of cell origin and switch to polar
                %coordinates
                    [y,x] = find(borderMask);
                    yCorr = y-stats.Centroid(2);
                    xCorr = x-stats.Centroid(1);
                    [theta,rho] = cart2pol(xCorr,yCorr);

                %Bin data to get smoother signal
                    thetaT = theta;
                    rhoT = rho;
                    list = linspace(min(theta),max(theta));
                    for bN = 1:length(list)
                        ind = thetaT <= list(bN);
                        bRho(bN) = mean(rhoT(ind));
                        thetaT(ind) = [];
                        rhoT(ind) = [];
                    end
                %Get fit from smoothed signal
                    pp = spline(list,bRho);
                    %--------------------------------------------------
                %Do the same for the edge of cSMAC
                %Make center of cell origin and switch to polar
                %coordinates
                    [y,x] = find(centerMask);
                    yCorr = y-stats.Centroid(2);
                    xCorr = x-stats.Centroid(1);
                    [theta,rho] = cart2pol(xCorr,yCorr);

                %Bin data to get smoother signal
                    thetaT = theta;
                    rhoT = rho;
                    list = linspace(min(theta),max(theta));
                    for bN = 1:length(list)
                        ind = thetaT <= list(bN);
                        bRho(bN) = mean(rhoT(ind));
                        thetaT(ind) = [];
                        rhoT(ind) = [];
                    end
                    
                %Get fit from smoothed signal
                    ppC = spline(list,bRho);
                    %-------------------------------------------------
                %Get track position
                    [thetaP,rhoP] = cart2pol(round(QP(:,2))-stats.Centroid(1),round(QP(:,1))-stats.Centroid(2));
                    fullDist = ppval(pp,thetaP);
                    centDist = ppval(ppC,thetaP);
                    nCDist(k) = (fullDist-rhoP)./(fullDist-centDist);
                    nDist(k) = (fullDist-rhoP)./(fullDist);
                    rCDist(k) = fullDist-(rhoP-centDist);
                    rDist(k) = fullDist-(rhoP);
                    cellEdge(k) = fullDist;
                    

            else %If detection not within mask, set values to NaN
                testL(k)=NaN;
                testO(k)=NaN;
                nDist(k)=NaN;
                rDist(k)=NaN;
                nCDist(k)=NaN;
                rCDist(k)=NaN;
                cellEdge(k) = NaN;
                testNoise(k) = NaN;
% %                 centEdge(k) = NaN;
            end
        end
    end
    multiChannelStruct(m).rawIntensityL = testL;
    multiChannelStruct(m).rawIntensityO = testO;
    multiChannelStruct(m).stdBackground = testNoise;
    testLN = testL./(testL(1)); 
    testON = testO./(testO(1));

    multiChannelStruct(m).latIntensity = testLN;
    multiChannelStruct(m).otherIntensity = testON;
    multiChannelStruct(m).normDistance = nDist;
    multiChannelStruct(m).realDistance = rDist;
    multiChannelStruct(m).normCSMACDistance = nCDist;
    multiChannelStruct(m).realCSMACDistance = rCDist;
    multiChannelStruct(m).cellEdge = cellEdge;
    multiChannelStruct(m).frames = T(m).f;
    multiChannelStruct(m).initialInt = movieInfoO(T(m).f(1)).ampCorr(T(m).tracksFeatIndxCG(1));

    
    clear testL testO nDist rDist nCDist rCDist cellEdge testNoise 
end
end

