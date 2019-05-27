% Script to correct for drift found in steady state & contraction time-lapses
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

%% Initialization
%Enter location of time lapses
fileLocation = '/project/biophysics/shared/RosenJaqamanCollab/code4github/exampleData';

%Supply channel number for LAT
latChannel = 1;

%Supply folder location where corrected LAT time-lapse will be saved
latFolder = 'DriftCorrected';

%Choose whether actin channeal will also be corrected
actinCheck = 1;

%Supply channel number for actin
actinChannel = 2;

%Supply folder location where corrected actin time-lapse will be saved
actinFolder = 'DriftCorrected';

%Get list of time lapses in folder
filePre = dir([fileLocation,filesep,'*.tif']);

%% Part 1- Correct Drift in LAT Channel
for m = 1:size(filePre,1)
    
    %Create stack with only time frames from LAT Channel
    j=1;
    imgFile = [filePre(m).name];    
    imgInfo = imfinfo(imgFile);
    endTime = length(imgInfo);
    
    %Note: Assumption that there are two channels
    for k = latChannel:2:endTime
        image = imread(imgFile,k);
        imStack(:,:,j)=image;
        j =j+1;
    end


    %Create name for new time-lapse
     outputFolder = [latFolder filesep filePre(m).name];
     
    %Apply drift correction
     [imStackCorr,maxCorrL,imTformL] = driftCorrectImageStack(imStack);
    % Save new time-lapse 
    i=1;
    I = uint16(imStackCorr(:,:,i));
    imwrite(I,outputFolder)
    for i = 2:size(imStack,3)
        I = uint16(imStackCorr(:,:,i));
        imwrite(I,outputFolder,'WriteMode','append');
    end
    clear imStack

%% PART 2- Actin Correction
    if actinCheck ==1
    %To use correction from one channel on the other
        j=1;
        imgFile = [filePre(m).name ];
        imgInfo = imfinfo(imgFile);
        endTime = length(imgInfo);
        for k = actinChannel:2:endTime
            image = imread(imgFile,k);
            imStack(:,:,j)=image;
            j =j+1;
        end
        imStackCorr = imStack;%%this needs to be actin
        imSize = size(imStack);
        refInd=1:size(imStack,3)-1;
        corrInd=2:size(imStack,3);
        nCorr=numel(corrInd);
        for j = 1:nCorr

            currRef = double(imStackCorr(:,:,refInd(j)));    
            currCorr = double(imStack(:,:,corrInd(j)));

            %Get cross-orrelation between images to do an initial translation since imregister is slow (optimization
            %based)    

            [Xi,Yi] = meshgrid((1:imSize(2))+imSize(2)-maxCorrL(corrInd(j),2),(1:imSize(1))+imSize(1)-maxCorrL(corrInd(j),1));
            currCorr = interp2(currCorr,Xi,Yi,'nearest');%Lazy way to do shifting + fill NaNs
            currCorr(isnan(currCorr)) = 0;


                %We can't call imregister directly because we want to specify fill
                %values

                Rfixed = imref2d(size(currRef));
                currCorr = imwarp(currCorr,imTformL{corrInd(j)},'OutputView',Rfixed,'FillValues',NaN);        
                currCorr(isnan(currCorr)) = 0;


            imStackCorr(:,:,corrInd(j)) = currCorr;


        end

         outputFolder = [actinFolder filesep filePre(m).name];
        i=1;
        I = uint16(imStackCorr(:,:,i));
        imwrite(I,[outputFolder(1:end-4) 'DriftCorrActin.tif' ])
        for i = 2:size(imStack,3)
            I = uint16(imStackCorr(:,:,i));

            imwrite(I,[outputFolder(1:end-4) 'DriftCorrActin.tif' ],'WriteMode','append');
        end
        clear imStackCorr imStack
    end
 end
