function [trackInfo] = multiChannelPlotTime(cond,typeM,controlM,plotM,boxM)
%MULTICHANNELPLOTTIME plots results from multi-channel track analysis
%SYNOPSIS:[trackInfo] = multiChannelPlotTime(cond,typeM,controlM,plotM,boxM)
%
%Note: should be run in folder where analyses are
%
%INPUT:
%
%       cond: Name of condition being analysed, refer to switch condition
%       below for proper names
%
%       typeM: 1 if Nck is being analyzed, 2 if Grb2. This is because the
%       LAT channel number is different for these conditions
%
%       controlM: 1 if looking at control tracks, zero otherwise
%
%       plotM: 1 to create new plot, 2 to use existing plot (no figure command)
%       zero for no plot
%
%       boxM: 1 to create new plot, 2 to use existing plot(no figure command)
%       zero for no plot
%
%OUTPUT: Plots
%
%       trackInfo: structure output of different track features before and
%       after actin engagement
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

%Based on condition, appropriate movie list is used 
switch cond
case 'Grb'
    listM = [3,6,10:15,18,19,22];%27,28
case 'Nck'
    listM = [3:10,12:18,20:29];
case  'DB'
    listM = [7:10,12:14,16,18,24,25,27];
case 'SMIFH'
    listM = [3,5,7:11,13,16:19,21,23];
case 'DMSO'
    listM = [3,13,15:17,18,20:24];
end

%Get names of movies
y = dir;
y1 = y(find(cellfun(@isdir,{y(:).name})));

%Initialize variables
testLatDist = []; %Cluster radial position
testR = []; %Track deviations from straight line
testPos = []; %Max deviation position within track
testD = [];%Max deviation position within cell
testGlobalTime = []; %Time point attributed to each cluster position
actinCheck = 0; %Use default actin distance threshold (0), or determine(1)
testIdx= []; %time point of cluster after actin engagement
testIntLat = []; %Lat intensity
testIntOther = [];%Other labelled molecule intensity


trackNum = 1;

%Go through list of movies
for m = listM
     %Load analysis from muliChannelTrackAnalysis
     load([y1(m).name '/multiChannelAnalysis1.mat'])   

     %Use results from actinPoisition to get distance threshold for each cell or use
     %default of 0.4 for every cell
    if actinCheck == 1
        %Attempt to get actin threshold
        load([y1(m).name '/lineFrame.mat'])
        avgLineFrame = NaN(length(lineFrame),11);
        for k = 1:length(lineFrame)
            avgLineFrame(k,:) = nanmean(lineFrame{k}(:,:));
        end
        finalLine = nanmean(avgLineFrame);
        [~,pos] = max(finalLine);
        actinThresh = pos/10;
%         threshDist = [threshDist; actinThresh];%Not used currently but
%         can look distribution of thresholds
    else
        actinThresh = 0.4;
    end
    
    %Load correct data whether looking at Nck or Grb2
    if typeM ==1
        multiChannelStruct = multiChannelStructNck;

    else
        multiChannelStruct = multiChannelStructGrb;
    end
    
    
    %% Filter tracks
    
    %Track length
    %Only use tracks of sufficient length
    lt = getTrackSEL(tracksFinal);
    lifeInd = find(lt(:,3)>=5);
    
    %Track asymmetry
    %go over all of these segments and determine if they are
    %asymmetric
    probDim=2;
    alphaAsym = 0.1;
    asymInd = NaN(size(lt,1),1);
    for k = lifeInd'
        startPoint = lt(k,1);
        endPoint  = lt(k,2);
        
        %get the particle positions along the track
        coordX = (tracksFinal(k).tracksCoordAmpCG(8*(startPoint-1)+1:8:8*(endPoint-(startPoint-1))))';
        coordY = (tracksFinal(k).tracksCoordAmpCG(8*(startPoint-1)+2:8:8*(endPoint-(startPoint-1))))';
        coordXY = [coordX coordY];
        
        %determine whether the track is sufficiently asymmetric
        [~,asymFlag] = asymDeterm2D3D(coordXY(:,1:probDim),alphaAsym);
        asymInd(k) = asymFlag;
    end
    
    %For a control look at symmetric tracks
    if controlM ==1
        distInd = find(asymInd==0);
    else
        distInd = find(asymInd==1);
    end
    
    % Intitial intensity (relevant for Nck)
    checkInd = zeros(size(lt,1),1);
    for k = distInd'
%     Is there any Nck initially?
         start= nanmean(multiChannelStruct(k).rawIntensityO(1:3));
         startNoise= nanmean(multiChannelStruct(k).stdBackground(1:3));
         
        if start > startNoise 
            checkInd(k) = 1;
        end
    end

    
    
    %Filter by distance travelled
    finalIndex = find(checkInd)';
    finalInd = zeros(size(lt,1),1);
    controlInd= zeros(size(lt,1),1);
    
    for k = finalIndex
    % Where does the track start and end?
        birth = multiChannelStruct(k).normDistance(1);
        death = multiChannelStruct(k).normDistance(end);

        if birth < actinThresh && death >= actinThresh
            %If track starts before threshold and ends after threshold,
            %keep index and store some data
            finalInd(k) = 1;
            multiChannelStruct(k).lifeTime = 1:length(multiChannelStruct(k).normDistance);
        elseif sum(multiChannelStruct(k).normDistance > actinThresh) == 0 && ~isnan(multiChannelStruct(k).latIntensity(1)) %birth > actinThresh && death > actinThresh
            %If track never crosses threshold, store as control
            controlInd(k) = 1;
            multiChannelStruct(k).lifeTime = 1:length(multiChannelStruct(k).normDistance);
        elseif sum(multiChannelStruct(k).normDistance < actinThresh) == 0 && ~isnan(multiChannelStruct(k).latIntensity(1))%birth < actinThresh && death < actinThresh
        %Currently ignore any tracks that start after threshold
%             controlInd(k) = 1;
%             multiChannelStruct(k).lifeTime = 1:length(multiChannelStruct(k).normDistance);
        end
        
    end
    
    %Keep relevant tracks
    if controlM ==1
        finalIndex = find(controlInd)';
    else
        finalIndex = find(finalInd)';
    end
    test = multiChannelStruct(finalIndex);
     
    %% Analysis
    if ~isempty(finalIndex) %&& controlM ~=1
        
        %go through each one
        for rs = 1:length(finalIndex)
            
        %Track coordinates
            x = tracksFinal(finalIndex(rs)).tracksCoordAmpCG(1:8:end)';
            y = tracksFinal(finalIndex(rs)).tracksCoordAmpCG(2:8:end)';
            
        %Initial track info
            intTrackLat = test(rs).rawIntensityL;
            intTrackOther = test(rs).rawIntensityO;
            
        %Track distances
            distTrackP = test(rs).normDistance;
            xT=x;
            yT=y;
            
        %update to only include points inside mask
        %Check that this means what you think it does
        %currently allowing all positions
        %{
            xT = x(test(rs).normCSMACDistance<1);
            yT = y(test(rs).normCSMACDistance<1);
            if isempty(xT)
                continue
            end
        % Only include positions outside of the cSMAC
            distTrack = distTrack(test(rs).normCSMACDistance<1);
            distTrackP = distTrackP(test(rs).normCSMACDistance<1);
            intTrackLat = intTrackLat(test(rs).normCSMACDistance<1);
            intTrackOther = intTrackOther(test(rs).normCSMACDistance<1);
            %}
        %% Find actin engagement
        %Set initial position to origin
            xCorr = xT - xT(1);
            yCorr = yT - yT(1);
            
        %Conver to polar coordinates, easier this way
            [th,rho] = cart2pol(xCorr,yCorr);
            
        %We want to put the last point on the x-axis
            theta = th(end)*(180/pi);
            tform = affine2d([cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; 0 0 1]);
            [xx,yy] = transformPointsForward(tform,xCorr,yCorr);
            
        %At this point take out all NaN's
            distTrackP = distTrackP(~isnan(yy));
            intTrackLat = intTrackLat(~isnan(yy));
            intTrackOther = intTrackOther(~isnan(yy));
            yy = yy(~isnan(yy));
            xx = xx(~isnan(xx));
            
            
            xCorr = xx - xx(end);
            yCorr = yy;
            [th,rho] = cart2pol(xCorr,yCorr);
            
        %New method to get angles between positions
            vector2AvgX = xCorr(2:end);
            vector2AvgY = yCorr(2:end);
            vector1AvgY = yCorr(1:end-1);
            vector1AvgX = xCorr(1:end-1);
            
            top = ((vector1AvgX.*vector2AvgX)+(vector1AvgY.*vector2AvgY));
            bottom = (sqrt((vector1AvgX).^2+(vector1AvgY).^2)).*(sqrt((vector2AvgX).^2+(vector2AvgY).^2));
            angle = acos(top./bottom);
            angle(angle==0) = NaN;
            dTheta = diff(th(~isnan(th)));
            
        %Get change in distances between adjacent positions 
            dRho = diff(rho(~isnan(rho)));
        %Finally get ratio of these changes to find "actin engagement"
            dRatio = abs(dRho)./abs(angle);
            level = prctile(dRatio,90);
            peaks = find(dRatio>=level);
            checkPeak = distTrackP(peaks)<=actinThresh;
            if sum(checkPeak)>0
            else
                level = prctile(dRatio(~isnan(dRatio)),75);
                peaks = find(dRatio>=level);
%                 checkPeak = distTrackP(peaks)<=(actinThresh);
%                 if sum(checkPeak)>0
%                 else
% %                     peaks=1;
%                 end
            end

        %Now correct everything by actin engagement
            newStart = peaks(1);
            trackIdx = (1:length(yy))-newStart;
            lifeTime = (1:length(yy));
            %% Comparing cluster properties before and after actin engagement
            %This should potentially be amended to look at actin threshold
            %instead...
            if newStart >0

                    
%                     testRawIntLat = [testRawIntLat (intTrackLat(1:newStart))];
%                     testRawIntOther = [testRawIntOther (intTrackOther(1:newStart))];
                testSpeedBefore =  (distTrackP(newStart)-distTrackP(1))./length(distTrackP(1:newStart));
                testSpeedAfter =  (distTrackP(end)-distTrackP(newStart))./(length(distTrackP(newStart:end))-1);

%                 latB=(intTrackLat(1:newStart));%./intTrackLat(1));
%                 latA=(intTrackLat(newStart:end));%./intTrackLat(newStart));
%                 otherB=(intTrackOther(1:newStart));
%                 otherA=(intTrackOther(newStart:end));
%                 ratioB=(intTrackOther(1:newStart))./latB;
%                 ratioA=(intTrackOther(newStart:end))./latA;

                %Fill in track information
                trackInfo(trackNum).latIntensity = [intTrackLat(1),intTrackLat(newStart),intTrackLat(end)];
                trackInfo(trackNum).otherIntensity = [intTrackOther(1),intTrackOther(newStart),intTrackOther(end)];
                trackInfo(trackNum).timeToEngage = newStart;
                trackInfo(trackNum).timeToEnd = length(intTrackLat)-newStart;
                trackInfo(trackNum).trackSpeedB= testSpeedBefore;
                trackInfo(trackNum).trackSpeedA= testSpeedAfter;
                trackInfo(trackNum).latIntChangeB = (intTrackLat(newStart)-intTrackLat(1))./(length(intTrackLat(1:newStart))-1);
                trackInfo(trackNum).latIntChangeA = (intTrackLat(end)-intTrackLat(newStart))./(length(intTrackLat(newStart:end))-1);
                trackInfo(trackNum).varianceLatB= var(intTrackLat(1:newStart));
                trackInfo(trackNum).varianceLatA= var(intTrackLat(newStart:end));
                trackInfo(trackNum).varianceOtherB= var(intTrackOther(1:newStart));
                trackInfo(trackNum).varianceOtherA= var(intTrackOther(newStart:end));
                trackNum = trackNum+1;


            end
            %% Storing everything for plotting
            intTrackLat = intTrackLat./mean(intTrackLat(1:3));%indexAvg:newStart
            intTrackOther = intTrackOther./mean(intTrackOther(1:3));
            testIdx= [testIdx trackIdx];%
            testIntLat = [testIntLat intTrackLat];%
            testIntOther = [testIntOther intTrackOther];%
            testLatDist = [testLatDist distTrackP];%
            testGlobalTime= [testGlobalTime lifeTime];%

                    

            
            %% Measure track deviation from straight line
            %Here we place the segment of the track after actin engagement
            %along the x-axis and use y values to assess deviation from a
            %straight line
            distTrackP = distTrackP(newStart:end);
            nX = xCorr(newStart:end);
            nY = yCorr(newStart:end);
            yCorrN = nY-nY(1);
            xCorrN = nX-nX(1);
            [th,rho] = cart2pol(xCorrN,yCorrN);
            theta = th(end)*(180/pi);
            tform = affine2d([cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; 0 0 1]);
            [xxN,yyN] = transformPointsForward(tform,xCorrN,yCorrN);
            yy= yyN;
            % Make majority of data points positive so all major deviations
            % are positive
            indMaj = find(yy<0);
            fracNeg = length(indMaj)./length(yy);
            
            if fracNeg > 0.5
                yy = -yy;
            end
            
            %Collect the deviations, max deviation (in track and in cell)
            %Note: not all are currently being used. Interchangeable in
            %boxM if statement
            [q,w] = max(abs(yy));
            testR = [testR yy']; %All deviations
            testPos = [testPos w./length(yy)];%Where is max in track?
            testD = [testD 1-distTrackP(1)];%Where is max in cell?
            %}

        end
    end

end

%% Plotting deviation from straight line
if boxM ==1
    figure;
end

%----------------------------------------------------
if boxM >0
    hold on
    edges = -10:20;
    h = histogram(testR,edges,'Normalization','probability');
    histP(:,1) = h.BinEdges(1:end-1);
    histP(:,2) = h.Values;
else
    histP = [];
end



%% Plotting change in intensity across synapse
%Here we are plotting a minimal boxplot by collecting the median, upper,
%and lower bounds of the data at different distances and times
%Set negative intensity values to zero
indNeg = find(testIntOther<0);
testIntOther(indNeg) =0;

if plotM == 1
    figure;
    hold on
    %Make three different plots
    for f = 1:3
        latInfo = testIntLat;
        nckInfo = testIntOther;
        switch f
            case 1 %Plot according to distance
                latInfo(2,:) = testLatDist;
                nckInfo(2,:) = testLatDist;
                tickC = 0.1:0.1:1;
            case 2 %Plot according to time
                latInfo(2,:) = testGlobalTime;
                nckInfo(2,:) = testGlobalTime;
                tickC = min(nckInfo(2,:)):max(nckInfo(2,:));
            case 3 %Plot according to time from radial start
                latInfo(2,:) = testIdx;
                nckInfo(2,:) = testIdx;
                tickC = min(nckInfo(2,:)):max(nckInfo(2,:));
        end
        %Get intensity from all track positions
        indN = isnan(latInfo(1,:));
        latInfo(:,indN) = [];
        nckInfo(:,indN) = [];
        
        %Looking at ration of Grb2 or Nck to LAT
        nckInfoT = nckInfo;
        nckInfoT(1,:) = nckInfoT(1,:)./latInfo(1,:);
        
        %Initialize some variables
        medNck = NaN(1,length(tickC));
        avgNck = NaN(1,length(tickC));
        nUpper = NaN(1,length(tickC));
        nLower = NaN(1,length(tickC));
        count = 1;

        hold on
        subplot(1,3,f)
        hold on
        
        % Go through bins and collect data
        for m =tickC 
            ind = find(nckInfoT(2,:)<=m);
            indOrg = ind;

            if ~isempty(ind) && length(ind)>=10
                if m == 0.2
                    firstDist = nckInfoT(1,ind);
                elseif m > 0.2 
                   newDist = nckInfoT(1,ind);
                    p = ranksum(firstDist,newDist)
                    m
                end
                medNck(1,count) = nanmedian(nckInfoT(1,ind)); 
                avgNck(1,count) = nanmean(nckInfoT(1,ind));
                nUpper(1,count) = prctile(nckInfoT(1,ind),50)-(1.57*(prctile(nckInfoT(1,ind),75)-prctile(nckInfoT(1,ind),25))./sqrt(length(ind)));
                nLower(1,count) = prctile(nckInfoT(1,ind),50)+(1.57*(prctile(nckInfoT(1,ind),75)-prctile(nckInfoT(1,ind),25))./sqrt(length(ind)));
                plot([m,m],[nLower(1,count),nUpper(1,count)],'Color',[0    0.4470    0.7410])
                nckInfoT(:,indOrg) = [];%[0.9290    0.6940    0.1250]
            else
                medNck(1,count) = NaN;
                avgNck(1,count) = NaN;
                nUpper(1,count) = NaN;
                nLower(1,count) = NaN;

                nckInfoT(:,indOrg) = [];
            end
            count = count+1;
        end
        scatter(tickC,medNck,'o','LineWidth',2,'MarkerEdgeColor','none','MarkerFaceColor',[0    0.4470    0.7410]);
        scatter(tickC,nUpper,'^','LineWidth',2,'MarkerEdgeColor','none','MarkerFaceColor',[0    0.4470    0.7410]);
        scatter(tickC,nLower,'v','LineWidth',2,'MarkerEdgeColor','none','MarkerFaceColor',[0    0.4470    0.7410]);
        
        %Change labelling and axes based on figure
        switch f
            case 1 
                set(gca,'XTick',0.05:0.1:1.15)
                set(gca,'XTickLabel',{'1','0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1',})
                xlim([0 1.1])
                xlabel('Normalized Radial position')
                ylabel('Intensity Ratio')
            case 2
                xlabel('Time from track start')
            case 3
                xlabel('Time from radial movement start')
        end
    end

elseif plotM ==2
    
    hold on

    for f = 1:3
        latInfo = testIntLat;
        nckInfo = testIntOther;
        switch f
            case 1 
                latInfo(2,:) = testLatDist;
                nckInfo(2,:) = testLatDist;
                tickC = 0.1:0.1:1;
            case 2
                latInfo(2,:) = testGlobalTime;
                nckInfo(2,:) = testGlobalTime;
                tickC = min(nckInfo(2,:)):max(nckInfo(2,:));
            case 3
                latInfo(2,:) = testIdx;
                nckInfo(2,:) = testIdx;
                tickC = min(nckInfo(2,:)):max(nckInfo(2,:));
        end
        indN = isnan(latInfo(1,:));
        latInfo(:,indN) = [];
        nckInfo(:,indN) = [];

        nckInfoT = nckInfo;
        nckInfoT(1,:) = nckInfoT(1,:)./latInfo(1,:);

        medNck = NaN(1,length(tickC));
        avgNck = NaN(1,length(tickC));
        stdNck = NaN(1,length(tickC));
        nUpper = NaN(1,length(tickC));
        nLower = NaN(1,length(tickC));
        count = 1;
        % figure; 
        hold on
        subplot(1,3,f)
        hold on
        for m =tickC %0.1:0.1:1
            ind = find(nckInfoT(2,:)<=m);
%             [~,inlierIdx] = detectOutliers(nckInfoT(1,ind),3);
            indOrg = ind;
%             ind= ind(inlierIdx);
            if ~isempty(ind) && length(ind)>=10
        %         boxplot(nckInfoT(1,ind),'plotstyle','compact','notch','on','positions',m*10,'widths',0.5)
                medNck(1,count) = nanmedian(nckInfoT(1,ind)); %medNck(1,round(m*10))
                avgNck(1,count) = nanmean(nckInfoT(1,ind));
            %     stdNck(1,round(m*10)) = nanstd(nckInfoT(1,ind));
                nUpper(1,count) = prctile(nckInfoT(1,ind),50)-(1.57*(prctile(nckInfoT(1,ind),75)-prctile(nckInfoT(1,ind),25))./sqrt(length(ind)));
                nLower(1,count) = prctile(nckInfoT(1,ind),50)+(1.57*(prctile(nckInfoT(1,ind),75)-prctile(nckInfoT(1,ind),25))./sqrt(length(ind)));
                if f >1
                    plot([m+0.3,m+0.3],[nLower(1,count),nUpper(1,count)],'Color',[0.9290    0.6940    0.1250])
                else
                    plot([m+0.025,m+0.025],[nLower(1,count),nUpper(1,count)],'Color',[0.9290    0.6940    0.1250])
                end
            nckInfoT(:,indOrg) = [];%[0.9290    0.6940    0.1250]
            else
                medNck(1,count) = NaN;
                avgNck(1,count) = NaN;
            %     stdNck(1,round(m*10)) = nanstd(nckInfoT(1,ind));
                nUpper(1,count) = NaN;
                nLower(1,count) = NaN;

                nckInfoT(:,indOrg) = [];
            end
            count = count+1;
        end
        if f >1
            scatter(tickC+0.3,medNck,'o','LineWidth',2,'MarkerEdgeColor','none','MarkerFaceColor',[0.9290    0.6940    0.1250]);
            scatter(tickC+0.3,nUpper,'^','LineWidth',2,'MarkerEdgeColor','none','MarkerFaceColor',[0.9290    0.6940    0.1250]);
            scatter(tickC+0.3,nLower,'v','LineWidth',2,'MarkerEdgeColor','none','MarkerFaceColor',[0.9290    0.6940    0.1250]);
        else
            scatter(tickC+0.025,medNck,'o','LineWidth',2,'MarkerEdgeColor','none','MarkerFaceColor',[0.9290    0.6940    0.1250]);
            scatter(tickC+0.025,nUpper,'^','LineWidth',2,'MarkerEdgeColor','none','MarkerFaceColor',[0.9290    0.6940    0.1250]);
            scatter(tickC+0.025,nLower,'v','LineWidth',2,'MarkerEdgeColor','none','MarkerFaceColor',[0.9290    0.6940    0.1250]);
        end
        %}
        ylim([-0.5 2.5])
        switch f
            case 1 
                set(gca,'XTick',0.065:0.1:1.15)
                set(gca,'XTickLabel',{'1','0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1',})
                xlim([0 1.1])
                xlabel('Normalized Radial position')
                ylabel('Intensity Ratio')
            case 2
                xlabel('Time from track start')
            case 3
                xlabel('Time from radial movement start')
        end
    end
%-----------------------------------------------------------------------------------------
%NOT IMPLEMENTED CURRENTLY
% Similar plot to plotM =1 but boxplot is colored by actin engagement
elseif plotM==3
    %{
    figure;
    hold on

    for f = 1:3
        latInfo = testIntOther;
        nckInfo = testIntOther;
        
        lat2Info = testIntLat;
        latPInfo = testIntLat;
        latNInfo = testIntLat;
        distInfo = testLatDist;
        
        indN = isnan(latInfo(1,:));
        latInfo(:,indN) = [];
        nckInfo(:,indN) = [];
        distInfo(:,indN) = [];
        lat2Info(:,indN) = [];
        latPInfo(:,indN) = [];
        latNInfo(:,indN) = [];        

        switch f
            case 1 
                latInfo(2,:) = distInfo;
                nckInfo(2,:) = distInfo;
%                 nckInfo(1,:) = nckInfo(1,:)./latInfo(1,:);
                indP= find(testIdx>=0);
                nckInfo = nckInfo(:,indP);
                latPInfo = latPInfo(:,indP); 
%                 latInfo = nckInfo;
%                 latInfo(1,:) = nckInfo(1,:)./latInfo(1,:);
                indN= find(testIdx<0);
                latInfo = latInfo(:,indN);
                latNInfo = latNInfo(:,indN);
                tickC = 0.1:0.1:1;
                nckInfoT = nckInfo;
                nckInfoT(1,:) = nckInfo(1,:)./latPInfo(1,:);
                latInfoT = latInfo;
                latInfoT(1,:) = latInfo(1,:)./latNInfo(1,:);
                indEr = find(latInfoT(2,:)>0.4);
                latInfoT(:,indEr) = []; 
            case 2
                latInfo(2,:) = testGlobalTime;
                nckInfo(2,:) = testGlobalTime;
                tickC = min(nckInfo(2,:)):max(nckInfo(2,:));
                nckInfoT = nckInfo;
                nckInfoT(1,:) = nckInfo(1,:)./lat2Info(1,:);
                latInfoT = latInfo;
                latInfoT(1,:) = latInfo(1,:)./lat2Info(1,:);
            case 3
                latInfo(2,:) = testIdx;
                nckInfo(2,:) = testIdx;
                tickC = min(nckInfo(2,:)):max(nckInfo(2,:));
                nckInfoT = nckInfo;
                nckInfoT(1,:) = nckInfo(1,:)./lat2Info(1,:);
                latInfoT = latInfo;
                latInfoT(1,:) = latInfo(1,:)./lat2Info(1,:);
                indP= find(testIdx>=0);
                indN= find(testIdx<0);
                latInfoT= latInfoT(:,indN);
                nckInfoT= nckInfoT(:,indP);
        end


        medNck = NaN(1,length(tickC));
        avgNck = NaN(1,length(tickC));

        nUpper = NaN(1,length(tickC));
        nLower = NaN(1,length(tickC));
        
        medLat = NaN(1,length(tickC));
        avgLat = NaN(1,length(tickC));
        nUpperL = NaN(1,length(tickC));
        nLowerL = NaN(1,length(tickC));
        count = 1;
        figure; 
        hold on
        subplot(1,3,f)
        hold on
        for m =tickC %0.1:0.1:1
            ind = find(nckInfoT(2,:)<=m);
%             [~,inlierIdx] = detectOutliers(nckInfoT(1,ind),3);
            indOrg = ind;
%             ind= ind(inlierIdx);
            if ~isempty(ind) && length(ind)>=10

                medNck(1,count) = nanmedian(nckInfoT(1,ind)); %medNck(1,round(m*10))
                avgNck(1,count) = nanmean(nckInfoT(1,ind));
                
                nUpper(1,count) = prctile(nckInfoT(1,ind),50)-(1.57*(prctile(nckInfoT(1,ind),75)-prctile(nckInfoT(1,ind),25))./sqrt(length(ind)));
                nLower(1,count) = prctile(nckInfoT(1,ind),50)+(1.57*(prctile(nckInfoT(1,ind),75)-prctile(nckInfoT(1,ind),25))./sqrt(length(ind)));
               
                if f ==3 
                    plot([m,m],[nLower(1,count),nUpper(1,count)],'Color',[0    0.4470    0.7410])

                elseif f == 2
                    plot([m,m],[nLower(1,count),nUpper(1,count)],'Color',[0.8500    0.3250    0.0980])
                else
                    plot([m+0.015,m+0.015],[nLower(1,count),nUpper(1,count)],'Color',[0    0.4470    0.7410])
                end

                nckInfoT(:,indOrg) = [];%[0.9290    0.6940    0.1250]
            else
                medNck(1,count) = NaN;
                avgNck(1,count) = NaN;
                nUpper(1,count) = NaN;
                nLower(1,count) = NaN;

                nckInfoT(:,indOrg) = [];
                
            end
            ind = find(latInfoT(2,:)<=m);
%             [~,inlierIdx] = detectOutliers(latInfoT(1,ind),3);
            indOrg = ind;
%             ind= ind(inlierIdx);
            if ~isempty(ind) && length(ind)>=10

                medLat(1,count) = nanmedian(latInfoT(1,ind)); %medNck(1,round(m*10))
                avgLat(1,count) = nanmean(latInfoT(1,ind));

                nUpperL(1,count) = prctile(latInfoT(1,ind),50)-(1.57*(prctile(latInfoT(1,ind),75)-prctile(latInfoT(1,ind),25))./sqrt(length(ind)));
                nLowerL(1,count) = prctile(latInfoT(1,ind),50)+(1.57*(prctile(latInfoT(1,ind),75)-prctile(latInfoT(1,ind),25))./sqrt(length(ind)));                
                if f ==3 
                    plot([m,m],[nLowerL(1,count),nUpperL(1,count)],'Color',[0.9290    0.6940    0.1250])
                elseif f == 2
                    
                else
                    plot([m-0.015,m-0.015],[nLowerL(1,count),nUpperL(1,count)],'Color',[0.9290    0.6940    0.1250])
                end
                latInfoT(:,indOrg) =[];
            else
                
                medLat(1,count) = NaN;
                avgLat(1,count) = NaN;
                nUpperL(1,count) = NaN;
                nLowerL(1,count) = NaN;

                latInfoT(:,indOrg) =[];
            end
            count = count+1;
        end
        if f ==3 
            scatter(tickC,medNck,'o','LineWidth',2,'MarkerEdgeColor','none','MarkerFaceColor',[0    0.4470    0.7410]);
            scatter(tickC,nUpper,'^','LineWidth',2,'MarkerEdgeColor','none','MarkerFaceColor',[0    0.4470    0.7410]);
            scatter(tickC,nLower,'v','LineWidth',2,'MarkerEdgeColor','none','MarkerFaceColor',[0    0.4470    0.7410]);

            scatter(tickC,medLat,'o','LineWidth',2,'MarkerEdgeColor','none','MarkerFaceColor',[0.9290    0.6940    0.1250]);
            scatter(tickC,nUpperL,'^','LineWidth',2,'MarkerEdgeColor','none','MarkerFaceColor',[0.9290    0.6940    0.1250]);
            scatter(tickC,nLowerL,'v','LineWidth',2,'MarkerEdgeColor','none','MarkerFaceColor',[0.9290    0.6940    0.1250]);
        elseif f ==2
            scatter(tickC,medLat,'o','LineWidth',2,'MarkerEdgeColor','none','MarkerFaceColor',[0.8500    0.3250    0.0980]);
            scatter(tickC,nUpperL,'^','LineWidth',2,'MarkerEdgeColor','none','MarkerFaceColor',[0.8500    0.3250    0.0980]);
            scatter(tickC,nLowerL,'v','LineWidth',2,'MarkerEdgeColor','none','MarkerFaceColor',[0.8500    0.3250    0.0980]);
        else
            scatter(tickC+0.015,medNck,'o','LineWidth',2,'MarkerEdgeColor','none','MarkerFaceColor',[0    0.4470    0.7410]);
            scatter(tickC+0.015,nUpper,'^','LineWidth',2,'MarkerEdgeColor','none','MarkerFaceColor',[0    0.4470    0.7410]);
            scatter(tickC+0.015,nLower,'v','LineWidth',2,'MarkerEdgeColor','none','MarkerFaceColor',[0    0.4470    0.7410]);

            scatter(tickC-0.015,medLat,'o','LineWidth',2,'MarkerEdgeColor','none','MarkerFaceColor',[0.9290    0.6940    0.1250]);
            scatter(tickC-0.015,nUpperL,'^','LineWidth',2,'MarkerEdgeColor','none','MarkerFaceColor',[0.9290    0.6940    0.1250]);
            scatter(tickC-0.015,nLowerL,'v','LineWidth',2,'MarkerEdgeColor','none','MarkerFaceColor',[0.9290    0.6940    0.1250]);
        end
        %
        ylim([-1 3])
        switch f
            case 1 
                set(gca,'XTick',0.05:0.1:1.15)
                set(gca,'XTickLabel',{'1','0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1',})
                xlim([0 1.1])
                xlabel('Normalized Radial position')
                ylabel('Intensity Ratio')
            case 2
                xlabel('Time from track start')
            case 3
                xlabel('Time from radial movement start')
        end
    end
    %}

end



end

