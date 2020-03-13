function KiNeT( scores, dt )
%KiNeT() Runs basic KiNeT analysis from

% Remington, E. D., Narain, D., Hosseini, E. A., & Jazayeri, M. (2018).
% Flexible sensorimotor computations through rapid reconfiguration of
% cortical dynamics. Neuron, 98(5), 1005-1019.

% X should be an N (dimensions/neurons/factors) x C (conditions) x T (time
% points) matrix. Conditions should be ordered in increasing duration.
% Conditions for which the length is less than the max length should be
% padded with NaN values. The length of each condition should be the same
% for all N.
%
% E. Remington 2019

[nDim, nCond, ~] = size(scores);

% Select middle length condition
midIdx = round(mean(1:nCond));

for i = nCond:-1:1
    ntAll(i) = find( ...
        ~isnan(scores(1,i,:)),1,'last');
end

nt = ntAll(midIdx);
distancesAll = nan(nt,nCond);
minDistInd = nan(nt,nCond);
vectorAngleMarginal = nan(nt,nCond - 2);
normVectorMarginal = nan(nt,nDim,nCond - 1);
normIntraManifoldVector = nan(nt,nDim);
distancesMarginal = nan(nt,nCond-1);
spaceSimilarity = nan(nt,1);
spaceAngle = nan(nt,1);

cRgb100 = [0.75 0.75 0.75; 0.25 0.25 0.25];
lineColor = InterpColorMap(cRgb100,linspace(0,1,nCond));

for j = 1:nt
    %% Match time points by distance (Speed, Distance)
    for i = nCond:-1:1
        distances = squeeze(sqrt(sum((scores(:,midIdx,j) - ...
            scores(:,i,:)).^2)));
        [minDistTemp(i),minDistInd(j,i)] = min(distances);
        
        normVector(:,i) = (scores(:,midIdx,j) - ...
            scores(:,i,minDistInd(j,i)))/ ...
            norm((scores(:,midIdx,j) - ...
            scores(:,i,minDistInd(j,i))));
    end
    
    
    %% Vectors across trajectories and angles between vectors (Direction)
    for i = nCond:-1:2
        vectorMarginal = scores(:,i,minDistInd(j,i)) - ...
            scores(:,i-1,minDistInd(j,i-1));
        
        normVectorMarginal(j,:,i-1) = vectorMarginal./norm(vectorMarginal);
        
        distancesMarginal(j,i-1) = squeeze(sqrt(sum((scores(:,i,minDistInd(j,i)) - ...
            scores(:,i-1,minDistInd(j,i-1))).^2)));
    end
    
    normIntraManifoldVector(j,:) = mean(normVectorMarginal(j,:,:),3)/ ...
        norm(mean(normVectorMarginal(j,:,:),3));
    
    for i = nCond-1:-1:2
        vectorAngleMarginal(j,i-1) = acos(dot(normVectorMarginal(j,:,i), ...
            normVectorMarginal(j,:,i-1)));
    end
    
    
    %% Assign conditions as being closer to first or last (for Distance).
    %Irrelevant for three conditions. Tested for 5 or more.
    for i = nCond:-1:1
        angleVectorShort(i) = acos(dot(normVector(:,1), ...
            normVector(:,i)));
        angleVectorLong(i) = acos(dot(normVector(:,end), ...
            normVector(:,i)));
        
        if angleVectorShort(i) < angleVectorLong(i)
            distancesAll(j,i) = -minDistTemp(i);
        elseif angleVectorShort(i) > angleVectorLong(i)
            distancesAll(j,i) = minDistTemp(i);
        else
            distancesAll(j,i) = 0;
        end
    end
    
    %% Calculate subspace similarity relative to first time point
    for i = nCond:-1:1
        fullActivitySnap(i,:) = squeeze(scores(:,i,minDistInd(j,i)));
        originalFullSnap(i,:) = squeeze(scores(:,i,minDistInd(1,i)));
    end
    %Similarity index
    spaceSimilarity(j)  = ...
        DatasetSimilarity(fullActivitySnap,originalFullSnap);
    
    %Rotation index (uses first dimension of dataset)
    if j > 1
        spaceAngle(j) = acos(dot(normIntraManifoldVector(1,:), ...
            normIntraManifoldVector(j,:)));
    else
        spaceAngle(j) = 0;
    end
   
end
%% 'Speed', 'Direction', and 'Distance'
drawnow()
figure('position',[190   141   413   912]);

limIndX = find(any( ...
    ~isnan(scores(1,midIdx,:)),2),1,'last');
limIndY = find(any(any( ...
    ~isnan(scores),2),1),1,'last');
yLim = [0 limIndY*dt];

timeVect = dt:dt:nt*dt;
xLim = [0 limIndX*dt];
indVect = 1:nt;
xLabel = 'Ref. time from set (s)';
yAxisLocation = 'Left';
rampMult = 1;

subplot(2,1,1)
pos1 = get(gca,'Position');
set(gca,'Position',[pos1(1) pos1(2)+0.2*pos1(4) pos1(3) 0.8*pos1(4)])
for i = nCond:-1:1
    line(timeVect,minDistInd(indVect,i)*dt, ...
        'Color',lineColor(i,:))
    line(timeVect,rampMult*timeVect*(ntAll(i)/nt), ...
        'Color',lineColor(i,:), ...
        'LineStyle','--')
end

ylabel({'Speed: Time of','minimum distance (s)'})
set(gca,'XLim',xLim, ...
    'YLim',yLim, ...
    'YAxisLocation',yAxisLocation, ...
    'DataAspectRatio',[1 1 1], ...
    'XTick',[])
ar = get(gca,'PlotBoxAspectRatio');
pos1 = get(gca,'Position');

subplot(4,1,3)

meanVectorAngleMarginal = nanmean(vectorAngleMarginal');

plot(timeVect, ...
    meanVectorAngleMarginal(indVect), ...
    'Color',lineColor(midIdx,:));
hold on;

line([timeVect(1) timeVect(end)],[pi/2 pi/2], ...
    'Color','k', ...
    'LineStyle','--');

pos2 = get(gca,'Position');
ylabel({'Directoion:', 'Angle (rad)'})
set(gca,'XLim',xLim, ...
    'YLim',[0 pi], ...
    'YTick',pi/2, ...
    'YTickLabel','pi/2', ...
    'YTickLabelRotation',90, ...
    'YAxisLocation',yAxisLocation, ...
    'PlotBoxAspectRatio',[ar(1)*pos1(4)/pos2(4) ar(2) ar(3)], ...
    'XTick',[])

box off

subplot(4,1,4)

for i = nCond:-1:1
    line(timeVect,distancesAll(indVect,i), ...
        'Color',lineColor(i,:));
end

pos3 = get(gca,'Position');
ylabel({'Distance: from','ref. trajectory'})
set(gca,'XLim',xLim, ...
    'YAxisLocation',yAxisLocation, ...
    'PlotBoxAspectRatio',[ar(1)*pos1(4)/pos3(4) ar(2) ar(3)])


%% Similarity index and rotation angle
drawnow()
figure('position',[190   141   413   912]);

limIndX = find(any( ...
    ~isnan(scores(1,midIdx,:)),2),1,'last');
limIndY = find(any(any( ...
    ~isnan(scores),2),1),1,'last');
yLim = [0 limIndY*dt];

timeVect = dt:dt:nt*dt;
xLim = [0 limIndX*dt];
indVect = 1:nt;
xLabel = 'Ref. time from set (s)';
yAxisLocation = 'Left';

%subplot(2,1,1) just for scaling x similarity plots
subplot(2,1,1)
pos1 = get(gca,'Position');
set(gca,'Position',[pos1(1) pos1(2)+0.2*pos1(4) pos1(3) 0.8*pos1(4)])
line(timeVect,nan(size(timeVect)), ...
    'Color',lineColor(1,:), ...
    'LineStyle','--')

set(gca,'XLim',xLim, ...
    'YLim',yLim, ...
    'YAxisLocation',yAxisLocation, ...
    'DataAspectRatio',[1 1 1], ...
    'XTick',[])
ar = get(gca,'PlotBoxAspectRatio');
pos1 = get(gca,'Position');


subplot(4,1,3)
plot(timeVect,spaceAngle(indVect), ...
    'Color',lineColor(midIdx,:))

line([timeVect(1) timeVect(end)],[pi/2 pi/2], ...
    'Color','k', ...
    'LineStyle','--');

pos2 = get(gca,'Position');
ylabel({'Rotation index:','Angle wrt t = 0'})
set(gca,'XLim',xLim, ...
    'YLim',[0 pi], ...
    'YTick',pi/2, ...
    'YTickLabel','pi/2', ...
    'YTickLabelRotation',90, ...
    'YAxisLocation',yAxisLocation, ...
    'PlotBoxAspectRatio',[ar(1)*pos1(4)/pos2(4) ar(2) ar(3)], ...
    'XTick',[])

box off


subplot(4,1,4)

plot(timeVect,spaceSimilarity(indVect), ...
    'Color',lineColor(midIdx,:))

ylabel({'Similarity index', 'wrt t = 0'})
pos4 = get(gca,'Position');
set(gca,'XLim',xLim, ...
    'YAxisLocation',yAxisLocation, ...
    'YLim',[0 1], ...
    'PlotBoxAspectRatio',[ar(1)*pos1(4)/pos4(4) ar(2) ar(3)])
xlabel(xLabel)
box off

drawnow()
set(gcf,'Position',[190   141   413   912])

end

function interpColorMap = InterpColorMap(cRgb,xq)
%%
assert(all(xq >= 0 & xq <= 1), ...
    'Interpolated points must be between 0 and 1')

if size(xq,1) > 1
    xq = xq';
end

cHsv = rgb2hsv(cRgb);

hInterp = interp1(linspace(0,1,size(cHsv,1)), cHsv(:,1), xq);
sInterp = interp1(linspace(0,1,size(cHsv,1)), cHsv(:,2), xq);
vInterp = interp1(linspace(0,1,size(cHsv,1)), cHsv(:,3), xq);
cHsv = [hInterp' sInterp' vInterp'];

interpColorMap = hsv2rgb(cHsv);


end

function [s1, principleAngle] = DatasetSimilarity( input1, input2, varargin )
%% DatasetSimilarity A measure of the similarity between two datasets
%
% based on:
%
% Garcia C. A simple procedure for the comparison of covariance matrices.
% BMC Evol Biol. 2012;12: 222.
%
% Note: variance explained is normalized, and high numbers mean more
% similarity, unlike in the publication.

if ~all(size(input1,2) == size(input2,2))
    error('Input dimension mismatch')
end

[coef1,~,latent1] = pca(input1);
[coef2,~,latent2] = pca(input2);
scores21 = input2*coef1;
scores12 = input1*coef2;

var11 = latent1/sum(latent1);
var22 = latent2/sum(latent2);
var21 = (nanvar(scores21)/sum(nanvar(scores21)))';
var12 = (nanvar(scores12)/sum(nanvar(scores12)))';

%Total similarity
s1 = 1 - 1/4*sum((abs(var11 - var21) + abs(var22 - var12)));

principleAngle = acos(dot(coef1(:,1),coef2(:,1)));

end
