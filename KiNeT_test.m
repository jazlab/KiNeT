dt = 0.01;
nBin = 5;
midBin = median(1:nBin);
responseDur = linspace(0.5,1,nBin);
slope = 1./responseDur;

%% Trajectories in 3D, multiple gains
figure('Position',[462 356 1079 556])
ax1 = subplot(1,2,2);
hold on;
ax2 = subplot(3,2,1);
hold on;
ax3 = subplot(3,2,3);
hold on;
ax4 = subplot(3,2,5);
hold on;

rotation = pi/2;

scores = nan(3,nBin,ceil(max(responseDur)/dt));

for j = 1:2
    if j == 1
        lineColor = InterpColorMap( ...
            [0.75 0.75 0.75; 0.25 0.25 0.25],linspace(0,1,numel(responseDur)));
    else
        cHsv150 = rgb2hsv([1 0 0]);
        cHsv150 = [cHsv150(1) 0.25 1; cHsv150(1) 1 0.5];
        cRgb150 = hsv2rgb(cHsv150);
        lineColor = InterpColorMap(cRgb150,linspace(0,1,numel(responseDur)));
    end
    
    for i = numel(responseDur):-1:1
        t = dt:dt:responseDur(i);
        nt = numel(t);
        scores(1,i,1:nt) = 3*slope(i)*t;
        scores(2,i,1:nt) = sin(pi*t/responseDur(i))/responseDur(i) + wrev(slope(i)^2*t);
        scores(3,i,1:nt) = 0.5*j*ones(size(t)) + wrev(j*slope(i)*t);
        
        for k = 1:numel(t)
            Theta = rotation*t(k)/t(end);
            rotationMat = [cos(Theta) -1*sin(Theta); ...
                sin(Theta) cos(Theta)];
            temp = [scores(2,i,k) scores(3,i,k)]*rotationMat;
            scores(2,i,k) = temp(1);
            scores(3,i,k) = temp(2);
        end
            
        
        plot3(ax1,squeeze(scores(2,i,:)),squeeze(scores(1,i,:)), ...
            squeeze(scores(3,i,:)),'Color',lineColor(i,:));
        plot(ax3,t,squeeze(scores(1,i,1:nt)),'Color',lineColor(i,:));
        plot(ax2,t,squeeze(scores(2,i,1:nt)),'Color',lineColor(i,:));
        plot(ax4,t,squeeze(scores(3,i,1:nt)),'Color',lineColor(i,:));
        
        hTemp = plot3(ax1,scores(2,i,1),scores(1,i,1),scores(3,i,1), ...
            'Marker','o', ...
            'Color',lineColor(i,:), ...
            'LineStyle','none');
        if j == 1 && i == numel(responseDur)
            hLeg(1) = hTemp;
        end
        
        hTemp = plot3(ax1,scores(2,i,nt),scores(1,i,nt),scores(3,i,nt), ...
            'Marker','x', ...
            'Color',lineColor(i,:), ...
            'LineStyle','none');
        if j == 1 && i == numel(responseDur)
            hLeg(2) = hTemp;
        end
        
    end
end
set(ax1,'DataAspectRatio',[1 1 1], ...
    'CameraPosition',[-1.5189  -13.9889   36.6751])
grid(ax1)

legend(hLeg,'Set','Go', ...
    'Location','Northwest')

xlabel(ax1,'Rate, Unit 1')
ylabel(ax1,'Rate, Unit 2')
zlabel(ax1,'Rate, Unit 3')
xlabel(ax3,'Time')
ylabel(ax2,{'Rate,','Unit 1'})
ylabel(ax3,{'Rate,','Unit 2'})
set(ax1,'XTick',[],'YTick',[],'ZTickLabel',[]);
set(ax2,'XTick',[],'YTick',[]);
set(ax3,'XTick',[],'YTick',[]);

%%
KiNeT(scores,dt)
