%RVELMODEVALSUP Suppliment to rvelmodeval.m for plotting PWS events
% Three options for plotting extra results from rvelmodeval.m. First is
% focusd earthquake tracking for PWS events from best 1% of models, second
% is tracking PWS events for all random models, and third is a specific
% event chosen by number by the user.

%% Specific inputs
pws=5;    % How many PWS events in your .CNV, assume all at beginning
chosenevent=1;
printplot=1;

%% Focused EQ tracking - Best Velocity Mods PWS

EQF=figure('units','normalized','outerposition',[0 0 1 1]);
%setup plot
h(1)=surf(mX,mY,mZ,'FaceAlpha',0.4,'LineStyle','none');
hold on

% Add contours
contour3(mX,mY,mZ,'k','LineWidth',2);

% Add stations and mean corrections
avcor=mean(ptcorplotbest,2);
stamkrc(1:length(avcor),1)='r';
stamkrc(avcor<0)='b';
for in = 1:length(stationname)
    staplot=plot3(stlon(in,1),stlat(in,1),stele(in,1)+zshift,...
        '^','MarkerSize',abs(avcor(in,1))+8,...
        'MarkerFaceColor',stamkrc(in,1),'Color','k');
end

%plot event starting locations
[~,elat,elon,edep,~,~,~]=readcnv([parentdir outmoddir inmod]);
seqplot=plot3(elon(1:pws,1),elat(1:pws,1),-edep(1:pws,1)*1000+zshift,'g.','MarkerSize',20);

%plot event locations from best rms models
for eqn=length(best_mods):-1:1
[~,e1lat,e1lon,e1dep,~,~,~]=readcnv([parentdir velmoddir best_mods{eqn,1} filesep cnvfilename]);
plot3([elon(1:pws,1)';e1lon(1:pws,1)'],[elat(1:pws,1)';e1lat(1:pws,1)'],-[edep(1:pws,1)';e1dep(1:pws,1)']*1000+zshift,':k')
eeqplot=plot3(e1lon(1:pws,1),e1lat(1:pws,1),-e1dep(1:pws,1)*1000+zshift,'.','Color',[(length(best_mods)+1-eqn)/length(best_mods),0,0],'MarkerSize',(21-eqn*2));
end

[~,e2lat,e2lon,e2dep,~,~,~]=readcnv([parentdir velmoddir modfolder{iminsta,1} filesep cnvfilename]);
plot3([elon(1:pws,1)';e2lon(1:pws,1)'],[elat(1:pws,1)';e2lat(1:pws,1)'],-[edep(1:pws,1)';e2dep(1:pws,1)']*1000+zshift,':k')
beqplot=plot3(e2lon(1:pws,1),e2lat(1:pws,1),-e2dep(1:pws,1)*1000+zshift,'.','Color','b','MarkerSize',(20));

%set(gca,'XLim',[-90.9329 -90.8287],'YLim',[14.4238 14.5280])
set(gca,'XLim',[-90.9329 -90.8287],'YLim',[14.4238 14.5280],'Zlim',[-5000,4000])% Fuego zoom in
hold off
axis square
legend([staplot,seqplot,eeqplot,beqplot],{'Station (Red >= 0)',...
    ['                    ';'Starting Earthquakes';...
    '                    '],...
    ['Earthquakes from best RMS models';...
    'Red is best model by RMS        ';...
    'Black is last best model by RMS '],...
    'Blue is lowest RMS and Sta Corr '},'Location','SouthOutside')
title('PWS Event locations of best 1% random models by RMS')

%% Save Plots

if printplot == 1;
    
    img2 = getframe(EQF);
    imwrite(img2.cdata, [parentdir runtitle 'PWSbestplot.png']);
    view([90,0])
    img3 = getframe(EQF);
    imwrite(img3.cdata, [parentdir runtitle 'PWSbestplotfromE.png']);
    view([0,0])
    img4 = getframe(EQF);
    imwrite(img4.cdata, [parentdir runtitle 'PWSbestplotfromS.png']);
    view([0,90])
    img5 = getframe(EQF);
    imwrite(img5.cdata, [parentdir runtitle 'PWSbestplotfromT.png']);
    
end
%% Focused EQ tracking - All Velocity Mods PWS

EQA=figure('units','normalized','outerposition',[0 0 1 1]);

h(1)=surf(mX,mY,mZ,'FaceAlpha',0.4,'LineStyle','none');
hold on

% Add contours
contour3(mX,mY,mZ,'k','LineWidth',2);

% Add stations and mean corrections
avcor=mean(ptcorplotbest,2);
stamkrc(1:length(avcor),1)='r';
stamkrc(avcor<0)='b';
for in = 1:length(stationname)
    staplot=plot3(stlon(in,1),stlat(in,1),stele(in,1)+zshift,...
        '^','MarkerSize',abs(avcor(in,1))+8,...
        'MarkerFaceColor',stamkrc(in,1),'Color','k');
end

%plot event starting locations
[~,elat,elon,edep,~,~,~]=readcnv([parentdir outmoddir inmod]);
seqplot=plot3(elon(1:pws,1),elat(1:pws,1),-edep(1:pws,1)*1000+zshift,'g.','MarkerSize',20);

%plot event locations from all models
for eqn=length(modfolder):-1:1
[~,e1lat,e1lon,e1dep,~,~,e1rms,~,~]=readcnv([parentdir velmoddir modfolder{rms_rank(eqn,1),1} filesep cnvfilename]);
plot3([elon(1:pws,1)';e1lon(1:pws,1)'],[elat(1:pws,1)';e1lat(1:pws,1)'],-[edep(1:pws,1)';e1dep(1:pws,1)']*1000+zshift,':k')
eeqplot=plot3(e1lon(1:pws,1),e1lat(1:pws,1),-e1dep(1:pws,1)*1000+zshift,'.','Color',[(nummod+1-eqn)/nummod,0,0],'MarkerSize',(20-eqn*0.01));
end

%set(gca,'XLim',[-90.9329 -90.8287],'YLim',[14.4238 14.5280])
set(gca,'XLim',[-90.9329 -90.8287],'YLim',[14.4238 14.5280],'Zlim',[-5000,4000])% Fuego zoom in
hold off
axis square
legend([staplot,seqplot,eeqplot],{'Station (Red >= 0)',...
    ['                    ';'Starting Earthquakes';...
    '                    '],...
    'Earthquakes from all models'},'Location','SouthOutside')
title('PWS Event locations of all random models')

%% Save Plots

if printplot == 1;
    
    img2 = getframe(EQA);
    imwrite(img2.cdata, [parentdir runtitle 'PWSallplot.png']);
    view([90,0])
    img3 = getframe(EQA);
    imwrite(img3.cdata, [parentdir runtitle 'PWSallplotfromE.png']);
    view([0,0])
    img4 = getframe(EQA);
    imwrite(img4.cdata, [parentdir runtitle 'PWSallplotfromS.png']);
    view([0,90])
    img5 = getframe(EQA);
    imwrite(img5.cdata, [parentdir runtitle 'PWSallplotfromT.png']);
    
end


%% Focused EQ tracking - All Velocity Mods of Specific Event
for chosenevent=1:pws
EQS=figure('units','normalized','outerposition',[0 0 1 1]);

h(1)=surf(mX,mY,mZ,'FaceAlpha',0.4,'LineStyle','none');
hold on

% Add contours
contour3(mX,mY,mZ,'k','LineWidth',2);

% Add stations and mean corrections
avcor=mean(ptcorplotbest,2);
stamkrc(1:length(avcor),1)='r';
stamkrc(avcor<0)='b';

for in = 1:length(stationname)
    staplot=plot3(stlon(in,1),stlat(in,1),stele(in,1)+zshift,...
        '^','MarkerSize',abs(avcor(in,1))+8,...
        'MarkerFaceColor',stamkrc(in,1),'Color','k');
end

%plot event starting locations
[eorig,elat,elon,edep,~,~,~]=readcnv([parentdir outmoddir inmod]);
seqplot=plot3(elon(chosenevent,1),elat(chosenevent,1),-edep(chosenevent,1)*1000+zshift,'g.','MarkerSize',20);
eorig=eorig(chosenevent);

%plot event locations from all models
for eqn=length(modfolder):-1:1
[~,e1lat,e1lon,e1dep,~,~,~]=readcnv([parentdir velmoddir modfolder{rms_rank(eqn,1),1} filesep cnvfilename]);
plot3([elon(chosenevent,1)';e1lon(chosenevent,1)'],[elat(chosenevent,1)';e1lat(chosenevent,1)'],-[edep(chosenevent,1)';e1dep(chosenevent,1)']*1000+zshift,':k')
eeqplot=plot3(e1lon(chosenevent,1),e1lat(chosenevent,1),-e1dep(chosenevent,1)*1000+zshift,'.','Color',[(nummod+1-eqn)/nummod,0,0],'MarkerSize',(20-eqn*0.01));
end

set(gca,'XLim',[-90.9329 -90.8287],'YLim',[14.4238 14.5280])
% set(gca,'XLim',[min(xm),max(xm)],'YLim',[min(ym),max(ym)])
% set(gca,'XLim',[-90.9329 -90.8287],'YLim',[14.4238 14.5280],'Zlim',[-5000,4000])% Fuego zoom in
hold off
axis square
legend([staplot,seqplot,eeqplot],{'Station (Red >= 0)',...
    ['                      ';'Initial Event Location';...
    '                      '],...
    'Resulting location from all models'},'Location','SouthOutside')
titletext=['Event ' num2str(chosenevent) ': ' datestr(eorig,'yy-mm-dd HH:MM') ' - Locations from all random models'];
title(titletext)

%% Save Plots

if printplot == 1;
    
    img2 = getframe(EQS);
    imwrite(img2.cdata, [parentdir runtitle 'event' num2str(chosenevent,'%02d') 'plot.png']);
    view([90,0])
    img3 = getframe(EQS);
    imwrite(img3.cdata, [parentdir runtitle 'event' num2str(chosenevent,'%02d') 'plotfromE.png']);
    view([0,0])
    img4 = getframe(EQS);
    imwrite(img4.cdata, [parentdir runtitle 'event' num2str(chosenevent,'%02d') 'plotfromS.png']);
    view([0,90])
    img5 = getframe(EQS);
    imwrite(img5.cdata, [parentdir runtitle 'event' num2str(chosenevent,'%02d') 'plotfromT.png']);
    
end

%% Focused EQ tracking - Best Velocity Mods PWS

EQB=figure('units','normalized','outerposition',[0 0 1 1]);
%setup plot
h(1)=surf(mX,mY,mZ,'FaceAlpha',0.4,'LineStyle','none');
hold on

% Add contours
contour3(mX,mY,mZ,'k','LineWidth',2);

% Add stations and mean corrections
avcor=mean(ptcorplotbest,2);
stamkrc(1:length(avcor),1)='r';
stamkrc(avcor<0)='b';
for in = 1:length(stationname)
    staplot=plot3(stlon(in,1),stlat(in,1),stele(in,1)+zshift,...
        '^','MarkerSize',abs(avcor(in,1))+8,...
        'MarkerFaceColor',stamkrc(in,1),'Color','k');
end

%plot event starting locations
[eorig,elat,elon,edep,~,~,~]=readcnv([parentdir outmoddir inmod]);
seqplot=plot3(elon(chosenevent,1),elat(chosenevent,1),-edep(chosenevent,1)*1000+zshift,'g.','MarkerSize',20);
eorig=eorig(chosenevent);

%plot event locations from best rms models
for eqn=length(best_mods):-1:1
[~,e1lat,e1lon,e1dep,~,~,~]=readcnv([parentdir velmoddir best_mods{eqn,1} filesep cnvfilename]);
plot3([elon(chosenevent,1)';e1lon(chosenevent,1)'],[elat(chosenevent,1)';e1lat(chosenevent,1)'],-[edep(chosenevent,1)';e1dep(chosenevent,1)']*1000+zshift,':k')
eeqplot=plot3(e1lon(chosenevent,1),e1lat(chosenevent,1),-e1dep(chosenevent,1)*1000+zshift,'.','Color',[(length(best_mods)+1-eqn)/length(best_mods),0,0],'MarkerSize',(21-eqn*2));
text(e1lon(chosenevent,1),e1lat(chosenevent,1),-e1dep(chosenevent,1)*1000+zshift,num2str(eqn))
end

[~,e2lat,e2lon,e2dep,~,~,~]=readcnv([parentdir velmoddir modfolder{iminsta,1} filesep cnvfilename]);
plot3([elon(chosenevent,1)';e2lon(chosenevent,1)'],[elat(chosenevent,1)';e2lat(chosenevent,1)'],-[edep(chosenevent,1)';e2dep(chosenevent,1)']*1000+zshift,':k')
beqplot=plot3(e2lon(chosenevent,1),e2lat(chosenevent,1),-e2dep(chosenevent,1)*1000+zshift,'.','Color','b','MarkerSize',(20));
text(e2lon(chosenevent,1),e2lat(chosenevent,1),-e2dep(chosenevent,1)*1000+zshift,'minsta\rms')

%set(gca,'XLim',[-90.9329 -90.8287],'YLim',[14.4238 14.5280])
set(gca,'XLim',[-90.9329 -90.8287],'YLim',[14.4238 14.5280],'Zlim',[-5000,4000])% Fuego zoom in
hold off
axis square
legend([staplot,seqplot,eeqplot],{'Station (Red >= 0)',...
    ['                   ';'Starting Earthquake';...
    '                   '],...
    ['Earthquake from best RMS models';...
    'Red is best model by RMS       ';...
    'Black is last best model by RMS']},'Location','SouthOutside')
titletext=['Event ' num2str(chosenevent) ': ' datestr(eorig,'yy-mm-dd HH:MM') ' - Locations from best 1% random models by RMS'];
title(titletext)

%% Save Plots

if printplot == 1;
    
    img2 = getframe(EQB);
    imwrite(img2.cdata, [parentdir runtitle 'event' num2str(chosenevent,'%02d') 'bestplot.png']);
    view([90,0])
    img3 = getframe(EQB);
    imwrite(img3.cdata, [parentdir runtitle 'event' num2str(chosenevent,'%02d') 'bestplotfromE.png']);
    view([0,0])
    img4 = getframe(EQB);
    imwrite(img4.cdata, [parentdir runtitle 'event' num2str(chosenevent,'%02d') 'bestplotfromS.png']);
    view([0,90])
    img5 = getframe(EQB);
    imwrite(img5.cdata, [parentdir runtitle 'event' num2str(chosenevent,'%02d') 'bestplotfromT.png']);
    
end
end