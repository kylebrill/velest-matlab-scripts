%SEMODEEVAL - run after single event mode has been run

printvar=0;  % If you want to save the output of vel_read_sta
printplot= 1;  % If you want plots saved of the statistics (.png files)

parentdir='/local/kabrill_res/velest_adjeqshpws/'; % Top VELEST Directory
outmoddir='data/'; % Directory of VELEST outputs
inmod='sim02rev.CNV'; %events being relocated in single event mode
runtitle='sem'; %title of single event mode run
outfiletag='.OUT'; % VELEST .OUT summary file from each run
stafiletag='.sta'; % station correction file from each run
xyzfile=[parentdir outmoddir 'fgtopo-2Kdatum.xyz']; % velest xyz file
datum=2000; % if datum change in velest
zshift=0; % actual velest Z shift from .cmn file
kmtolon=1.7968; %found in .OUT file of Single Event Mode Run
kmtolat=1.8440; %found in .OUT file of Single Event Mode Run

locfile='sem.loc'; % loc file from single event mode output
eventinfo=readseloc([parentdir outmoddir locfile]);

%account for deleted events during single event mode relocation
deleter=arrayfun(@(s) ~isempty(s.origin),eventinfo);
eventinfo=eventinfo(deleter);

% eventinfo(1:65).depth=-1000.*[eventinfo.depth]+zshift;

% add stations
[stationname,ptcor,stalat,stalon,staele] = ...
    vel_read_sta([parentdir outmoddir runtitle stafiletag],printvar);    
% get starting locations
[~,elatin,elonin,edepin,~,~,~]=readcnv([parentdir outmoddir inmod]);
edepin=-1000*edepin+datum+zshift;
elatin=elatin(deleter);
elonin=elonin(deleter);
edepin=edepin(deleter);

xe=cell(length(eventinfo),1);
ye=cell(length(eventinfo),1);
ze=cell(length(eventinfo),1);
% Event Error Ellipses
for in=1:length(eventinfo)
    [xe{in,1},ye{in,1},ze{in,1}]=ellipsoid(eventinfo(in,1).lon,...
        eventinfo(in,1).lat,...
        eventinfo(in,1).depth*-1000 + datum,...
        eventinfo(in,1).erx/kmtolon/60,...
        eventinfo(in,1).ery/kmtolat/60,...
        eventinfo(in,1).erz*1000);
end

%% Plots
% Basemap
EQP=figure('units','normalized','outerposition',[0 0 1 1]);

[xm,ym,zm]=xyzread(xyzfile);
mX=unique(xm); 
mY=unique(ym); mY=flipud(mY);

%mX=xm(1:length(find(xm==(xm(1)))));
%mY=ym(1:length(find(ym==(ym(1)))):end);
mZ=reshape(zm+datum,length(mX),length(mY))';

h(1)=surf(mX,mY,mZ,'FaceAlpha',0.5,'LineStyle','none');
hold on

% Plot Events with Error Ellipses cut by RMS limits
for in=1:length(eventinfo)
    if eventinfo(in,1).rms < 1 %&& eventinfo(in,1).gap < 180 && eventinfo(in,1).erz < 1.5
        h(2)=surf(xe{in,1},ye{in,1},ze{in,1},'FaceColor','r','FaceAlpha',(0.5-eventinfo(in,1).rms),'LineStyle',':');
    end
end
for in=1:length(eventinfo)
    if eventinfo(in,1).rms < 1 && eventinfo(in,1).gap < 180
        cp=plot3([eventinfo(:,1).lon],[eventinfo(:,1).lat],[eventinfo(:,1).depth]*-1000+datum,'.k');
    end
end

% Add contours
contour3(mX,mY,mZ,'k','LineWidth',2);

% Add stations and Corrections
stamkrc(1:length(ptcor),1)='b';
stamkrc(ptcor<0)='r';

for in = 1:length(stationname)
    plot3(stalon(in,1),stalat(in,1),staele(in,1)+datum,'^','MarkerSize',abs(ptcor(in,1))*4+8,...
        'MarkerFaceColor',stamkrc(in,1),'Color','k')
end

set(gca,'XLim',[-90.9329 -90.8287],'YLim',[14.4238 14.5280],'Zlim',[-5000,4000])% Fuego zoom in

axis square
%% Save plots
if printplot == 1;
    
    img2 = getframe(EQP);
    imwrite(img2.cdata, [parentdir runtitle 'eliplot.png']);
    view([90,0])
    img3 = getframe(EQP);
    imwrite(img3.cdata, [parentdir runtitle 'eliplotfromE.png']);
    view([0,0])
    img4 = getframe(EQP);
    imwrite(img4.cdata, [parentdir runtitle 'eliplotfromS.png']);
    view([0,90])
    img5 = getframe(EQP);
    imwrite(img5.cdata, [parentdir runtitle 'eliplotfromT.png']);
    
end

%% Copy .mod file for next step. 
% need a function which reads the velest.cmn file to figure out where the
% model is coming from. This function could also be used to streamline
% changes in the command line interactions with velest and allow slightly
% more rapid iterations and changes. 