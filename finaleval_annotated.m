%finaleval - Plot final velest output

%% Initialize variables
printplot= 1;  % If you want plots saved of the statistics (.pdf file)
inc_error= 0;
%% Direct script to files important files and AWK Scripts
% parentdir='/local/kabrill_res/velest_adjeqshpws/'; % Top VELEST Directory
parentdir='D:\Google Drive\data-archive\velest_archives\velest_adjeqshpws\';
outmoddir='data/'; % Directory of VELEST outputs
% runtitle='sem';
% shotfile='sim02sh.CNV'; % File with shots - if using
runtitle='sim03';
shotfile='sim03sh.CNV';
outfiletag='.OUT'; % VELEST .OUT summary file from each run
stafiletag='.sta'; % station correction file from each run
cnvfiletag='.CNV'; % earthquake data from each run
modfiletag='.mod'; %'.NW1_.mod'
xyzfile=[parentdir outmoddir 'fgtopo-2Kdatum.xyz']; % velest xyz file
vlat=14.474471; % actual volcano coordinates
vlon=-90.880940;
datumshift=2; % if used in velest topofile
zshift=0; %if used in velest (negative meters)

% awk script locations
if isunix
    scriptfolder = ['/mnt/data/MATLAB''' ' Drive''/volcanoseismo/KABrillCode/velest-matlab-scripts/'];
    vz_RMS_mod = [scriptfolder 'vz_RMS_mod.awk'];
    vz_output = [scriptfolder 'vz_output.awk'];
    vz_input = [scriptfolder 'vz_input.awk'];
    filestring = [parentdir outmoddir runtitle outfiletag];
elseif ispc
    scriptfolder = '"D:\MATLAB Drive\volcanoseismo\KABrillCode\velest-matlab-scripts\';
    vz_RMS_mod = [scriptfolder 'vz_RMS_mod.awk"'];
    vz_output = [scriptfolder 'vz_output.awk"'];
    vz_input = [scriptfolder 'vz_input.awk"'];
    filestring = ['"' parentdir outmoddir runtitle outfiletag '"'];
end


%% Access .OUT file and import variables to MATLAB using AWK scripts

% Get RMS from Models - AWK
commandstringrms = ['awk -f ' vz_RMS_mod ' ' filestring];
[~,rmsval]=system(commandstringrms);
rmsval=regexprep(rmsval,'\r\n|\n|\r','');%get rid of hidden characters

% Get Velocity Models - AWK
% commandstringoutvelmod = ['awk -f ' vz_output ' ' filestring ];
% [~,outvelmod]=system(commandstringoutvelmod);
% outvelmod=str2num(outvelmod); %#ok<ST2NM>
commandstringinvelmod = ['awk -f ' vz_input ' ' filestring ];
[~,outvelmod]=system(commandstringinvelmod);
outvelmod=str2num(outvelmod); %#ok,<ST2NM>

% Get Station Corrections
[stationname,ptcor,stalat,stalon,staele] = ...
    vel_read_sta([parentdir outmoddir runtitle stafiletag],0);

% Get Final Eq Locations, Event GAPs and Pick Weights
[~,eqlatout,eqlonout,eqdepout,~,gapout,rmsout,pweightsout,picksout]=...
    readcnv([parentdir outmoddir runtitle cnvfiletag],0);
eqdepout=-1*eqdepout+datumshift;

% Load shot data if included
if exist('shotfile','var')
    if exist([parentdir outmoddir shotfile],'file') == 2
        [~,shlatin,shlonin,shdepin,~,shgapin,shweightsin]=...
            readcnv([parentdir outmoddir shotfile]); % load shots
        shdepin=-1*shdepin+datumshift;
    end
end

% Load sigmavalues
if inc_error
    locfile=['D:\Documents\Visual Studio 2012\Projects\velest\velest\' ...
        'Release\velestdata\sem2.loc'];
    eventinfo = readseloc(locfile);
end
%% Plots 
%% Set up figure
CHK=figure('Units','Inches','Position',[7, 2, 6,7.75],...
    'PaperUnits', 'Inches', 'PaperSize', [6, 7.75],...
    'DefaultAxesFontSize',10,'DefaultAxesFontName','Arial');

%% Set up basemap data 
[xm,ym,zm]=xyzread(xyzfile);
mX=unique(xm);
mY=unique(ym); mY=flipud(mY);
mZ=reshape(zm./1000+datumshift,length(mX),length(mY))';

% get center of volcano
seeklat=abs(mY-vlat);
[~,inlat]=min(seeklat);
mvlat=mY(inlat);
seeklon=abs(mX-vlon);
[~,inlon]=min(seeklon);
mvlon=mX(inlon);

% get N-S profile from W
NSz=mZ(:,inlon);

% get E-W profile from S
EWz=mZ(inlat,:);

% plot contour map view
contour(mX,mY,mZ,'k')
hold on
cp=gca;
set(cp,'Unit','Inches','Position',[0.25 3.5 4 4])
set(cp,'xaxislocation','top')
set(cp,'YTickLabelRotation',90)
set(cp,'FontName','Arial')
axis equal
% Interactive option available here to get a good set of limits

%% organize eq data for better looking plots

% outliers will be used to keep track of eqlocations that have been moved
% to the edge of the plot for a different marker
outliersout=logical((eqlatout > max (mY))+(eqlatout < min (mY))+...
    (eqlonout > max (mX))+(eqlonout < min(mX)));

% adjust eq locations to keep them on the map
eqlatout(eqlatout>max(mY))=max(mY);
eqlatout(eqlatout<min(mY))=min(mY);
eqlonout(eqlonout>max(mX))=max(mX);
eqlonout(eqlonout<min(mX))=min(mX);

% different markers for outliers
markersout(1:length(outliersout),1)='.';
markersout(outliersout)='s';

%% Plot mapview eqs based on gap color
for in=1:length(gapout)
    plot(eqlonout(in,1),eqlatout(in,1),markersout(in),'Color',[0,0,0],'MarkerSize',15)%[gapout(in,1)/360,0,0])
end
for in=1:5
    plot(eqlonout(in,1),eqlatout(in,1),'sk','MarkerSize',10)%,'MarkerSize',in*3)
    clust=in;
%     if in>2
%         clust=in+1;
%     end
    if in==1 || in==4
    text(eqlonout(in,1),eqlatout(in,1),[num2str(clust) '\rightarrow   ' ],...
        'FontName','Arial','HorizontalAlignment','right')
    else
    text(eqlonout(in,1),eqlatout(in,1),['  \leftarrow' num2str(clust)],...
        'FontName','Arial')
    end
end
if exist('shotfile','var')
    plot(shlonin,shlatin,'*k')
end
if inc_error
    loncenter=eqlonout(5,1);
    latcenter=eqlatout(5,1);
    zcenter=eqdepout(5,1);
    thta=0:0.01:2*pi;
    elons=0.1*km2deg(eventinfo(5).erx)*cos(thta)+loncenter;
    elats=0.1*km2deg(eventinfo(5).ery)*sin(thta)+latcenter;
    plot(elons,elats,'r--')
end
%% plot stations with corrections
stamkrc(1:length(ptcor),1)={[0.5 0.5 0.5]};
stamkrc(ptcor<0)={[0.85 0.85 0.85]};

for in = 1:length(stationname)
    plot(stalon(in,1),stalat(in,1),'^','MarkerSize',abs(ptcor(in,1))*4+8,...
        'MarkerFaceColor',stamkrc{in,1},'Color','k','LineWidth',1.25)
end

% set(cp,'XLim',[-90.9329 -90.8287],'YLim',[14.4238 14.5280])% Fuego zoom in
% set(cp,'XLim',[-90.9417  -90.8245],'YLim',[14.4329   14.5501])%Acatenango Included
set(cp,'XLim',[-90.9034 -90.8570],'YLim',[14.4496 14.4960])%only events
% set(cp,'XLim',[-90.915 -90.845],'YLim',[14.44 14.51])
% set(cp,'XLim',[min(mX) max(mX)],'YLim',[min(mY) max(mY)])% Reset zoom
%%
lonlab=get(cp,'XTickLabel');
% lontag=repmat({[char(0176) 'W']},size(lonlab));
% lonlab=strcat(lonlab,lontag);
lonlab(2:2:8)={[]};
lonlab(9)={[]};
set(cp,'XTickLabel',lonlab)
latlab=get(cp,'YTickLabel');
% lattag=repmat({[char(0176) 'N']},size(latlab));
% latlab=strcat(latlab,lattag);
latlab(1)={[]};
latlab(2:2:10)={[]};
set(cp,'YTickLabel',latlab)
%%
sb1=line([-90.875 (-90.875+km2deg(1))],[14.453 14.453],'LineWidth',6,'Color','k');
sb2=line([-90.87501 (-90.875+km2deg(0.5))],[14.453 14.453],'LineWidth',4,'Color','w');
text((-90.875+km2deg(0.9)),14.455,'1 km')
text((-90.875+km2deg(0.4)),14.455,'0.5')
hold off

%% plot N-S profile from W
nsp=axes('Unit','Inches','Position',[4.25 3.5 1.5 4]);
plot(NSz,mY,'k')
set(nsp,'xdir','reverse')
xlim([1,4])
ylim([min(mY) max(mY)])
set(nsp,'YTickLabel','')
set(nsp,'XTickLabel','')
% set(nsp,'xaxislocation','top')
% set(nsp,'XTickLabelRotation',90)
set(nsp,'FontName','Arial')
hold on

% plot xsection eqs
for in=1:length(gapout)
    plot(eqdepout(in,1),eqlatout(in,1),markersout(in),'Color',[0,0,0],'MarkerSize',15)%[gapout(in,1)/360,0,0])
end
for in=1:5
    plot(eqdepout(in,1),eqlatout(in,1),'sk','MarkerSize',10)
    clust=in;
%     if in>2
%         clust=in+1;
%     end
    if in == 1 
        eqlatout(in,1)=eqlatout(in,1)+0.0005;
    elseif in == 2
        eqlatout(in,1)=eqlatout(in,1)-0.0005;
    end
    text(eqdepout(in,1),eqlatout(in,1),['  \leftarrow' num2str(clust)],...
        'FontName','Arial')
end
if exist('shotfile','var')
    plot(shdepin,shlatin,'*k')
end
if inc_error
    depcenter=eqdepout(5,1);
    edepslat=0.1*eventinfo(5).erz*cos(thta)+depcenter;
    elatsdep=0.1*km2deg(eventinfo(5).ery)*sin(thta)+latcenter;
    plot(edepslat,elatsdep,'r--')
end
hold off

%% plot E-W profile from S
ewp=axes('Unit','Inches','Position',[0.25 2 4 1.5]);
plot(mX,EWz,'k')
xlim([min(mX) max(mX)])
ylim([1,4])
set(ewp,'XTickLabel','')
set(ewp,'YTickLabel','')
set(ewp,'FontName','Arial')
hold on
% Plot xsection eqs
for in=1:length(gapout)
    plot(eqlonout(in,1),eqdepout(in,1),markersout(in),'Color',[0,0,0],'MarkerSize',15)%[gapout(in,1)/360,0,0])
end
for in=1:5
    plot(eqlonout(in,1),eqdepout(in,1),'sk','MarkerSize',10)
    clust=in;
%     if in>2
%         clust=in+1;
%     end
    if in==1
        text(eqlonout(in,1)-0.0003,eqdepout(in,1)-0.07,{'\uparrow';num2str(clust)},...
            'FontName','Arial','VerticalAlignment','top','HorizontalAlignment','center')
    elseif in==2 || in == 3 
        text(eqlonout(in,1),eqdepout(in,1)-0.07,{'\uparrow';num2str(clust)},...
            'FontName','Arial','VerticalAlignment','top','HorizontalAlignment','center')
    elseif in==4
        text(eqlonout(in,1),eqdepout(in,1),[num2str(clust) '\rightarrow  '],...
            'FontName','Arial','HorizontalAlignment','right')
    else
        text(eqlonout(in,1),eqdepout(in,1),['  \leftarrow' num2str(clust)],...
            'FontName','Arial')
    end
end
if exist('shotfile','var')
    plot(shlonin,shdepin,'*k')
end
if inc_error
    elonsdep=0.1*km2deg(eventinfo(5).erx)*cos(thta)+loncenter;
    edepslon=0.1*eventinfo(5).erz*sin(thta)+depcenter;
    plot(elonsdep,edepslon,'r--')
end
hold off

%% link axes for zooming
linkaxes([cp,ewp],'x')
linkaxes([cp,nsp],'y')
ewp.YTickLabelRotation=90;
ewp.YTickLabel={'','2km','3km',''};
% ewp.YAxisLocation='right';
nsp.XTickLabel={'','2km','3km',''};
nsp.XAxisLocation='top';
%% plot velocity model
vmp=axes('Unit','Inches','Position',[4.325 .5 1.25 2.75]);
plot(outvelmod(:,2),-1*outvelmod(:,1)+datumshift,'k')
hold on
% plot(invelmod(:,2),-1000*invelmod(:,1)+datumshift,'k--')
set(vmp,'yaxislocation','right')
% set(vmp,'YTickLabelRotation',90)
set(vmp,'FontName','Arial')
% vmpleg=legend('Out Model','Starting Model','Location','best');
% set(vmpleg,'FontSize',6);
% line([0.34 0.34],[-10000 4000],'Color','r')

for in=1:2:sum(outvelmod(:,1)<2)
    xpos=outvelmod(in,2)+0.2;
    if xpos>(ceil(max(outvelmod(:,2))))/2
        xpos=0.5;
    end
    ypos=-1*(outvelmod(in,1)+outvelmod(in+1))/2+datumshift;
    if ypos<-10
        ypos=-9;
    end
    text(xpos,ypos,[num2str(outvelmod(in,2),3) ' km/s'],'FontName','Arial')
end
% text(1.5,-9,[num2str(outvelmod(13,2),3) ' km/s'],'FontName','Arial')
% ylim([-14000 4000])
ylim([-2 4]) % These will be changed based on what looks good
xlim([0 ceil(max(outvelmod(:,2)))])
title('    P-wave Velocity')
ylabel('Elevation (km)','Position',[6.0611   -0.0208         0])
xlabel('Velocity (km/s)')
%% plot RMS
vmp=axes('Unit','Inches','Position',[.5 .5 1.25 1.25]);
edges=0:0.025:0.4;
histogram(rmsout,edges,'facecolor','k','facealpha',.5)
% set(vmp,'yaxislocation','right')
set(vmp,'XTick',[0.1,0.2,0.3])
set(vmp,'FontName','Arial')
ylabel('Event Count');%,'Rotation',270,'VerticalAlignment','Bottom','Position',[0.4444    4.0000         0])
xlabel('RMS error')
%% print RMS
% txt=uicontrol('style','text','Unit','Inches','Position',[.75 .5 3 1.25],...
%     'BackgroundColor','white','FontName','Arial',...
%     'String',{'Velocity Model for Fuego Volcano';...
%     ['Final RMS Error of Model = ' num2str(rmsval,4)];...
%     [num2str(length(eqdepout)) ' Events   ' num2str(length(shdepin)) ' Shots'];...
%     'Histogram is RMS of Each Event';...
%     'Velocities in km/s, Depths in m';
%     'Points are Earthquakes';'Asterisks are Shots';'Squares are PWS events'});

%% plot station corrections as boxplot
stationnameplot=repmat(cell2mat(stationname),[size(ptcor,2),1]);
ptcorplot=reshape(ptcor,[size(ptcor,1)*size(ptcor,2),1]);
stc=axes('Unit','Inches','Position',[2.25 0.5 2 1.25]);
boxplot(ptcorplot,stationnameplot,'Color',[0,0,0]);
% xlabel('Station Names')
set(gca,'XTickLabelRotation',45);
ylabel('Seconds','Position',[-0.7674    0.0550         0]);
title('Station Corrections')
%% Save as .pdf
savename=strrep(runtitle,'sim03','finalshallow-5PWSCEAnnotated');
if printplot
    s=hgexport('readstyle','default');
    s.Format='pdf';
    hgexport(CHK,[parentdir savename '.pdf'],s)
end