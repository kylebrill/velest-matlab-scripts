%simuleval - Use to Evaluate individual velest models from preliminary runs
clear
%% Initialize variables
printvar = 0;  % If you want locations_raw.dat, sta.dat, and gap.dat files
printplot= 1;  % If you want plots saved of the statistics (.png files)

%% Direct script to files important files and AWK Scripts
parentdir='/local/kabrill_res/velest_adjeqshpws/'; % Top VELEST Directory
outmoddir='data/'; % Directory of VELEST outputs
inmod='sim02rev.CNV';
runtitle='sim03';
shotfile='sim02sh.CNV'; % File with shots - if using
outfiletag='.OUT'; % VELEST .OUT summary file from each run
stafiletag='.sta'; % station correction file from each run
cnvfiletag='.CNV'; % earthquake data from each run
modfiletag='.mod'; %'.NW1_.mod'
xyzfile=[parentdir outmoddir 'fgtopo-2Kdatum.xyz']; % velest xyz file
vlat=14.474471; % actual volcano coordinates
vlon=-90.880940;
datumshift=2000; % if used in velest topofile
zshift=0; %if used in velest (negative meters)

% awk script locations
scriptfolder = ['/mnt/data/MATLAB''' ' Drive''/volcanoseismo/KABrillCode/velest-matlab-scripts/'];
vz_RMS_mod = [scriptfolder 'vz_RMS_mod.awk'];
vz_output = [scriptfolder 'vz_output.awk'];
vz_input = [scriptfolder 'vz_input.awk'];
% vz_EQ = [scriptfolder 'vz_EQ.awk'];
% ex_sim = [scriptfolder 'extract_simultaneousloc.awk'];

%% Access .OUT file and import variables to MATLAB using AWK scripts

% Get RMS from Models - AWK
filestring = [parentdir outmoddir runtitle outfiletag];
commandstringrms = ['awk -f ' vz_RMS_mod ' ' filestring];
[~,rmsval]=system(commandstringrms);
rmsval=regexprep(rmsval,'\r\n|\n|\r','');%get rid of hidden characters

% Get Velocity Models - AWK
commandstringoutvelmod = ['awk -f ' vz_output ' ' filestring ];
[~,outvelmod]=system(commandstringoutvelmod);
outvelmod=str2num(outvelmod); %#ok<ST2NM>
commandstringinvelmod = ['awk -f ' vz_input ' ' filestring ];
[~,invelmod]=system(commandstringinvelmod);
invelmod=str2num(invelmod); %#ok,<ST2NM>

% Get Station Corrections
[stationname,ptcor,stalat,stalon,staele] = ...
    vel_read_sta([parentdir outmoddir runtitle stafiletag],printvar);

% Get Initial Eq locations, Event GAPs and Pick Weights
[~,eqlatin,eqlonin,eqdepin,~,gapin,pweightsin]=...
    readcnv([parentdir outmoddir inmod],printvar);
eqdepin=-1000*eqdepin+datumshift+zshift;

% Get Final Eq Locations, Event GAPs and Pick Weights
[~,eqlatout,eqlonout,eqdepout,~,gapout,pweightsout]=...
    readcnv([parentdir outmoddir runtitle cnvfiletag],printvar);
eqdepout=-1000*eqdepout+datumshift;

% Check to see if final eq locations and initial eq locations are equal
if length(eqlatin)~=length(eqlatout)
    % Check to see if adding shot data will fix the problem
    if exist([parentdir outmoddir shotfile],'file') == 2
        [~,shlatin,shlonin,shdepin,~,shgapin,shweightsin]=...
            readcnv([parentdir outmoddir shotfile]); % load shots
        eqlatin=[eqlatin;shlatin];%append shot data to eq locations
        eqlonin=[eqlonin;shlonin];
        eqdepin=[eqdepin;-1000*shdepin+datumshift];
        gapin=[gapin;shgapin];
        pweightsin=pweightsin+shweightsin;
    else
        error('The shotfile cannot be found.')
    end
    if length(eqlatin)~=length(eqlatout)
        error('Earthquakes in do not match Earthquakes out. Something is wrong')
    end
end
    

%% Plots 
%% Set up figure
CHK=figure('Units','Inches','Position',[0, 0, 8.5,11],...
    'PaperUnits', 'Inches', 'PaperSize', [8.5, 11]);

%% Set up basemap data 
[xm,ym,zm]=xyzread(xyzfile);
mX=unique(xm);
mY=unique(ym); mY=flipud(mY);
mZ=reshape(zm+datumshift,length(mX),length(mY))';

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
set(cp,'Unit','Inches','Position',[0.75 3.75 5 5])
set(cp,'yaxislocation','right')
set(cp,'YTickLabelRotation',90)
axis equal
% Interactive option available here to get a good set of limits

%% organize eq data for better looking plots

% outliers will be used to keep track of eqlocations that have been moved
% to the edge of the plot for a different marker
outliersout=logical((eqlatout > max (mY))+(eqlatout < min (mY))+...
    (eqlonout > max (mX))+(eqlonout < min(mX)));
outliersin=logical((eqlatin > max (mY))+(eqlatin < min (mY))+...
    (eqlonin > max (mX))+(eqlonin < min(mX)));

% adjust eq locations to keep them on the map
eqlatout(eqlatout>max(mY))=max(mY);
eqlatout(eqlatout<min(mY))=min(mY);
eqlonout(eqlonout>max(mX))=max(mX);
eqlonout(eqlonout<min(mX))=min(mX);
eqlatin(eqlatin>max(mY))=max(mY);
eqlatin(eqlatin<min(mY))=min(mY);
eqlonin(eqlonin>max(mX))=max(mX);
eqlonin(eqlonin<min(mX))=min(mX);

% different markers for outliers
markersout(1:length(outliersout),1)='.';
markersout(outliersout)='s';
markersin(1:length(outliersin),1)='.';
markersin(outliersin)='s';

%% Plot mapview eqs based on gap color
for in=1:length(gapin)
    plot(eqlonin(in,1),eqlatin(in,1),markersin(in),'Color',[0,gapin(in,1)/360,0]);
end
plot([eqlonout';eqlonin'],[eqlatout';eqlatin'],':k')
for in=1:length(gapout)
    plot(eqlonout(in,1),eqlatout(in,1),markersout(in),'Color',[gapout(in,1)/360,0,0])
end
for in=1:5
    plot(eqlonout(in,1),eqlatout(in,1),'sk','MarkerSize',in*3)
end

%% plot stations with corrections
stamkrc(1:length(ptcor),1)='b';
stamkrc(ptcor<0)='r';

for in = 1:length(stationname)
    plot(stalon(in,1),stalat(in,1),'^','MarkerSize',abs(ptcor(in,1))*4+8,...
        'MarkerFaceColor',stamkrc(in,1),'Color','k')
end

set(cp,'XLim',[-90.9329 -90.8287],'YLim',[14.4238 14.5280])% Fuego zoom in
% set(cp,'XLim',[min(mX) max(mX)],'YLim',[min(mY) max(mY)])% Reset zoom
hold off

%% plot N-S profile from W
nsp=axes('Unit','Inches','Position',[6 3.75 1.5 5]);
plot(NSz,mY,'k')
set(nsp,'xdir','reverse')
xlim([-5000,4000])
ylim([min(mY) max(mY)])
set(nsp,'YTickLabel','')
set(nsp,'xaxislocation','top')
set(nsp,'XTickLabelRotation',90)
hold on

% plot xsection eqs based on gap color
for in=1:length(gapin)
    plot(eqdepin(in,1),eqlatin(in,1),markersin(in),'Color',[0,gapin(in,1)/360,0]);
end
plot([eqdepout';eqdepin'],[eqlatout';eqlatin'],':k')
for in=1:length(gapout)
    plot(eqdepout(in,1),eqlatout(in,1),markersout(in),'Color',[gapout(in,1)/360,0,0])
end
for in=1:5
    plot(eqdepout(in,1),eqlatout(in,1),'sk')
end
hold off

%% plot E-W profile from S
ewp=axes('Unit','Inches','Position',[0.75 2 5 1.5]);
plot(mX,EWz,'k')
xlim([min(mX) max(mX)])
ylim([-5000,4000])
set(ewp,'XTickLabel','')
hold on
% Plot xsection eqs based on gap color
for in=1:length(gapin)
    plot(eqlonin(in,1),eqdepin(in,1),markersin(in),'Color',[0,gapin(in,1)/360,0]);
end
plot([eqlonout';eqlonin'],[eqdepout';eqdepin'],':k')
for in=1:length(gapout)
    plot(eqlonout(in,1),eqdepout(in,1),markersout(in),'Color',[gapout(in,1)/360,0,0])
end
for in=1:5
    plot(eqlonout(in,1),eqdepout(in,1),'sk')
end
hold off

%% link axes for zooming
linkaxes([cp,ewp],'x')
linkaxes([cp,nsp],'y')

%% plot velocity model
vmp=axes('Unit','Inches','Position',[6 .5 1.5 3]);
plot(outvelmod(:,2),-1000*outvelmod(:,1)+datumshift,'k')
hold on
plot(invelmod(:,2),-1000*invelmod(:,1)+datumshift,'k--')
set(vmp,'yaxislocation','right')
vmpleg=legend('Out Model','Starting Model','Location','best');
set(vmpleg,'FontSize',6);
ylim([-10000 4000])
line([0.34 0.34],[-10000 4000],'Color','r')

%% plot old and new gaps
vmp=axes('Unit','Inches','Position',[4.5 .5 1.25 1.25]);
edges=0:15:360;
histogram(gapout,edges,'facecolor','r','facealpha',.5)
hold on
histogram(gapin,edges,'facecolor','g','facealpha',.5)
xlim([0,360])
set(vmp,'XTick',[90,180,270])

%% print RMS
txt=uicontrol('style','text','Unit','Inches','Position',[.75 .5 3 1.25],...
    'BackgroundColor','white',...
    'String',{['Run Title: ' runtitle];['Starting Model: ' inmod];...
    ['RMS of Out Model = ' rmsval];[num2str(length(gapin)) ' Events'];...
    'Histogram is GAP of Events'; 'Velocities in km/s, Depths in m';
    'Red is out model, green is starting model'});

%% Save as .pdf
if printplot
    s=hgexport('readstyle','default');
    s.Format='pdf';
    hgexport(CHK,[parentdir outmoddir runtitle '.pdf'],s)
end

%% Print new mod file

modstr = [parentdir outmoddir runtitle modfiletag];
if exist(modstr,'file') ~= 2
    inmodstr=[parentdir outmoddir strrep(inmod,'.CNV','') modfiletag];
    if strcmp(inmodstr(end-4),'v') % allow for inputs produced by cnvcutter
        inmodstr=strrep(inmodstr,'rev','');
    end
    
    fidin=fopen(inmodstr);
    if fidin == -1
        error('Can''t find original .mod file to copy. Please do this manually.')
    end
    tmpline=fgetl(fidin);
    fidout=fopen(modstr,'w');
    linnum=1;
    i=1:2:length(outvelmod);
    while ischar(tmpline)
        if ~isempty(tmpline)
            if linnum <= 2
                fprintf(fidout,[tmpline '\n']);
                linnum = linnum+1;
            elseif linnum == 3
                vdamp=str2double(tmpline(21:26));
                fprintf(fidout,' %4.2f      %5.1f    %06.2f            %s\n',outvelmod(1,2),outvelmod(1,1),vdamp,'P-VELOCITY MODEL');
                linnum = linnum+1;
            else
                vdamp=str2double(tmpline(21:26));
                fprintf(fidout,' %4.2f      %5.1f    %06.2f\n',outvelmod(i(linnum-2),2),outvelmod(i(linnum-2),1),vdamp);
                linnum = linnum+1;
            end
        end
        tmpline=fgetl(fidin);
    end
    fclose(fidin);
    fclose(fidout);
end
