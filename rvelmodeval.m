%RVELMODEVAL - Evaluate results from rvelmod.m
% Analyzes the outputs from rvelmod.m which produces random velocity models
% and then runs the VELEST program on each model to simultaneously invert
% for velocity and earthquake locations based on picks made in ANTELOPE
% using DBLOC2. 
%
% The script then plots a histogram of RMS of all runs, a histogram of the
% pick weights comming out of ANTELOPE, all the velocity models shaded by
% their RMS (lower RMS values showing up as brighter red), the velocity
% models with the lowest 1% of RMS values with mean, median, and mode of
% those lowest 1% RMS models, box and whisker plots for station corrections
% from each run, box and whisker plots for station corrections for the runs
% with the lowest 1% of RMS values, GAP values from all events from all
% modesl and finally GAP values from all events from the models with the
% lowest 1% of RMS values.
%
% Written by Kyle Brill 
% Contributions by Fede Lanza
% Last edited by Kyle Brill on 2016/06/01 at 22:00 EDT

%% Initialize variables
printvar = 0;  % If you want locations_raw.dat, sta.dat, and gap.dat files
onefig = 1; % Print stats as one fig or plot up figures in separate windows
printplot= 1;  % If you want plots saved of the statistics (.png files)

%% Direct script to files important files and AWK Scripts
parentdir='/local/kabrill_res/velest_adjeqshpws/'; % Top VELEST Directory
velmoddir='vmsSW2iadj/'; % Directory of VELEST random models
outmoddir='data/'; % Directory of VELEST outputs
runtitle='SW2irandadj';
inmod='irandadj.CNV';
outfilename=[runtitle '.OUT']; % VELEST .OUT summary file from each run
stafilename=[runtitle '.sta']; % station correction file from each run
cnvfilename=[runtitle '.CNV']; % earthquake data from each run
xyzfile=[parentdir outmoddir 'fgtopo-2Kdatum.xyz']; % velest xyz file
zshift=2000; % if datum change in velest

% awk script locations
scriptfolder = '/mnt/data/MATLAB\ Drive/volcanoseismo/KABrillCode/velest-matlab-scripts/';
% scriptfolder = '"D:/MATLAB Drive/volcanoseismo/KABrillCode/velest-matlab-scripts/"';
vz_RMS_mod=  [scriptfolder 'vz_RMS_mod.awk'];
vz_output= [scriptfolder 'vz_output.awk'];
vz_stadelay= [scriptfolder 'vz_stadelay.awk'];
vz_stadelayed= [scriptfolder 'vz_stadelayed.awk'];
% vz_input='/local/kabrill_res/velest_working/scripts/vz_put.awk';
% vz_EQ='/local/kabrill_res/velest_working/scripts/vz_EQ.awk';
% ex_sim='/local/kabrill_res/velest_working/scripts/extract_simultaneousloc.awk';

folderstruct=dir([parentdir velmoddir]); % Get directory structure
foldercheck=char({folderstruct.name}.') == 'm'; % Only folders beggining m
folderstruct=folderstruct(foldercheck); %account for . and .., DS_store...
nummod=length(folderstruct); % Get number of models
modfolder={folderstruct.name}'; % Get list of model folders

%% Access .OUT file and import variables to MATLAB using AWK scripts
rmsfinal=zeros(nummod,1);
velmods=cell(nummod,1);
eqlatout=cell(nummod,1);
eqlonout=cell(nummod,1);
eqdepout=cell(nummod,1);
gap=cell(nummod,1);
pweights=cell(nummod,1);

fprintf(1,'Reading Model Number     '); %used for progress notification
for ind = 1:nummod
    stacorname=[parentdir velmoddir modfolder{ind} filesep stafilename];
    if ~exist(stacorname,'file')
%         disp('.sta file does not exist')
        continue
    end
    % Get RMS from Models - AWK
    filestring = [parentdir velmoddir modfolder{ind} filesep outfilename];
    commandstringrms = ['awk -f ' vz_RMS_mod ' ' filestring];
    [~,rmsval]=system(commandstringrms);%,'-echo');
    rmsfinal(ind,1)=str2double(rmsval);
    
    % Get Velocity Models - AWK
    commandstringvelmod = ['awk -f ' vz_output ' ' filestring ];%' | awk  ''{ print $2,$1 }'];
 	[~,velmod]=system(commandstringvelmod);%,'-echo');		
    velmods{ind,1}=str2num(velmod); %#ok<ST2NM>
    
    % Get Station Corrections
    if ~exist('stlat','var')
        [stationname,ptcor1,stlat,stlon,stele] = ...
        vel_read_sta([parentdir velmoddir modfolder{ind} filesep stafilename],printvar);
        ptcor=zeros(length(stationname),nummod);
        ptcor(:,ind)=ptcor1;
    else 
        [stationname,ptcor(:,ind)] = ...
        vel_read_sta([parentdir velmoddir modfolder{ind} filesep stafilename],printvar);
    end
    %correct for ptcors greater than 10
    test=isnan(ptcor(:,ind));
    if sum(test)>0
        commandstringdelay = ['awk -f ' vz_stadelay ' ' filestring ];
        commandstringsta = ['awk -f ' vz_stadelayed ' ' filestring ];
        [~,pcorfix]=system(commandstringdelay);
        pcorfix=str2num(pcorfix); %#ok<ST2NM>
        [~,pcorfixed]=system(commandstringsta);
        pcorfixed=strsplit(pcorfixed,'\n');
        pcorfixed=pcorfixed(~cellfun(@isempty, pcorfixed))';
        badsta=stationname(test);
        for in=1:length(badsta)
            ptcor(strcmp(badsta(in,1),pcorfixed),ind)=...
                pcorfix(strcmp(pcorfixed,badsta(in,1)));
        end
    end

    % Get Final Eq Locations, Event GAPs and Pick Weights
    [~,eqlatout{ind,1},eqlonout{ind,1},eqdepout{ind,1},~,gap{ind,1},pweights{ind,1}]=...
        readcnv([parentdir velmoddir modfolder{ind} filesep cnvfilename],printvar);
    
    % Print Status
    fprintf(1,'\b\b\b\b%4.4d',ind); pause(.000001)
end
fprintf('\n')

%% Clean up

% Clean up for failed models
good=rmsfinal~=0;
eqdepout=eqdepout(good);
eqlatout=eqlatout(good);
eqlonout=eqlonout(good);
folderstruct=folderstruct(good);
gap=gap(good);
modfolder=modfolder(good);
% nummod=length(folderstruct);
ptcor=ptcor(:,good);
pweights=pweights(good);
rmsfinal=rmsfinal(good);
velmods=velmods(good);

% Clean up for unreasonable velocities
good=zeros(size(velmods));
for jnd=1:length(velmods)
    if velmods{jnd,1}(1,2)>=.2
        good(jnd,1)=1;
    end
end
good=logical(good);
eqdepout=eqdepout(good);
eqlatout=eqlatout(good);
eqlonout=eqlonout(good);
folderstruct=folderstruct(good);
gap=gap(good);
modfolder=modfolder(good);
nummod=length(folderstruct);
ptcor=ptcor(:,good);
pweights=pweights(good);
rmsfinal=rmsfinal(good);
velmods=velmods(good);


%% Pull out variables for plotting

% RMS
thresh=0.01; %percent of models you want to see
[threshget,rms_rank]=sort(rmsfinal); %get an order for best to worst RMS
threshset=ceil(thresh*nummod); %define how many models make the cut
threshrms=threshget(threshset+1); %isolate models above threshold
best_rms_ind=rms_rank(1:threshset); %get the indices of models above threshold

% Velocity Models
best_mods=modfolder(best_rms_ind); %folder name of best models
best_mods_rms=rmsfinal(best_rms_ind); %rms values of best models
best_vel_mods=velmods(best_rms_ind);%actual best models

layer_depths=best_vel_mods{1,1}(:,1); %get layer depths values
best_vel_mod_num=reshape(cell2mat(best_vel_mods),...
    [length(layer_depths),... 
    size(best_vel_mods,1)*size(best_vel_mods{1,1},2)]);
best_vel_mod_num=best_vel_mod_num(:,threshset+1:end);%step to get velocities out

best_vel_mod_mean=mean(best_vel_mod_num,2);
best_vel_mod_median=median(best_vel_mod_num,2);
best_vel_mod_mode=mode(best_vel_mod_num,2);

% Station Corrections
ptcorplot=reshape(ptcor,[size(ptcor,1)*size(ptcor,2),1]);
ptcorplotbesttmp=ptcor(:,best_rms_ind);
ptcorplotbest=reshape(ptcorplotbesttmp,[size(ptcorplotbesttmp,1)*...
    size(ptcorplotbesttmp,2),1]);
stationnameplot=repmat(cell2mat(stationname),[size(ptcor,2),1]);
stanameplotbest=repmat(cell2mat(stationname),[size(ptcorplotbesttmp,2),1]);

%% Plots

% Histogram of Model RMS
if onefig
    STATS=figure('units','normalized','outerposition',[0 0 1 1]);
    MHIST=subplot(2,4,1);
else
    MHIST=figure(1);
end

hist(rmsfinal,100)% histogram of rms values
xlabel('RMS Error')
ylabel('Number of Models')
title('RMS Values for Velocity Models')

% Histogram of pick weights
if onefig
    PWH=subplot(2,4,5);
else
    PWH=figure(2);
end

numpweights=[zeros(1,pweights{1,1}(1,1)),ones(1,pweights{1,1}(1,2)),...
    2*ones(1,pweights{1,1}(1,3)),...
    3*ones(1,pweights{1,1}(1,4)),4*ones(1,pweights{1,1}(1,5))];
if sum(pweights{1,1})~=length(numpweights)
    error('something went wrong with getting pick weights')
end
histogram(numpweights)
set(gca,'XTick',[0,1,2,3,4])
set(gca,'XLim',[-1, 5])
xlabel('Pick Weights')
ylabel('Number of Picks')
title(['Pick Weights of ' num2str(length(numpweights)) ' Picks'])

% Plot Velocity Models
% All Velocity Models shaded by Error from lowest to highest
if onefig
    AVM=subplot(2,4,2);
else
    AVM=figure(3);
end

for nnd=length(velmods):-1:1 %go backwards to put red(better) on top

    if nnd==1 %plot best velocity model in red and mark for legend
        bplot=plot(velmods{rms_rank(nnd,1),1}(:,2),...
            -velmods{rms_rank(nnd,1),1}(:,1),...
            'Color',[(nummod+1-nnd)/nummod,0,0]);
    elseif nnd==length(rms_rank) %plot worst rms velocity model in black 
        cplot=plot(velmods{rms_rank(nnd,1),1}(:,2),...
            -velmods{rms_rank(nnd,1),1}(:,1),...
            'Color',[(nummod+1-nnd)/nummod,0,0]);    
    else
        aplot=plot(velmods{rms_rank(nnd,1),1}(:,2),...
            -velmods{rms_rank(nnd,1),1}(:,1),...
            'Color',[(nummod+1-nnd)/nummod,0,0]);
    end
    hold on
end
hold off
xlim([0,6])
xlabel('Velocity (km/s)')
ylim([-10,-layer_depths(1)])
% ylim([-layer_depths(length(layer_depths)-1,1)-5,-layer_depths(1,1)])
ylabel('Depth (km)')
legend([bplot,cplot],{'Low RMS','High RMS'},'Location','ne')
title({'All Velocity Models'; 'Shaded by RMS Value'})

%Top velocity models with Mean, Median, and Mode from those 10
if onefig
    TVM=subplot(2,4,6);
else
    TVM=figure(4);
end

for jnd=length(best_vel_mods):-1:1
    rplot=plot(best_vel_mods{jnd,1}(:,2),-best_vel_mods{jnd,1}(:,1),...
        'Color',...
        [((length(best_vel_mods)+1-jnd)/length(best_vel_mods)),0,0]);
    hold on
end
gplot=plot(best_vel_mod_mean,-layer_depths,'g');
cplot=plot(best_vel_mod_median,-layer_depths,'b');
yplot=plot(best_vel_mod_mode,-layer_depths,'c');
hold off
xlim([0,6])
xlabel('Velocity (km/s)')
ylim([-10,-layer_depths(1)])
% ylim([-layer_depths(length(layer_depths)-1,1)-5,-layer_depths(1,1)])
ylabel('Depth (km)')
legend([rplot,gplot,cplot,yplot],{'Top 1%','Mean','Median','Mode'},'Location','ne')
title('Best 1% of Velocity Models')
pause(0.1)

%Station corrections of all runs
if onefig
    SCA=subplot(2,4,3);
else
    SCA=figure(5);
end

boxplot(ptcorplot,stationnameplot);
sca=get(gca);
xlabel('Station Names')
set(gca,'XTickLabelRotation',45);
ylabel(['Corrections Applied Over ' num2str(nummod) ' Runs']);
title({'Station Corrections Applied'; 'in All Runs'})
scaxlim=get(gca,'XLim');
scaylim=get(gca,'YLim');
pause(0.1)

%Station corrections for lowest RMS
if onefig
    SCB=subplot(2,4,7);
else
    SCB=figure(6);
end

boxplot(ptcorplotbest,stanameplotbest);
scb=get(gca);
xlabel('Station Names')
set(gca,'XTickLabelRotation',45);
ylabel(['Corrections Applied Over Best ' num2str(thresh*100) '% of Runs']);
title({'Station Corrections Applied'; 'for Best Runs'})
set(gca,'XLim',scaxlim);
set(gca,'YLim',scaylim);
pause(0.1)

%GAP histogram
if onefig
    GAP=subplot(2,4,4);
else
    GAP=figure(7);
end

histogram(cell2mat(gap),180)
xlabel('GAPs')
ylabel('Number of Events')
title({'GAPs from all events,'; 'all models'})
pause(1)

%GAP histogram for iterations with lowest RMS
if onefig
    GAPB=subplot(2,4,8);
else
    GAPB=figure(8); %#ok<*UNRCH>
end

histogram(cell2mat(gap(best_rms_ind)),180)
xlabel('GAPs')
ylabel('Number of Events')
title({'GAPs from all events,'; 'models with lowest RMS'})


%% Save Plots

if printplot == 1 && onefig ==1;
%     STATS.Position=[-1246,-25,1680,956];
    img = getframe(STATS);
    imwrite(img.cdata, [parentdir runtitle 'statsplot.png']);
    
end

%% EQ Location Tracking

EQP=figure('units','normalized','outerposition',[0 0 1 1]);

%setup plot
[xm,ym,zm]=xyzread(xyzfile);
mX=unique(xm); 
mY=unique(ym); mY=flipud(mY);

% mX=xm(1:length(find(xm==(xm(1)))));
% mY=ym(1:length(find(ym==(ym(1)))):end);
mZ=reshape(zm+zshift,length(mX),length(mY))';

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
seqplot=plot3(elon,elat,-edep*1000+zshift,'g.','MarkerSize',20);

%plot event locations from best rms models
for eqn=length(best_mods):-1:1
[~,e1lat,e1lon,e1dep,~,~,~]=readcnv([parentdir velmoddir best_mods{eqn,1} filesep cnvfilename]);
plot3([elon';e1lon(1:length(elon),1)'],[elat';e1lat(1:length(elat),1)'],-[edep';e1dep(1:length(edep),1)']*1000+zshift,':k')
eeqplot=plot3(e1lon(1:length(elon)),e1lat(1:length(elon)),-e1dep(1:length(elon))*1000+zshift,'.','Color',[(length(best_mods)+1-eqn)/length(best_mods),0,0],'MarkerSize',(21-eqn*2));
end

%set(gca,'XLim',[-90.9329 -90.8287],'YLim',[14.4238 14.5280])
set(gca,'XLim',[-90.9329 -90.8287],'YLim',[14.4238 14.5280],'Zlim',[-5000,4000])% Fuego zoom in
hold off
axis square
legend([staplot,seqplot,eeqplot],{'Station (Red >= 0)',...
    ['                    ';'Starting Earthquakes';...
    '                    '],...
    ['Earthquakes from best RMS models';...
    'Red is best model by RMS        ';...
    'Black is last best model by RMS ']},'Location','SouthOutside')
title('Earthquake locations of best 1% random models by RMS')

%% Save Plots

if printplot == 1;
    
    img2 = getframe(EQP);
    imwrite(img2.cdata, [parentdir runtitle 'eqplot.png']);
    view([90,0])
    img3 = getframe(EQP);
    imwrite(img3.cdata, [parentdir runtitle 'eqplotfromE.png']);
    view([0,0])
    img4 = getframe(EQP);
    imwrite(img4.cdata, [parentdir runtitle 'eqplotfromS.png']);
    view([0,90])
    img5 = getframe(EQP);
    imwrite(img5.cdata, [parentdir runtitle 'eqplotfromT.png']);
    
end

%% Print information on best mods

%Closest models to the mean median and mode of the top 1%
[~,imean]=min(sum(abs(best_vel_mod_num-repmat(best_vel_mod_mean,1,size(best_vel_mod_num,2)))));
[~,imedian]=min(sum(abs(best_vel_mod_num-repmat(best_vel_mod_median,1,size(best_vel_mod_num,2)))));
[~,imode]=min(sum(abs(best_vel_mod_num-repmat(best_vel_mod_mode,1,size(best_vel_mod_num,2)))));
[~,iminsta]=min(sqrt(range(ptcor)'.^2+rmsfinal.^2));

fprintf(['The model with the lowest RMS is ' best_mods{1,1} ' at ' num2str(best_mods_rms(1,1)) '\n'])
fprintf(['The model closest to the mean of the best 1%% of models is ' best_mods{imean,1} ' at ' num2str(best_mods_rms(imean,1)) '\n'])
fprintf(['The model with lowest station correction range and rms is ' modfolder{iminsta,1} ' at ' num2str(rmsfinal(iminsta,1)) '\n'])
fprintf(['Spread of Best 1%% of RMS is ' num2str(best_mods_rms(end,1)-best_mods_rms(1,1)) '\n'])
fprintf(['Difference between lowest station correction range and rms and best rms is ' num2str(rmsfinal(iminsta,1)-best_mods_rms(1,1)) '\n'])
disp([best_mods{1,1} '       ' best_mods{imean,1} '       ' modfolder{iminsta,1}])
iminstanum=velmods{iminsta,1};
iminstanum=iminstanum(1:2:end,2);
disp([best_vel_mod_num(1:2:end,[1,imean]),iminstanum])
fprintf('Station corrections for lowest station correction range and rms: \n')
disp(stationname')
fprintf([num2str(ptcor(:,iminsta)') '\n'])

%% Run rvelmodevalsup
rvelmodevalsup