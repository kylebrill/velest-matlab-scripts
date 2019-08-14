%ADJUSTPICKS uses residuals to identify and correct bad picks
% Retrieves event times and picks from .CNV
% Retrieves residuals from .OUT
% Loads GISMO waveforms from datasource
% Multiple channels multiple stations for picking
% Rearrange based on explosions vs events
% Save new picks and prepare to rerun velest.
% 
% By Kyle Brill 
% Last Edited 23 October, 2016

%% Inputs
%.CNV file
% cnvfile='SW2irand.CNV';L
cnvfile='NW1rand1rev.CNV';
shotfile='NW1rand1sh.CNV';
% cnvfile='SW1irandadjrev.CNV';

%.OUT file
% outfile='SW2irand.OUT';
outfile='NW1rand1.OUT';


%.CNV file out
cnvfileout='outrand1adj.CNV';
%PWS
numPWS=5;
%max acceptable residual
maxres=0.25;

% root directory of velest data
% resdirpath=whichresdirpath;
% rootdir=[resdirpath 'velest_eqshpws' filesep 'data' filesep];
rootdir='M:\kabrill_res\velest_adjeqshpws\data\';

%datasource
fuego2012path=whichfuego2012;
ds=datasource('file',[fuego2012path '%s/%s%s%04d%02d%02d.mat'],...
     'station','station','channel','year','month','day');
%% Open .CNV file and retrieve event times and picks

[orig,lat,lon,depth,mag,gap,rms,wts,picks]=readcnv([rootdir cnvfile]);
eqnumcnv=length(orig);
if exist('shotfile','var')
[sho,shla,shlo,shd,shm,shg,shr,shw,shp]=readcnv([rootdir shotfile]);
shnumcnv=length(sho);
orig=[orig;sho]; lat=[lat;shla]; lon=[lon;shlo]; depth=[depth;shd];
mag=[mag;shm];gap=[gap;shg]; rms=[rms;shr]; wts=wts+shw; picks=[picks;shp];
clear sho shla shlo shd shm shg shr shw shp
end

%% Open .OUT file and retrieve residuals for each pick at each station
[event,shotnumout,eqnumout]=get_sta_res([rootdir outfile]);

%% Trim events if needed to match input eq and sh
if eqnumcnv+shnumcnv ~= eqnumout+shotnumout
    eventkeep=ones(length(event),1);
    for ie=numPWS+1:length(event)
        if sum([event(ie).origtime]==orig) ~= 1
            eventkeep(ie)=0;
        end
    end
    eventkeep=logical(eventkeep);
    event=event(eventkeep,1);
end
%% Waveform object details
secbefore=15;
secafter=45;
spad=datenum(0,0,0,0,0,secbefore);
epad=datenum(0,0,0,0,0,secafter);
fo=filterobject('B',[1,25],2);
so=spectralobject(512,(512*0.8),25,[50 125],'s');
nw='--';
lo='--';
cha={'HHZ','HHN','HHE'};
numcha=length(cha);
updatepicks=picks;

% load('D:\Google Drive\velest_adjeqshpws\updatepicksprog.mat','updatepicks')

%% Start big loop over all events
% for ev=numPWS+1:length(event)
for ev=1:length(event)
%     ev=numPWS+1;
%% Find large residuals to address
stas=updatepicks(ev,1).sta;
res=zeros(length(stas),1);
for in=1:length(stas)
    res(in,1)=event(ev,1).(stas{in,1}).res;
end
fprintf('Working on event #%04d\n',ev)
%% If a residual, go into picking loop, otherwise leave alone

if sum(abs(res)>maxres)>0
    %% Get Current Picks
    cnvpks=secbefore+picks(ev,1).obt;
    cnvpickmatrix=nan(length(stas),6);
    cnvpickmatrix(1:length(stas),1)=1:1:length(stas);
    cnvpickmatrix(1:length(stas),2)=cnvpks;    
    
    %% Setup ChannelTags
    numsta=length(stas);
    sta=strrep(stas,'_','')'; %get rid of underscores in station names
    focuscheck=sta(abs(res)>maxres); %indices of stations with high residuals
    %reshaping for ChannelTags
    sta=repmat(sta,numcha,1);
    sta=reshape(sta,1,numsta*numcha);
    chan=repmat(cha,1,numsta);
    ct=ChannelTag.array(nw,sta,lo,chan);
    %% Make waveform object
    wo=waveform(ds,ct,orig(ev)-spad,orig(ev)+epad);
    wo=fix_data_length(wo);
    wo=demean(wo);
    wo=filtfilt(fo,wo);
    sps=get(wo(1,1),'freq');
    waves=wo(ct.matching([],[],[],'HHE'));%vertical channels for picking
    
    %% Skip to next if we have already picked these events
    pmf=sprintf('%02d',ev);
    
    if ~exist([pmf '.txt'],'file')
    %% Show all channels for stations with high residuals and picks
    focusgrab=zeros(size(wo));
    for kn=1:length(focuscheck)
        focusgrab=focusgrab+ct.matching([],focuscheck{1,kn},[],[]);
    end
    focusgrab=logical(focusgrab);
    helpwave=wo(focusgrab);%all channels for stations with resid > maxres
    plot_panels(helpwave,true);
    pp=gcf;
    pp.Position=[962 52 958 932];
    
    % add current picks??
    hAllAxes = findobj(pp,'type','axes');
    for mn=1:length(hAllAxes)
        if strfind(hAllAxes(mn,1).YLabel.String{2,1},'HHZ')
            whcsta=~cellfun(@isempty,strfind(stas,hAllAxes(mn,1).YLabel.String{1,1}));
            axes(hAllAxes(mn,1)); %#ok<LAXES>
            line([cnvpks(whcsta),cnvpks(whcsta)],get(gca,'YLim'))%,'LineWidth',2)
        end
    end
            
    
    %% Pick arrivals
    pickcts=get(waves,'ChannelTag');
    refdisp=[num2cell(1:1:length(res))',{pickcts.station}',num2cell(res)];
    disp(refdisp)
    pickmatrix=seis_pick(cell2mat(get(waves,'data')),1/sps,1);
    
    else
        pickmatrix=load([pmf '.txt'],'-ascii');
    end
    
    %% Get arrivals and errors for each pick from pickmatrix
    aperr=(pickmatrix(:,4)-pickmatrix(:,2))/2; %error=spread from midpoint
    updateobt=picks(ev,1).obt;
    stas=picks(ev,1).sta;
    phase=picks(ev,1).pha;
    wt=picks(ev,1).wgt;
    for nn=1:length(aperr)
        if ~isnan(aperr(nn,1))
            newpicks=pickmatrix(nn,2)+aperr(nn,1); %add midpoint to first to get arr
            
            %get absolute arrival time
            UTCpicks=get(waves(1,nn),'start')+datenum(0,0,0,0,0,newpicks);
            %get time relative to origin for new obt
            chksign=sign(UTCpicks-orig(ev));
            abstime=datestr(abs(UTCpicks-orig(ev)),'SS.FFF');
            updateobt(nn,1)=chksign*str2num(datestr(abs(UTCpicks-orig(ev)),'SS.FFF')); %#ok<ST2NM>
            
            %% Set up phase and weight information to pass to new .CNV
            %phase
            phase(nn,1)={'P'};
            
            %adjust pick weight - this is the scheme used by UUSS
            wt(nn,1)=0;
            if abs(aperr(nn))<=0.06
                wt(nn,1)=0;
            elseif abs(aperr(nn))<=.12
                wt(nn,1)=1;
            elseif abs(aperr(nn))<=.30
                wt(nn,1)=2;
            elseif abs(aperr(nn))<=.60
                wt(nn,1)=3;
            else
                wt(nn,1)=4;
            end
        end
    end
    updatepicks(ev,1).obt=updateobt;
    updatepicks(ev,1).sta=stas;
    updatepicks(ev,1).pha=phase;
    updatepicks(ev,1).wgt=wt;
end

closereq

end
%% Write out the new CNV that includes the new PWS picks and excludes the events used for the PWS
fprintf('\n')
fid=fopen([rootdir cnvfileout],'w');
accepted = 0;
for nevt=1:length(updatepicks)
% write out the hypocenter info
%(3i2,1x,2i2,1x,f5.2,1x,f7.4,a1,1x,f8.4,a1,1x,f7.2,2x,f5.2)
%a4,a1,i1,f6.2
%      840210 23 6 50.92 37.1653N 121.5500W    5.62   0.00     73  0.0 0.03  1.0  1.0
%      HSPMP0  1.31CADMP1  1.55HGSMP0  1.97JCBMP1  2.36CCOMP2  3.07HCAMP1  2.79
%      HGWMP1  3.19CAOMP1  3.62JSTMP1  3.83HCRMP1  3.83HFEMP0  3.98JHLMP2  4.44
%      HPLMP0  4.31JALMP0  4.53
% convert from epoch time to string to datenum to datevec
        
        [yr,mo,da,hr,mn,sc]=datevec(datenum(orig(nevt)));
        fprintf(fid,'%02d%02d%02d %02d%02d %05.2f %7.4fN %8.4fW%7.2f  %5.2f    %3d     %5.2f\n',yr-2000,mo,da,hr,mn,sc,lat(nevt),-1*lon(nevt),depth(nevt),mag(nevt),gap(nevt),rms(nevt));
        accepted = accepted + 1; %Count number of events printed
        % write out the pick info
        for npk=1:length(updatepicks(nevt).obt)
            %make sure stations are padded to 4 characters long
            tmpsta=char(updatepicks(nevt).sta{npk});
            if length(tmpsta)>4, updatepicks(nevt).sta{npk}(1:4)=tmpsta(1:4);end
            if length(tmpsta)==1, updatepicks(nevt).sta{npk}(1:4)='____';updatepicks(nevt).sta{npk}(1)=tmpsta(1);end
            if length(tmpsta)==2, updatepicks(nevt).sta{npk}(1:4)='____';updatepicks(nevt).sta{npk}(1:2)=tmpsta(1:2);end
            if length(tmpsta)==3, updatepicks(nevt).sta{npk}(1:4)='____';updatepicks(nevt).sta{npk}(1:3)=tmpsta(1:3);end
            if length(tmpsta)==4, updatepicks(nevt).sta{npk}(1:4)=tmpsta(1:4);end
            
            if npk==7 || npk==13 || npk==19
                fprintf(fid,'\n');
            end

            fprintf(fid,'%4s%1s%1d%6.2f',updatepicks(nevt).sta{npk},updatepicks(nevt).pha{npk},updatepicks(nevt).wgt(npk),updatepicks(nevt).obt(npk));
            
            if npk==length(updatepicks(nevt).obt)
                fprintf(fid,'\n\n');
            end
        end
end

fclose(fid);