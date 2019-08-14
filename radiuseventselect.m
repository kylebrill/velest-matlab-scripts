% pick out events which dbgenloc locates within seismic network

% run ('/opt/antelope/5.5/setup.m') %Initialize Antelope
%% Initialize Variables

%fullpath to db descriptor file
dbpath='/local/kabrill_res/Fuego/Fuego_2012/db/goodpicks-mag/goodpicks'; 

%station used to establish radius around vent
refsta='S1'; 

vollat=14.474471;%lat of fuego vent
vollon=90.880940;%lon of fuego vent
minpicks=6; % number of detecting stations needed 
mindepth=5; % depth in km of lowest earthquake locations accepted as located by antelope (no topography)
pwcut=2; %pick weight cutoff, values greater than this will not be printed

%outfile name
cnvfile='hypo6picksP2orbetter.CNV';

%% Open database and extract information
db = dbopen(dbpath,'r');
dbsites=dblookup_table(db,'site'); %open site table
allsta=dbgetv(dbsites,'sta');%cell array of station names
stalat=dbgetv(dbsites,'lat');%matrix of station latitudes
stalon=dbgetv(dbsites,'lon');%matrix of station longitudes

dbor=dblookup_table(db,'origin');
dbor_size = dbquery(dbor,'dbRECORD_COUNT');
disp(['sorting through ',num2str(dbor_size),' events'])

% origerr
dborerr = dblookup_table(db,'origerr');
orids=dbgetv(dbor,'orid');

% origin x assoc x arrival
dbassoc1 = dblookup_table(db,'assoc');
% phase=dbgetv(dbassoc1,'phase')
dbassoc = dbjoin(dbor,dbassoc1,'orid');

dbarr1 = dblookup_table(db,'arrival');
dbarr = dbjoin(dbassoc,dbarr1);
dbarr = dbsort(dbarr,'orid');
dbarr_size = dbquery(dbarr,'dbRECORD_COUNT');
disp(['with ', num2str(dbarr_size),' arrivals'])

[oo,ola,olo,ot,od,nass,mag]=dbgetv(dbor,'orid','lat','lon','time',...
    'depth','nass','ml');
[ao,at,as,ap,ach,aperr]=dbgetv(dbarr,'origin.orid','arrival.time',...
    'arrival.sta','arrival.iphase','arrival.chan','arrival.deltim');

%% Establish radius area based on max distance from vent

radlat=stalat(strcmp(refsta,allsta));
radlon=stalon(strcmp(refsta,allsta));

distmax=distance(vollat,vollon,radlat,radlon);

inside=distance(vollat,vollon,ola,olo) < distmax;

%% Establish pick weight cutoff
if pwcut == 0
    pwlim=0.06;
elseif pwcut == 1
    pwlim=0.12;
elseif pwcut == 2
    pwlim=0.30;
elseif pwcut == 3
    pwlim=0.60;
elseif pwcut == 4
    pwlim=9999;
end

%% Print Hypocenter info to file based on criteria
fid=fopen(cnvfile,'w');
a=find(oo>0);
accepted = 0;
%  oridstmp=orids(a)
for nevt=1:length(orids)
    if oo(nevt)>0
        tmp=find(oo(nevt)==ao);
        weightcut=find(aperr(tmp)<=pwlim);
        cutter=length(weightcut);
        if cutter >= minpicks %select events with number of detections
            if inside(nevt) %select events within radius of a station
                if od(nevt) < mindepth %select events within cone
                    % write out the hypocenter info
                    %(3i2,1x,2i2,1x,f5.2,1x,f7.4,a1,1x,f8.4,a1,1x,f7.2,2x,f5.2)
                    %a4,a1,i1,f6.2
                    %      840210 23 6 50.92 37.1653N 121.5500W    5.62   0.00     73  0.0 0.03  1.0  1.0
                    %      HSPMP0  1.31CADMP1  1.55HGSMP0  1.97JCBMP1  2.36CCOMP2  3.07HCAMP1  2.79
                    %      HGWMP1  3.19CAOMP1  3.62JSTMP1  3.83HCRMP1  3.83HFEMP0  3.98JHLMP2  4.44
                    %      HPLMP0  4.31JALMP0  4.53
                    % convert from epoch time to string to datenum to datevec
                    tsta=as(tmp(weightcut));
                    tpha=ap(tmp(weightcut));
                    twei=aperr(tmp(weightcut));
                    tat=at(tmp(weightcut));
                    
                    [yr,mo,da,hr,mn,sc]=datevec(datenum(strtime(ot(nevt))));
                    odatenum=datenum(strtime(ot(nevt)));
                    fprintf(fid,'%02d%02d%02d %02d%02d %05.2f %7.4fN %8.4fW %7.2f   %4.2f  360  9.9 9.99  9.9  9.9\n',yr-2000,mo,da,hr,mn,sc,ola(nevt),-1*olo(nevt),od(nevt),mag(nevt));
                    accepted = accepted + 1; %Count number of events printed
                    % write out the pick info
                    for npk=1:cutter
%                         wt=[]; diffsec=[];
                        % adjust pick time
                        [yrp,mop,dap,hrp,mnp,scp]=datevec(datenum(strtime(tat(npk))));
                        pdatenum=datenum(strtime(tat(npk)));
%                         diffsec(npk)=86400*(pdatenum-odatenum);
                        diffsec=86400*(pdatenum-odatenum);
                        % adjust pick weight
                        if twei(npk)<=0.06%0.1
%                             wt(npk)=0;
                            wt=0;
                        elseif twei(npk)<=0.12%.2
%                             wt(npk)=1;
                            wt=1;
                        elseif twei(npk)<=0.30
%                             wt(npk)=2;
                            wt=2;
                        elseif twei(npk)<=0.60%.4
%                             wt(npk)=3;
                            wt=3;
                        else
%                             wt(npk)=4;
                            wt=4;
                        end
                        phase=char(tpha(npk));
                        tmpsta=char(tsta(npk));
                        if length(tmpsta)>4, sta(1:4)=tmpsta(1:4);end
                        if length(tmpsta)==1, sta(1:4)='____';sta(1)=tmpsta(1);end
                        if length(tmpsta)==2, sta(1:4)='____';sta(1:2)=tmpsta(1:2);end
                        if length(tmpsta)==3, sta(1:4)='____';sta(1:3)=tmpsta(1:3);end
                        if length(tmpsta)==4, sta(1:4)=tmpsta(1:4);end
                        if npk==7 || npk==13
                            fprintf(fid,'\n');
                        end
%                         sprintf('%4s%1s%1d%6.2f',sta,phase,wt(npk),diffsec(npk))
%                         fprintf(fid,'%4s%1s%1d%6.2f',sta,phase,wt(npk),diffsec(npk));
%                         sprintf('%4s%1s%1d%6.2f',sta,phase,wt,diffsec)
                        fprintf(fid,'%4s%1s%1d%6.2f',sta,phase,wt,diffsec);

                        if npk==cutter
                            fprintf(fid,'\n\n');
                        end
                    end
                end
                clear sta phase wt diffsec
            end
        end
    else
        %bogus orid
    end
end
fclose(fid);
display([num2str(accepted) ' events in new .CNV'])