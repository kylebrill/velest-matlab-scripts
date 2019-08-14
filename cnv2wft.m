% function w = cnv2wfo(cnvfile,dbpath)
%CNV2WFO - get waveform objects from velest cnv files
%   W = CNV2WFT(CNVFILE) returns GISMO waveform object events from VELEST 
%   .CNVFILEs. 

dbpath='/local/kabrill_res/Fuego/Fuego_2012/db/goodpicks-mag/goodpicks';
cnvfile='/local/kabrill_res/velest_acamrecipereorg/velestdata/fuegonex-rand01m631.CNV';

[origtime,~,~,~,~,~,~]=readcnv(cnvfile);

db = dbopen(dbpath,'r');
dbor=dblookup_table(db,'origin');
orids=dbgetv(dbor,'orid');
otimes=dbgetv(dbor,'time');
otimes=datenum(strtime(otimes));

% ismember(otimes,origtime);

% dbassoc=dblookup_table(db,'assoc');
% dbarr=dblookup_table(db,'arrival');

% starttimes=artimes;

dbsites=dblookup_table(db,'site'); %open site table
allsta=dbgetv(dbsites,'sta');%cell array of station names

ds=datasource('antelope','/local/kabrill_res/Fuego/Fuego_2012/db/goodpicks/goodpicks');
scnl=scnlobject(allsta,'HHZ','XX','');
w=waveform(ds,scnl,origtime-datenum(0,0,0,0,5,00),origtime+datenum(0,0,0,0,10,0));
fo=filterobject('B',[1,10],2);
w=filtfilt(fo,w);
for in=1:length(w)
    plot(w(in))
    pause
end

