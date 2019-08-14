function [varargout]=readcnv(varargin)
%READCNV reads information from VELEST earthquake .CNV files
%   [ORIG,LAT,LON,DEPTH,MAG,GAP,RMS,WTS,PICKS] =
%       READCNV(CNVFILE) opens a .CNV file made for or produced by VELEST
%       and outputs event ORIGINS, LATITUDES, LONGITUDES, MAGNITUDES, GAPS,
%       RMS, number of pick WEIGHTS 1-4, and information about PICKS, like
%       station name (PICKS.STA), phase (PICKS.PHA), weight (PICKS.WGT),
%       and observation time (PICKS.OBT). If a given field is unwanted, use
%       ~ as a placesaver.
%
% Example: [orig,lat,lon,depth,mag,gap,rms,wts,picks]=readcnv('hypo.CNV');
%
%   [ORIG,LAT,LON,DEPTH,MAG,GAP,RMS,WTS,PICKS] =
%       READCNV(CNVFILE,NUMEVENTS) allows inputing the number of events in
%       the .cnv file to allow for preallocating. 
%
% Example: [orig,lat,lon,depth,mag,gap,rms,wts,picks] =
%               readcnv('hypo.CNV',127);
%
%   [ORIG,LAT,LON,DEPTH,MAG,GAP,RMS,WTS,PICKS] =
%       READCNV(CNVFILE,PRINTFLAG) allows for output files to be produced
%       in each directory which are useful for plotting with GMT or reading
%       with AWK. A total of 2 .dat files can be produced for each .CNV
%       input file:
%
%           1) Hypocenters in a "locations_raw_out.dat" file
%           2) GAPs for each event in a "GAP.dat" file
%       
%       Setting printflag = 1 prints both for each .CNV input file.
%       Printflag can also be a two element row vector of 1s and 0s where
%       1s mean that file will be created and 0s omit specific file
%       creation. The order of the elements determines which plots are
%       created.
%
% Example:  [orig,lat,lon,depth,mag,gap,rms,wts,picks] = 
%                 readcnv('hypo.CNV',1); 
%       Prints all three .dat files.
%
% Example:  [orig,lat,lon,depth,mag,gap,rms,wts,picks] = 
%                 readcnv('hypo.CNV',[1,0]);
%       Prints Hypocenter location files only.
%
%   [ORIG,LAT,LON,DEPTH,MAG,GAP,RMS,WTS,PICKS] =
%       READCNV(CNVFILE,PRINTFLAG,EVENTNUM) or (CNVFILE,EVENTNUM,PRINTFLAG)
%       allows inputing the number of events in the .cnv file to allow for 
%       preallocating as well as printing .dat files described above.
%
% Written and last edited by Kyle Brill 2016/08/09 12:00 EDT

%% Initialize variables.fileID = fopen(cnvfile,'r');
if ~exist('varargin','var')
    error('Input variables needed.')
end

%Check 1st input variable
if ~ischar(varargin{1})
    error('First input must be .CNV file string.')
end

%Make sure 1st input variable is loaded properly and give defaults.
cnvfile=varargin{1};
numevents=1;
printflag=[0,0];

fileID = fopen(cnvfile,'r');
finfo = dir(cnvfile);
if isempty(finfo)
    error('Unable to find .CNV file as entered.')
end

% Check 2nd input variable
if nargin == 2
    if ischar(varargin{2})
        error('Second input should either be a print flag or a scalar number of events.')
    end
    if length(varargin{2}) == 1
        if varargin{2} > 1
            printflag=[0,0];
            numevents=varargin{2};
        elseif varargin{2} == 1
            printflag=[1,1];
            numevents=1;
        elseif varargin{2} == 0;
            printflag=[0,0];
            numevents=1;
        end
    elseif length(varargin{2}) == 2
        printflag = varargin{2};
        numevents=1;
    else 
        error('Second input should either be a print flag or a scalar number of events.')
    end
end

%Check 2nd and 3rd input variables.
if nargin == 3
    if ischar(varargin{2}) || ischar(varargin{3})
        error('Second and third inputs should either be a print flag or a scalar number of events.')
    end
    
    if length(varargin{2}) > length(varargin{3})
        
        if length(varagin{2}) == 2 && sum(varargin{2}) <= 2 && sum(varargin{2}) >= 0
            printflag = varargin{2};
        else
            error('Print Flag (input 2) is incorrect.')
        end
        
        if length(varargin{3}) == 1 && isscalar(varargin{3})
            numevents = varagin{3};
        else
            error('Number of events (input 3) is incorrect.')
        end
        
    elseif length(varargin{2}) < length(varargin{3})
        if length(varagin{3}) == 2 && sum(varargin{3}) <= 2 && sum(varargin{3}) >= 0
            printflag = varargin{3};
        else
            error('Print Flag (input 3) is incorrect.')
        end
        
        if length(varargin{2}) == 1 && isscalar(varargin{2})
            numevents = varagin{2};
        else
            error('Number of events (input 2) is incorrect.')
        end
        
    elseif length(varargin{2}) == length(varargin{3})
        if varargin{2} > 1 && varargin{3} <= 1
            numevents=varargin{2};
            if varargin{3} == 1
                printflag = [1,1];
            elseif varargin{3} == 0
                printflag = [0,0];
            end
        elseif varargin{3} > 1 && varargin{2} <= 1
            numevents=varargin{3};
            if varargin{2} == 1
                printflag = [1,1];
            elseif varargin{2} == 0
                printflag = [0,0];
            end
        elseif varargin{2} > 1 && varargin{3} > 1
            error('Unacceptable inputs')
        end
        
     end
end
    
if nargin > 3
    error ('Check inputs, something is screwy.')
end

if nargout < 7
    error('Must have 7 output values. Use ~ as placesaver for unwanted outputs.')
    % This would be good to change to a case switch output sometime so it
    % could handle whatever output was required...
end

if ~exist('printflag','var') || ~exist('numevents','var')
    error('Need to debug code above this error.')
end

%% Preallocate Variables
orig=zeros(numevents,1);
lat=zeros(numevents,1);
lon=zeros(numevents,1);
depth=zeros(numevents,1);
mag=zeros(numevents,1);
gap=zeros(numevents,1);
rms=zeros(numevents,1);
picks(numevents,1).obt=[];

%% Read file into MATLAB

eventnum=1;

tmpline=fgetl(fileID);
while ischar(tmpline)
    if ~isempty(tmpline)
        if isstrprop(tmpline(1),'digit') %read event location info
            
            year2 = str2double(tmpline(1:2));
            month = str2double(tmpline(3:4));
            day = str2double(tmpline(5:6));
            hour = str2double(tmpline(8:9));
            mn = str2double(tmpline(10:11));
            sec = str2double(tmpline(13:17));
            
            orig(eventnum,1) = datenum(year2+2000,month,day,hour,mn,sec);
            
            latt = str2double(tmpline(19:25));
            latc= tmpline(26);
            if strcmp(latc,'S')
                latt=-latt;
            end
            lat(eventnum,1)= latt;
            
            lonn = str2double(tmpline(28:35));
            lonc= tmpline(36);
            if strcmp(lonc,'W')
                lonn=-lonn;
            end
            lon(eventnum,1)= lonn;
            
            depth(eventnum,1) = str2double(tmpline(37:44));
            
            mag(eventnum,1) = str2double(tmpline(46:51));
            
            gap(eventnum,1) = str2double(tmpline(55:57));
            
            rms(eventnum,1) = str2double(tmpline(62:67));
%             ifix = str2double(tmpline(63:67));            
            
        else %read pick info
            ll=length(tmpline);
            if ~exist('numpicks','var')
                numpicks=ll/12;
                for ps=1:numpicks
                    sh=(ps-1)*12;
                    picks(eventnum,1).sta{ps,1}=tmpline(1+sh:4+sh);
                    picks(eventnum,1).pha{ps,1}=tmpline(5+sh);
                    picks(eventnum,1).wgt(ps,1)=str2double(tmpline(6+sh));
                    picks(eventnum,1).obt(ps,1)=str2double(tmpline(8+sh:12+sh));
                end
            else
                numpicks=ll/12;
                for pt=1:numpicks
                    sh=(pt-1)*12;
                    picks(eventnum,1).sta{pt+ps,1}=tmpline(1+sh:4+sh);
                    picks(eventnum,1).pha{pt+ps,1}=tmpline(5+sh);
                    picks(eventnum,1).wgt(pt+ps,1)=str2double(tmpline(6+sh));
                    picks(eventnum,1).obt(pt+ps,1)=str2double(tmpline(8+sh:12+sh));
                end
            end
            
            
            
        end
    else
        eventnum=eventnum+1;
        clear numpicks
    end
    tmpline=fgetl(fileID);
end
fclose(fileID);

%% Tally up Pick Weights

wts=[0,0,0,0,0];

for in=1:length(picks)
    P0=sum(picks(in,1).wgt==0);
    P1=sum(picks(in,1).wgt==1);
    P2=sum(picks(in,1).wgt==2);
    P3=sum(picks(in,1).wgt==3);
    P4=sum(picks(in,1).wgt==4);
    
    wts(1,1)=wts(1,1)+P0;
    wts(1,2)=wts(1,2)+P1;
    wts(1,3)=wts(1,3)+P2;
    wts(1,4)=wts(1,4)+P3;
    wts(1,5)=wts(1,5)+P4;
end

%% Outputs

if nargout == 7
    varargout{1}=orig;
    varargout{2}=lat;
    varargout{3}=lon;
    varargout{4}=depth;
    varargout{5}=mag;
    varargout{6}=gap;
    varargout{7}=wts;
elseif nargout== 8
    varargout{1}=orig;
    varargout{2}=lat;
    varargout{3}=lon;
    varargout{4}=depth;
    varargout{5}=mag;
    varargout{6}=gap;
    varargout{7}=wts;
    varargout{8}=picks;
elseif nargout== 9
    varargout{1}=orig;
    varargout{2}=lat;
    varargout{3}=lon;
    varargout{4}=depth;
    varargout{5}=mag;
    varargout{6}=gap;
    varargout{7}=rms;
    varargout{8}=wts;
    varargout{9}=picks;
end

%% Printing Outputs

if sum(printflag)>0
    
    eventnum=eventnum-1;
    event=(1:1:eventnum)';
    
    if printflag(1,1) == 1
        outlocfilename=strrep(cnvfile,...
            cnvfile(1,end-3:end),'locations_raw_out.dat');
        
        fid=fopen(outlocfilename,'w');
        fprintf(fid,'#Event Latitude Longitude Depth\n');
        for in=1:eventnum
            if in==eventnum
                fprintf(fid,'%4d %7.4f %7.4f  %5.2f',event(in,1),...
                    lat(in,1),lon(in,1),depth(in,1));
            else
                fprintf(fid,'%4d %7.4f %7.4f  %5.2f\n',event(in,1),...
                    lat(in,1),lon(in,1),depth(in,1));
            end
        end
        fclose(fid);
    end
    
    if printflag(1,2) == 1
        gapfile=strrep(cnvfile,cnvfile(1,end-3:end),'GAP.dat');
        fid=fopen(gapfile,'w');
        fprintf(fid,'# GAPs were between %3d and %3d (average GAP was %3d)\n',...
            min(gap),max(gap),ceil(mean(gap)));
        fprintf(fid,'# Event#  GAP\n');
        for in=1:eventnum
            if in==eventnum
                fprintf(fid,'%4d %4d',event(in,1),gap(in,1));
            else
                fprintf(fid,'%4d %4d\n',event(in,1),gap(in,1));
            end
        end
        fclose(fid);
    end
    
end
        