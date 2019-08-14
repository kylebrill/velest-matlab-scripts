function [event,shotnum,eqnum]=get_sta_res(outfilename)
%GET_STA_RES reports station residuals for each event from velest .OUT file
%   [EVENT] = GET_STA_RES (OUTFILENAME) opens velest .OUT file and returns
%   residual information for each EVENT in the file OUTFILENAME, which
%   includes both events and shots. This can be useful when trying to
%   eliminate bad data from a dataset. Output is a structure with origin
%   times and whatever information is available for each station at each
%   event time.
%
% Example: [event,shotnum,eqnum]=get_sta_res('velest.OUT');
%
%Written by Kyle Brill
%Last Edited by Kyle Brill - 11 October, 2016 11:50 EDT

%% Read file into MATLAB
fout = fopen(outfilename,'r');
outfile = textscan(fout,'%s','delimiter','\n');
outfile = outfile{1,1};
fclose(fout);

%% Load information from .OUT File

% Establish Beginning and End of Information
startRow=find(~cellfun(@isempty,strfind(outfile,'~~~ output final station residuals ...')))+8;
endRow=find(~cellfun(@isempty,strfind(outfile,'station statistics, remember nsp was set to: 1')))-10;
stares = outfile(startRow(1,1):endRow(1,1),1);

% Number of events in .OUT file
eqandshots=strsplit(outfile{38,1},' ');
eqnum=str2double(eqandshots{1});
shotnum=str2double(eqandshots{2});
eventnum=eqnum+shotnum;

% Number of stations in .OUT file
numsta=outfile{find(~cellfun(@isempty,strfind(...
    outfile,'Total number of stations with readings:'))),1}; %#ok<FNDSB>
numsta=str2double(numsta(40:end));

% Station names
stacells=outfile(...
    find(~cellfun(@isempty,strfind(outfile,'stn latitude longitude elev')))+1:...
    find(~cellfun(@isempty,strfind(outfile,'stn latitude longitude elev')))+numsta,...
    1);
staname=cell(numsta,1);
for in=1:numsta
    if in<10
        staname{in,1}=stacells{in,1}(3:6);
    else
        staname{in,1}=stacells{in,1}(4:7);
    end
end
staname=char(staname);

%% Rearrange data into useful numbers

%preallocate structure
sta=struct('ph',char(),'wt',[],'res',[],'ttime',[],'delta',[]);
event=struct('origtime',0);
for in=1:numsta
    event.(staname(in,:))=sta;
end
event=repmat(event,eventnum,1);

%fill out structure
for in=1:length(stares)
    rdr=stares{in,1};
    if ~isempty(rdr)%something is in the cell
        if strcmp(rdr(1:7),'Station') %Start of each event
            ev=str2double(rdr(29:32));
            event(ev,1).origtime=datenum(rdr(36:end),'yymmdd HHMM SS.FFF');
        elseif strcmp(rdr(1:7),'sta ph ')
            %do nothing
        else
            event(ev,1).(rdr(1:4)).ph = rdr(6);
            event(ev,1).(rdr(1:4)).wt = str2double(rdr(8:9));
            event(ev,1).(rdr(1:4)).res = str2double(rdr(11:17));
            event(ev,1).(rdr(1:4)).ttime = str2double(rdr(17:22));
            event(ev,1).(rdr(1:4)).delta = str2double(rdr(23:28));
            if length(rdr)>29
                event(ev,1).(rdr(34:37)).ph = rdr(39);
                event(ev,1).(rdr(34:37)).wt = str2double(rdr(41:42));
                event(ev,1).(rdr(34:37)).res = str2double(rdr(44:49));
                event(ev,1).(rdr(34:37)).ttime = str2double(rdr(50:55));
                event(ev,1).(rdr(34:37)).delta = str2double(rdr(56:61));
            end
        end
    end
    
end

%% Clean up
% Eliminate elements of structure with no information
for in=1:length(event)
    for jn=1:numsta
        if isempty(event(in).(staname(jn,1:4)).ttime)
            event(in).(staname(jn,1:4))=[];
        end
    end
end

%% Outputs