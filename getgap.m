function [event,gap] = getgap(outfile,printflag)
%GETGAP Imports GAP of each EVENT from final iteration of VELEST run
%   GETGAP(OUTFILE,PRINTFLAG) Extracts GAP data from .OUT file from final
%   iteration of a VELEST run and prints it to a *gap.dat file in the same
%   directory as the .OUT file if PRINTFLAG == 1
%
% Example:
%   getgap('fuegorandmods01.OUT',1)
%
%   [EVENT,GAP] = GETGAP(OUTFILE) Imports GAP data and Event Number from
%   .OUT file from final iteration of VELEST run.
%
% Example:
%   [event,gap] = getgap('fuegorandmods01.OUT')
%
% Written and last edited by Kyle Brill 2016/04/07 00:30 EDT

%% Initialize variables.
if nargin<=1
    printflag = 0;
end

%% Read file into MATLAB
fileID = fopen(outfile,'r');
% c=textread(outfile,'%s','delimiter','\n');
c=textscan(fileID,'%s','delimiter','\n');
c=c{1,1};
fclose(fileID);

%% Establish Beginning and End of Information
startRow=find(~cellfun(@isempty,strfind(c,'Event# ->')))+1;
endRow=find(~cellfun(@isempty,strfind(c,'GAPs were between')))-2;
gaptext=c(startRow:endRow,1);

%% Rearrange data into useful numbers
gapblocks=cell(length(gaptext),15);
for in=1:length(gaptext)
    tmpblocks=strsplit(gaptext{in,1},' ');%Added to avoid Subscripted assignment dimension mismatch when each row of gaps is not the same length.
    gapblocks(in,1:length(tmpblocks))=strsplit(gaptext{in,1},' ');
end

event=str2double(gapblocks(:,1:3:15));
gap=str2double(gapblocks(:,3:3:15));

[row,col]=size(event);
numevents=row*col;

event=reshape(event',[numevents,1]);
gap=reshape(gap',[numevents,1]);

event=event(~isnan(event)); %get rid of extra space at the bottom from nonsquare gap matrices
gap=gap(~isnan(gap));
numevents=length(event); %get real number of events

%% Print if desired
if printflag == 1
    gapfile=strrep(outfile,outfile(1,end-3:end),'GAP.dat');
    fid=fopen(gapfile,'w');
    fprintf(fid,['# ' c{endRow+2,1} ' ' c{endRow+3,1}]);
    fprintf(fid,'\nEvent#  GAP\n');
    for in=1:numevents
        if in==numevents
            fprintf(fid,'%4d %4d',event(in,1),gap(in,1));
        else
            fprintf(fid,'%4d %4d\n',event(in,1),gap(in,1));
        end
    end
    fclose(fid);
end