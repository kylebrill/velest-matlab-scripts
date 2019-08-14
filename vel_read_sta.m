function varargout = vel_read_sta(stafile,printflag,startRow,endRow)
%VEL_READ_STA Import contents of .STA output from VELEST
%   VEL_READ_STA(STAFILE, PRINTFLAG) Reads data from .STA file STAFILE and
%   outputs a .DAT file for use with GMT plotting scripts if PRINTFLAG is
%   set equal to 1.
%
% function [station,ptcor,latitude,longitude,elevation] = vel_read_sta(stafile, startRow, endRow)
%
%   [STATION,PTCOR]
%   = VEL_READ_STA(STAFILE,PRINTFLAG, STARTROW, ENDROW) Reads data from 
%   rows STARTROW through ENDROW of .sta file STAFILE and outputs STATION 
%   names as a Nx1 cell array and PTCOR P-wave arrival station corrections 
%   as a Nx1 numerical array.
%
% Example:
%   [station,ptcor]
%   = vel_read_sta('velest.sta', 1, 2, 11);
%
%   [STATION,PTCOR,LATITUDE,LONGITUDE,ELEVATION]
%   = VEL_READ_STA(STAFILE,PRINTFLAG, STARTROW, ENDROW) Reads data from 
%   rows STARTROW through ENDROW of .sta file STAFILE and outputs STATION 
%   names as a Nx1 cell array, PTCOR P-wave arrival station corrections as
%   a Nx1 numerical array, and station LATITUDE, LONGITUDE, and ELEVATION
%   values as Nx1 numerical arrays.
%
% Example:
%   [station,ptcor,latitude,longitude,elevation]
%   = vel_read_sta('velest.sta', 1, 2, 11);
%
%    See also TEXTSCAN.
%    Function can be expanded to output more information if need arises.

% Auto-generated by MATLAB on 2016/04/04 10:35:53
% Expanded and Last Edited by Kyle Brill on 2016/04/04 18:00:00 EDT

%% Initialize variables.
if ~any([nargout==0,nargout==2,nargout==5])
    error('Number of output variables should be none, 2, or 5.')
end
if nargin<=3
    startRow = 2;
    endRow = inf;
end
if nargin==2
    if printflag>1
        error('1 is the only acceptable Print Flag. Check inputs.');
    end
end
if nargin<=1
    printflag = 0;
end

%% Format string for each line of text:
%   column1: text (%s)
%	column2: double (%f)
%   column3: text (%s)
%	column4: double (%f)
%   column5: text (%s)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%4s%7f%1s%9f%1s%5f%2f%3f%6s%7f%[^\n\r]';

%% Open the text file.
fileID = fopen(stafile,'r');
finfo = dir(stafile);
if isempty(finfo)
    error('Unable to find .sta file as entered.')
end

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
station = dataArray{:, 1};
latitude = dataArray{:, 2};
latNS = dataArray{:, 3};
longitude = dataArray{:, 4};
lonEW = dataArray{:, 5};
elevation = dataArray{:, 6};
model = dataArray{:, 7};
icc = dataArray{:, 8};
ptcor = dataArray{:, 9};
stcor = dataArray{:, 10};


%% Adjust outputs
for ind = 1:length(latitude)
    if strcmp(latNS{ind},'S')
        latitude(ind)=-latitude(ind);
    end
    if strcmp(lonEW{ind},'W')
        longitude(ind)=-longitude(ind);
    end
end

% This accounts for station corrections greater than 9.99
ptcor=str2double(ptcor);

%% Outputs

if nargout == 2
    varargout{1} = station;
    varargout{2} = ptcor;
elseif nargout == 5
    varargout{1} = station;
    varargout{2} = ptcor;
    varargout{3} = latitude;
    varargout{4} = longitude;
    varargout{5} = elevation;
% elseif nargout == ??  %Add cases here for expanded outputs.
end


%% Print output for GMT plotting

if printflag == 1
    stadatfilename=strrep(stafile,stafile(1,end-3:end),'sta.dat');
    fid=fopen(stadatfilename,'w');
    for in=1:length(icc)
        if in==length(icc)
            fprintf(fid,'%s %7.4f %7.4f  %4d %d %2d %5.2f  %5.2f',...
                station{in,1},latitude(in,1),longitude(in,1),...
                elevation(in,1),model(in,1),icc(in,1),...
                ptcor(in,1),stcor(in,1));
        else
            fprintf(fid,'%s %7.4f %7.4f  %4d %d %2d %5.2f  %5.2f\n',...
                station{in,1},latitude(in,1),longitude(in,1),...
                elevation(in,1),model(in,1),icc(in,1),...
                ptcor(in,1),stcor(in,1));
        end
    end
    fclose(fid);
end