function varargout = get_initl_eq(outfilename,printflag)
%GET_INITL_EQ Read .OUT file from VELEST simultaneous run for final EQ locs
%   GET_INITL_EQ(OUTFILE,PRINTFLAG) Reads .OUT file for earthquake
%   locations from final iteration of VELEST simultaneous run and outputs a
%   locations_raw_in.dat file if PRINTFLAG is set to 1.
%
%   [LATITUDE, LONGITUDE, DEPTH] =
%   = GET_INITL_EQ(OUTFILENAME,PRINTFLAG) Reads .OUT file for earthquake
%   locations from final iteration of VELEST simultaneous run and outputs
%   LATITUDE, LONGITUDE, and DEPTH of thos final earthquake locations and
%   outputs locations_raw_in.dat file if PRINTFLAG is set to 1.
%
% Example:
%   [latitude,longitude,depth] =
%   get_initl_eq('fuegorandtrial01.OUT',1);

% Written and last edited by Kyle Brill 2016/04/08 13:45 EDT

%% Initialize variables.
if nargin<=1
    printflag = 0;
end

%% Read file into MATLAB
fout = fopen(outfilename,'r');
outfile = textscan(fout,'%s','delimiter','\n');
outfile = outfile{1,1};
fclose(fout);

%% Establish Beginning and End of Information
startRow=find(~cellfun(@isempty,strfind(outfile,'I N P U T - D A T A')))+3;
endRow=find(~cellfun(@isempty,strfind(outfile,'readings for station')))-2;
eqinfo = outfile(startRow(1,1):endRow(1,1),1);

%% Rearrange data into useful numbers
eqcolumns=cell(length(eqinfo),13);
for in=1:length(eqinfo)
    eqcolumns(in,:)=strsplit(eqinfo{in,1},' ');
end

eventnum = cellfun(@str2double,eqcolumns(:,1));
% date=datenum([cell2mat(eqcolumns(:,2)) cell2mat(eqcolumns(:,3)) num2str(str2double(eqcolumns(:,4)),'%05.2f\n')],'yymmddHHMMSS.FFF');

latitude=zeros(length(eqcolumns),1);
longitude=zeros(length(eqcolumns),1);
for in=1:length(eqcolumns)
    latitude(in,1)=str2double(eqcolumns{in,5}(1,1:7));
    if strcmp('S',eqcolumns{in,5}(1,8))
        latitude(in,1)=-latitude(in,1);
    end
    longitude(in,1)=str2double(eqcolumns{in,6}(1,1:7));
    if strcmp('W',eqcolumns{in,6}(1,8))
        longitude(in,1)=-longitude(in,1);
    end
end

depth=str2double(eqcolumns(:,7));
% mag=str2double(eqcolumns(:,11));
% no=str2double(eqcolumns(:,13));
% ifx=str2double(eqcolumns(:,12));
% x=str2double(eqcolumns(:,8));
% y=str2double(eqcolumns(:,9));
% z=str2double(eqcolumns(:,10));

%% Outputs

if nargout == 3
    varargout{1} = latitude;
    varargout{2} = longitude;
    varargout{3} = depth;
% elseif nargout == ??  %Add cases here for expanded outputs.

end

%% Print output for GMT plotting

if printflag == 1
    
    outlocfilename=strrep(outfilename,...
        outfilename(1,end-3:end),'locations_raw_out.dat');
    
    fid=fopen(outlocfilename,'w');
    fprintf(fid,'#Event Latitude Longitude Depth\n');
    for in=1:length(eventnum)
        if in==length(eventnum)
            fprintf(fid,'%4d %7.4f %7.4f  %5.2f',...
                eventnum(in,1),latitude(in,1),longitude(in,1),...
                depth(in,1));
        else
            fprintf(fid,'%4d %7.4f %7.4f  %5.2f\n',...
                eventnum(in,1),latitude(in,1),longitude(in,1),...
                depth(in,1));
        end
    end
    fclose(fid);
end