%load DEM for Guate Region, and save it as a .mat. 
clear,close all

%load the station coordinates, only for plotting purposes
filename = 'D:\Google Drive\Research\Fuego2012\2012stations.dat';
delimiter = ',';
startRow = 2;

formatSpec = '%s%f%f%f%[^\n\r]';

fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

fclose(fileID);
name = dataArray{:, 1};
stalat = dataArray{:, 2};
stalon = dataArray{:, 3};
elev = dataArray{:, 4};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%% Write station list out for use with velest
stalonout=stalon*-1;
stalonout=num2str(stalonout,'%#03.4f');
stalonout=[repmat('0',size(stalon)),stalonout,repmat('W',size(stalon))];
stalatout=num2str(stalat,'%02.4f');
stalatout=[stalatout,repmat('N',size(stalat))];
stanameout=repmat('____',size(name));
for in=1:length(name)
    if length(name{in,1})==1
        stanameout(in,:)=[name{in,1} '___'];
    elseif length(name{in,1})==2
        stanameout(in,:)=[name{in,1} '__'];
    elseif length(name{in,1})==3
        stanameout(in,:)=[name{in,1} '_'];
    elseif length(name{in,1})>3
        stanameout(in,:)=name{in,1}(1:4,1);
    end
end
staelevout=stanameout;
for in=1:length(elev)
    staelevout(in,:)=sprintf('%04d',elev(in,1)-2000);
end


inout=num2str(ones(10,1));
icc=(1:1:length(elev))';
for in=1:length(elev)
    iout(in,:)='1';
    iccout(in,:)=sprintf('%02d',icc(in,1));
    ipdel(in,:)='0.00';
    isdel(in,:)='0.00';
    spaces(in,:)=' ';
end

header='(a4,f7.4,a1,1x,f8.4,a1,1x,i4,1x,i1,1x,i2,1x,f5.2,2x,f5.2)';
body=[stanameout,spaces, stalatout, spaces, stalonout, spaces, staelevout,spaces,iout,spaces,iccout,spaces,ipdel,spaces,spaces,isdel];

dlmwrite('fuego0.1.sta',header,'delimiter','')
dlmwrite('fuego0.1.sta',body,'-append','delimiter','')
%%
% load GuateGeoTiff
[A R]=geotiffread('D:\Google Drive\GIS\Guate_Regional5.TIF');
pixeldeg=R.CellExtentInLatitude;

% clean up 
yy=linspace(R.LatitudeLimits(1,1),R.LatitudeLimits(1,2),R.RasterSize(1,1));
xx=linspace(R.LongitudeLimits(1,1),R.LongitudeLimits(1,2),R.RasterSize(1,2));
A=double(A);
A(A<-2000)=NaN;

figure(1)
imagesc(xx,fliplr(yy),A)
hold on
plot(stalon,stalat,'*k')
axis xy
axis equal
%%
xwin=12002-6900:12002-6400; 
ywin=4700:5200;

figure(2)
imagesc(A(6400:6900,4700:5200))
axis equal

latlist=R.LatitudeLimits(1,1):R.CellExtentInLatitude:R.LatitudeLimits(1,2);
lonlist=R.LongitudeLimits(1,1):R.CellExtentInLongitude:R.LongitudeLimits(1,2);
latlist=latlist(1,2:end);
lonlist=lonlist(1,2:end);
latlist=latlist(xwin);
latlist=fliplr(latlist);
lonlist=lonlist(ywin); 

% figure(3)
% contour(lonlist,latlist,A(6400:6900,4700:5200),100)
A=A(6400:6900,4700:5200);

%convert from decimal degrees to UTM
for i=1:length(latlist)
    lonlist2(i,1:length(lonlist))=lonlist;
end
for i=1:length(lonlist)
    latlist2(1:length(latlist),i)=latlist;
end
for i=1:length(latlist)
    [x(i,:) y(i,:) z]=deg2utm(latlist2(i,:),lonlist2(i,:));
end

figure(4)
contour(x,y,A,100)
axis equal
hold on

%load station coordinates and convert them also in UTM
[stax stay staz]=deg2utm(stalat,stalon);
%figure(5)
plot(stax,stay,'*K')

%save('Pacaya_resized_DEM.mat','stax','stay','x','y','A') 

%% save topo file for VELEST input

test=[latlist2(1,1) lonlist2(1,1) A(1,1)];

count=1;
for i=1:length(latlist)
    
    for j=1:length(lonlist)
        testtopo(count,:)=[latlist2(i,j) lonlist2(i,j) A(i,j)];
        count=count+1;
    end
end

fgtopo=[testtopo(:,2) testtopo(:,1) testtopo(:,3)];

%use a newdatum (at 2000 m) so you can use 2 velocity layes in the cone
%without violating velest rules. 

newdatum=testtopo(:,3)-2000;
fg1topo=[testtopo(:,2) testtopo(:,1) newdatum];
%%
tic
dlmwrite('fgtopo-2Kdatum.xyz',fg1topo,'delimiter',' ','precision','%.7f')
toc
%%
tic
dlmwrite('fgtopo-orig.xyz',fgtopo,'delimiter',' ','precision','%.7f')
toc

