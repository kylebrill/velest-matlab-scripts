clear,close all

%load digital elevation map, resized from default
% addpath('~/Documents/Pacaya_2013/make_topo_vel/')

% %define the corners of the area of interest
% %Upper left coordinates
% lonul=-90.75004575;
% latul=14.50004746;
% %Lower left coordinates
% lonll=-90.75004575;
% latll=14.16662146;
% %upper right coordinates
% lonur=-90.49997625;
% latur=14.50004746;
% %lower right coordinates
% lonlr=-90.49997625;
% latlr=14.16662146;
% %center, not useful?
% % lonc=-72.4999999;
% % latc=-37.4999996;

% load GuateGeoTiff
[A R]=geotiffread('/mnt/data/Google Drive/GIS/Guate_Regional5.TIF');
pixeldeg=R.CellExtentInLatitude;

% %WGS84
% %Model Pixel scale tag:  9.1500e-05
% pixeldeg=9.1500e-05;
% [A R]=geotiffread('Elevation_2058_I_and_2059_II_IGN_2006_Pacaya_10m_v1.tif');
% %pause

yy=linspace(R.LatitudeLimits(1,1),R.LatitudeLimits(1,2),R.RasterSize(1,1));
xx=linspace(R.LongitudeLimits(1,1),R.LongitudeLimits(1,2),R.RasterSize(1,2));
A=double(A);
A(A<-2000)=NaN;

% yy=linspace(14.166621464,14.50004746,2733);
% xx=linspace(-90.75004575,-90.49997625,3644);

% figure(1)
% imagesc(xx,fliplr(yy),A)
% axis xy
% axis equal

xwin=12002-6900:12002-6400; 
ywin=4700:5200;
% xwin=700:1950; 
% ywin=950:2420;

%3400-3900 %5000-6000
%imagesc(A(5200:5450,3600:3800))
% figure(2)
% imagesc(A(6400:6900,4700:5200))
% axis equal
% figure(2)
% imagesc(A(xwin,ywin))
% axis equal

A=A(6400:6900,4700:5200);
% A=A(xwin,ywin);

latlist=R.LatitudeLimits(1,1):R.CellExtentInLatitude:R.LatitudeLimits(1,2);
lonlist=R.LongitudeLimits(1,1):R.CellExtentInLongitude:R.LongitudeLimits(1,2);
latlist=latlist(1,2:end);
lonlist=lonlist(1,2:end);
latlist=latlist(xwin);
latlist=fliplr(latlist);
lonlist=lonlist(ywin); 

% latlist=latul:-pixeldeg:latll;
% lonlist=lonll:pixeldeg:lonlr;
% latlist=latlist(xwin);
% lonlist=lonlist(ywin); 

figure(3)
contour(lonlist,latlist,A,50)

% get corners of small bounding box for velest to attempt first 
boxed=ginput(2);

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

% figure(4)
% contour(x,y,A,100)
% axis equal
% hold on

%load station coordinates and convert them also in UTM
% [stax stay staz]=deg2utm(stalat,stalon);
%figure(5)
% plot(stax,stay,'*K')

% 
% %convert from decimal degrees to UTM
% for i=1:length(latlist)
%     lonlist2(i,1:length(lonlist))=lonlist;
% end
% for i=1:length(lonlist)
%     latlist2(1:length(latlist),i)=latlist;
% end
% for i=1:length(latlist)
%     [x(i,:) y(i,:) z]=deg2utm(latlist2(i,:),lonlist2(i,:));
% end
% 
% figure(4)
% contour(x,y,A,100)
% axis equal
% hold on

% save ('Fuego_DEM.mat','x','y','A')
% 
% %save('Pacaya_resized_DEM_2.mat','x','y','A') 

%% save topo file for VELEST input

count=1;
for i=1:length(latlist)
    
    for j=1:length(lonlist)
        testtopo(count,:)=[lonlist2(i,j) latlist2(i,j) A(i,j)];
        count=count+1;
    end
end

% fgtopo=[testtopo(:,1) testtopo(:,2) testtopo(:,3)];

%use a newdatum (at 2000 m) so you can use 2 velocity layes in the cone
%without violating velest rules. 

newdatum=testtopo(:,3)-2000;
fgtopo=[testtopo(:,1) testtopo(:,2) newdatum];


% dlmwrite('fgtopo.xyz',fg1topo,'delimiter',' ','precision','%.7f')
% dlmwrite('pac1topo-orig.xyz',ftopo,'delimiter',' ','precision','%.7f')

%rearrange topo file lines so velest can work faster
mvtotop=find(testtopo(:,1) > boxed(1,1) & testtopo(:,1) < boxed(2,1) & testtopo(:,2) < boxed(1,2) & testtopo(:,2) > boxed(2,2));
lookfirst=testtopo(mvtotop,1);
lookfirst(:,2)=testtopo(mvtotop,2);
lookfirst(:,3)=newdatum(mvtotop,1);

mX=unique(lookfirst(:,1));
mY=unique(lookfirst(:,2)); mY=flipud(mY);
mZ=reshape(lookfirst(:,3),length(mX),length(mY))';
contour(mX,mY,mZ)

addtobot=[testtopo(:,1:2),newdatum];
addtobot(mvtotop,:)=[];

speedupvel=[lookfirst;addtobot];

% dlmwrite('fgtopo-2Kfast.xyz',speedupvel,'delimiter',' ','precision','%.7f')
%% make grd for gmt plotting
% %load DEM for Pacaya, cut it to the area of interest and save it as a .mat. 
% load('Pacaya_resized_DEM_2.mat')
% 
newdem=double(A);
% 
clear A
modelspacing=90;

%parameters set to make the grid a box, instead of an inclined box
modx=floor(x(1,1)):modelspacing:ceil(x(1,501)+90);
mody=floor(y(501,1)):modelspacing:ceil(y(1,1))-90*13;
% modx=752600:modelspacing:765000;
% mody=1584600:modelspacing:1597000;
% 

topomodel=griddata(x,y,newdem,modx,mody');
% 
figure(1)
contour(topomodel)
% 
% zone='15 P';
% utmzone=zeros(1,length(modx));
% 
for i=1:length(modx)
utmzone(i,:)='15 P';
end
% 
[xx,yy]=utm2deg(modx,mody,utmzone);

grdwrite2(yy,xx,topomodel,'fg_grid.grd');