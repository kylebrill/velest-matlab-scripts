function xsec_gmt(grdfile)
% XSEC_GMT opens DEM created for GMT and determine cross section lines.
%
%   XSEC_GMT (GRDFILE) allows user to graphically select cross-section
%   lines and then prints lines to be used in bash script for printing
%   event locations from VELEST
%
% Example:
%   xsec_gmt('/local/kabrill_res/velest_template/gmtdata/fg_grid.grd');
%
% Requires: grdread2.m and grdinfo2.m by Kelsey Jordahl, available on 
%           Matlab Exchange
%
% Written by Kyle Brill
% Last updated 3/30/2016 at 20:00 EDT

%% Use Kelsey Jordahl functions to load .grd file and extract information

[X,Y,Z]=grdread2(grdfile);
D=grdinfo2(grdfile);
%D=(xmin, xmax, ymin, ymax, zmin, zmax, format, xinc, yinc). Format is
%  1 for pixel registration and 0 for grid node registration.
xmin=D(1);xmax=D(2);ymin=D(3);ymax=D(4);xinc=D(8);yinc=D(9);

%% Plot .grd file and select area where cross section will be drawn through

FIG=figure(7);
imagesc(X,Y,Z);
axis xy

disp('Zoom in on target area and press space when view is ready.')
pause;
disp('Click on the point you wish cross-sections to meet.')
pause(3);
figure(7)
[xcor,ycor]=ginput(1);

%% Output info needed to plot velest results with cross-section profiles

fprintf('EastWest line center (for if line) =\n%16.13f\n',ycor);
fprintf('E-W line =\n%8.4f/%8.4f\n\n',ycor-yinc, ycor+yinc);
fprintf('NorthSouth line center (for if line) =\n%16.13f\n',xcor);
fprintf('N-S line =\n%8.4f/%8.4f\n\n',xcor-xinc, xcor+xinc);