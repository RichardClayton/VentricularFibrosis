function MakeImages2D(froot,start,stop,step)

% Author        R.H.Clayton (r.h.clayton@sheffield.ac.uk)
% 
% This file is part of VentricularFibrosis.
%
% Copyright (c) Richard Clayton, 
% Department of Computer Science, 
% University of Sheffield, 2003
%
% VentricularFibrosis is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

% MakeImages2D : Makes series of images from 2D stf files
%   with filenames <froot>XXXX.stf

% initialise figure and movie parameters

fig=figure;
set(fig,'DoubleBuffer','on');
set(gca,'NextPlot','replace');

for FileNum = start:step:stop
    
% put together current filename
    fname = sprintf(strcat(froot,'%04d.stf'),FileNum);
    disp(fname);

% read data from file, and plot to fig window
    % PlotStf2D(fname);
    DataToPlot = ReadStf(2,fname);
    [ny,nx] = size(DataToPlot);

    % remove boundary layer
    %DataToPlot = Data(2:ny-1,2:nx-1);

    % set up 2D view
    pcolor(DataToPlot);
    colormap(hot);
    shading interp;

    % tinker with axes
    axis equal; axis tight;
    set(gcf,'Color','white');
    set(gca,'Visible','off');
    set(gca,'XTickLabel',[]); set(gca,'XTick',[]);
    set(gca,'YTickLabel',[]); set(gca,'YTick',[]);
    caxis([-90,30]);  % TP06
%     caxis([-80,10]);  % ? CRN
    % PlotStf2Dfloat(fname,-90,50);
    % PlotStf2Dfloat(fname,100,550);

% add time to figure, assuming data files at 1 ms intervals
%     TimeString = sprintf('%4d ms',FileNum);
%     text(10,10,TimeString,'Color','blue','FontWeight','bold')
    
% print to file
    plotfilename = sprintf(strcat(froot,'%04d.jpg'),FileNum);
    print(gcf, '-djpeg', plotfilename);
    
    % close all
    
end