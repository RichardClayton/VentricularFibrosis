function MakePatchyScar_S1234_isthmus_2

% Author        R.H.Clayton (r.h.clayton@sheffield.ac.uk)
% 
% This file is part of VentricularFibrosis.
%
% Copyright (c) Richard Clayton, 
% Department of Computer Science, 
% University of Sheffield, 2006, 2009, 2016, 2023
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

% creates directory structure for 20 samples of 'scar'
% at each of four length scales

% this version includes a central isthmus for figure of 8 re-entry

% this code makes use of the circulant embedding code included in 
% stationary_Gaussian_process.m, available from 
% https://uk.mathworks.com/matlabcentral/fileexchange/38880...
% -circulant-embedding-method-for-generating-stationary-gaussian-field/...
% content/stationary_Gaussian_process.m
%
% Reference:
% Kroese, D. P., & Botev, Z. I. (2015). Spatial Process Simulation.
% In Stochastic Geometry, Spatial Statistics and Random Fields(pp. 369-404)
% Springer International Publishing, DOI: 10.1007/978-3-319-10064-7_12

    maxD = 0.1;
    minD = 0.0; % 0.025;
    
    lengthScale = [5, 10, 20, 40];
    nSamples = 20;
    
    for l = 1:length(lengthScale)
        for sample = 1:nSamples
            newDir = sprintf('lengthScale%d_isthmus_2/sample%d',lengthScale(l),sample);
            if ~exist(newDir, 'dir')
                system(sprintf('mkdir %s',newDir));
            end
            PatchyDiffusion(maxD, minD, sample, lengthScale(l));
        end
        close all
    end

end

function PatchyDiffusion(maxD, minD, sample, lengthScale)
        
% generates diffusion field for simulating scar and border zone

    nRows = 400;
    nCols = 400;
        
    % set up covariance function
    rho=@(h)(exp(-((h(1)^2 + h(2)^2)/lengthScale^2)));

    % now obtain GRF with mean 0 and SD 1 -- so ranges between -3 and +3
    [f1,f2,tx,ty] = stationary_Gaussian_process(nRows,nCols,rho);
  
    D1 = calculateD(maxD, minD, f1, nRows, nCols);
    filePath = sprintf('lengthScale%d_isthmus/sample%d/',lengthScale,sample);
    writeToFile(filePath,D1,nRows,nCols);
    
end

function D = calculateD(maxD, minD, f1, nRows, nCols) 
    % apply sigmoid function to GRF
    rInner = 50;                % radius of scar core
    rOuter = 120;               % radius of border zone
    slope = -0.075;             % slope of sigmoid function

    
    for row = 1:nRows
        for col = 1:nCols           
            DGRF(row,col) = maxD * (f1(row,col) + 2.0)/4.0;

            % infarct field 1           
            x1 = col - nCols/3;
            y1 = row - 2*nRows/3;
            % y1 = row - nRows/2;
            r1(row,col) = sqrt(x1*x1 + y1*y1);
            infarct1(row,col) = 1.0/(1.0 + exp(slope * (r1(row,col) - rInner)));
            border1(row,col) = 1.0/(1.0 + exp(slope * (r1(row,col) - rOuter)));
            
            % infarct field 2 
            x2 = col - 2*nCols/3;
            % y2 = row - nRows/2;
            y2 = row - nRows/3;
            r2(row,col) = sqrt(x2*x2 + y2*y2);
            infarct2(row,col) = 1.0/(1.0 + exp(slope * (r2(row,col) - rInner)));            
            border2(row,col) = 1.0/(1.0 + exp(slope * (r2(row,col) - rOuter)));

            D(row,col) = maxD * (border1(row,col) * border2(row,col)) + (infarct1(row,col) * infarct2(row,col)) * DGRF(row,col);
                        
            if (D(row,col) < minD)
                D(row,col) = minD-0.001;
            end
            if (D(row,col) > maxD)
                D(row,col) = maxD;
            end           
        end
    end
 
end

function writeToFile(filePath, D, nRows, nCols)

        % display
        figure
        pcolor(D); shading flat
        colorbar;
        colormap(hot);
        
        % and write data and image to file
        fid = fopen(sprintf('%s/DiffusionCoefficient.txt',filePath),'w');
        for i = 1:nRows
            for j = 1:nCols
                fprintf(fid,'%5.3f ',D(i,j));
            end
            fprintf(fid,'\n');
        end
        fclose(fid);

        plotFile = sprintf('%s/GRF.png',filePath);
        print(gcf,'-dpng',plotFile);
        
        % generate region of random removal
        fid = fopen(sprintf('%s/DiffusionCoefficient_random.txt',filePath),'w');
        for i = 1:nRows
            for j = 1:nCols
                if (rand >(20.0*D(i,j)))
                    D(i,j) = -1;
                end
                fprintf(fid,'%5.3f ',D(i,j));
            end
            fprintf(fid,'\n');
        end
        fclose(fid);
        
        figure
        pcolor(D); shading flat
        colorbar;
        colormap(hot);
        caxis([0 0.1]);
        
        plotFile = sprintf('%s/GRF_random.png',filePath);
        print(gcf,'-dpng',plotFile);
        
        
end