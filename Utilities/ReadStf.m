function StfData=ReadStf(dims,fname)

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

% ReadStf : Reads and returns data from
%   file <fname> that is in simple text
%   (stf) format. For 2D data (dims = 2) 
%   the returned variable will be a 2D 
%   matrix, and for 3D data (dims = 3)
%   a 3D matrix.

% read in data

if (dims < 2) || (dims > 3)
    disp('ReadStf error: dimension must be either 2 or 3');
else
    % check that fname exists
    if exist(fname) > 0
        stf_fid = fopen( fname,'rt' );
        gzipFlag = 0;
    % otherwise try opening gzipped version
    else
        fnamegz = strcat(fname, '.gz');
        % eval(['!copy ' fnamegz ' t.gz']);
        gunzip(fnamegz);
        %!gunzip t.gz
        %stf_fid = fopen( 't', 'rt' );
        stf_fid = fopen( fname, 'rt');
        gzipFlag = 1;
    end

% read header information
    line_1 = fgetl(stf_fid);
    line_2 = fgetl(stf_fid);
    line_3 = fgetl(stf_fid);
    if dims == 2
        Dimensions = sscanf(line_3, 'DIMENSIONS %i %i');
    else 
        Dimensions = sscanf(line_3, 'DIMENSIONS %i %i %i');
    end
    line_4 = fgetl(stf_fid);
    if dims == 2 
        Bounds = sscanf(line_4, 'BOUNDS %i %i %i %i');
    else 
        Bounds = sscanf(line_4, 'BOUNDS %i %i %i %i %i %i');
    end
    line_5 = fgetl(stf_fid);
    line_6 = fgetl(stf_fid);

    if dims == 2
        StfData = zeros(Dimensions(2),Dimensions(1));
    else
        StfData = zeros(Dimensions(2),Dimensions(1),Dimensions(3));
    end
    
    % small piece of code to offset data if bounds are < 0
    Offset = [0 0 0 0 0 0];
    for n = 1:2:2*dims
        if (Bounds(n) < 0)
            Offset(n) = Bounds(n)*-1;
            Offset(n+1) = Bounds(n)*-1;
        end
    end
                
% now read data
    if dims == 2
        %h = waitbar(0,'Reading data...');
        for row = Bounds(3)+Offset(3)+1:Bounds(4)+Offset(4)+1
            %waitbar(row/(Bounds(4)+Offset(4)+1-Bounds(3)+Offset(3)+1),h);
            data_row = fgetl(stf_fid);
            TmpVector = strread(data_row);
            StfData(row,:) = TmpVector(1:Dimensions(1));
        end
        %close(h);
    else 
        h = waitbar(0,'Reading data...');
        for layer = Bounds(5)+Offset(5)+1:Bounds(6)+Offset(6)
            waitbar(layer/(Bounds(6)+Offset(6)-Bounds(5)+1+Offset(5)),h);
            for row = Bounds(3)+Offset(3)+1:Bounds(4)+Offset(4)
                data_row = fgetl(stf_fid);
                if (length(data_row) == 0)
                    data_row = fgetl(stf_fid);
                end
                TmpVector = strread(data_row);
                StfData(row,1:Dimensions(1),layer) = TmpVector(1:Dimensions(1));
            end
            % this line should be commented out for lr data
            data_row = fgetl(stf_fid);
        end
        close(h);
    end

    fclose( stf_fid );
    if (gzipFlag == 1)
        %delete ('t')
        gzip(fname);
        rmFile = sprintf('rm %s',fname);
        system(rmFile);
    end
end
