function out = getSequimPointData(xq,yq)
% xq, yq is the query coordinates

% Code to get vels at a certain point
% determin 6 nearest elements to given location, load in their files

% calculate z vel bins from seafloor to minimum surface height

% iterate through time. At each time, interpolate vels onto z vel
% bins. 

N = 5; % Number of surrounding elements to use in interpolation

% load in element positions, water level timeseries, timestamps, and
% bathymetry
mds = load('MetaData.mat');

% calculate 10 z query locations, spaced evenly from the bottom to the minimum
% water level
seabed = mds.bathySIO(xq,yq);
surface = mds.WLmin;
H = surface - seabed;
zq = fliplr((H/20):(H/10):(H-H/20))';

% calculate element distances from query point
d = sqrt((xq-mds.x).^2+(yq-mds.y).^2);

[mind,inds] = sort(d,'ascend');

if mind(1) == 0
    % then the query point is exactly on a node, do not need interpolation
    % in x and y
    
else
    inds = inds(1:N); % indicies of 6 nearest elements to use for interpolation
    
    for i = 1:N % load in surrounding points
        vEl(i) = load(sprintf('VelFiles/loc_%04.0f',inds(i)));
        vEl(i).u = double(vEl(i).u);
        vEl(i).v = double(vEl(i).v);
        vEl(i).w = double(vEl(i).w);
        vEl(i).z = vEl(i).z';
    end

    
    xt = []; yt = [];
    for i = 1:N
        xt = [xt ; ones(10,1).*vEl(i).x];
        yt = [yt ; ones(10,1).*vEl(i).y];
    end
    xq = repmat(xq,10,1);
    yq = repmat(yq,10,1);
    
        
   
     ut = zeros(10*N,1); vt = zeros(10*N,1); wt = zeros(10*N,1); zt = zeros(10*N,1);

     for j = 1:length(mds.t) % iterate through time
       
        
        for i = 1:N
            ut((i-1)*10+1:i*10) = vEl(i).u(:,j);
            vt((i-1)*10+1:i*10) = vEl(i).v(:,j);
            wt((i-1)*10+1:i*10) = vEl(i).w(:,j);
            zt((i-1)*10+1:i*10) = vEl(i).z(:,j);
        end

   
        SIO = scatteredInterpolant(xt,yt,zt,ut);
        out.u(:,j) = SIO(xq,yq,zq);
        SIO.Values = vt;
        out.v(:,j) = SIO(xq,yq,zq);
        SIO.Values = wt;
        out.w(:,j) = SIO(xq,yq,zq);
    end
end

out.zq = zq;
out.seabed = seabed;
out.minSurface = surface;
out.t = mds.t;