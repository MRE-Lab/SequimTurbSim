% create location specific files
% time series for u,v,w at all z for each element position
% time series of water height
% water depth

for i = 2346:4461
    i
    out.u = squeeze(u(:,i,:));
    out.v = squeeze(v(:,i,:));
    out.w = squeeze(w(:,i,:));
    out.x = element.x(i);
    out.y = element.y(i);
    out.d = element.d(i);
    [out.z,out.seabed,out.surface] = getHeightInfo(out.x,out.y,Tmatlab,WLm,bathySIO);
    save(sprintf('VelFiles/loc_%04.0f',i),'-struct','out')
end
    





