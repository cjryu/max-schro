% For now, this only averages epsilon. If averaged mu is needed on nodes,
% for each point, all the mu edges (12 in total) have to be summed up.
ep = (cat(1,epEx,epEx(Nx,:,:)) + cat(1,epEx(1,:,:),epEx) + cat(2,epEy,epEy(:,Ny,:)) + cat(2,epEy(:,1,:),epEy) + cat(3,epEz,epEz(:,:,Nz)) + cat(3,epEz(:,:,1),epEz))/6;

% sigx = (cat(1,sigxEx,sigxEx(Nx,:,:)) + cat(1,sigxEx(1,:,:),sigxEx) + cat(2,sigxEy,sigxEy(:,Ny,:)) + cat(2,sigxEy(:,1,:),sigxEy) + cat(3,sigxEz,sigxEz(:,:,Nz)) ...
%        + cat(3,sigxEz(:,:,1),sigxEz))/6;
% sigy = (cat(1,sigyEx,sigyEx(Nx,:,:)) + cat(1,sigyEx(1,:,:),sigyEx) + cat(2,sigyEy,sigyEy(:,Ny,:)) + cat(2,sigyEy(:,1,:),sigyEy) + cat(3,sigyEz,sigyEz(:,:,Nz)) ...
%        + cat(3,sigyEz(:,:,1),sigyEz))/6;
% sigz = (cat(1,sigzEx,sigzEx(Nx,:,:)) + cat(1,sigzEx(1,:,:),sigzEx) + cat(2,sigzEy,sigzEy(:,Ny,:)) + cat(2,sigzEy(:,1,:),sigzEy) + cat(3,sigzEz,sigzEz(:,:,Nz)) ...
%        + cat(3,sigzEz(:,:,1),sigzEz))/6;