% PML requires ep, mu, and sig on H grid. At this point, we only have ep on
% E grid and mu on H grid. (We don't need mu on E grid.)

epEx = epx; epEy = epy; epEz = epz;
clear epx epy epz
epHx = (ep(1:Nx+1,1:Ny,1:Nz)+ep(1:Nx+1,2:Ny+1,1:Nz)+ep(1:Nx+1,1:Ny,2:Nz+1)+ep(1:Nx+1,2:Ny+1,2:Nz+1))/4;
epHy = (ep(1:Nx,1:Ny+1,1:Nz)+ep(2:Nx+1,1:Ny+1,1:Nz)+ep(1:Nx,1:Ny+1,2:Nz+1)+ep(2:Nx+1,1:Ny+1,2:Nz+1))/4;
epHz = (ep(1:Nx,1:Ny,1:Nz+1)+ep(2:Nx+1,1:Ny,1:Nz+1)+ep(1:Nx,2:Ny+1,1:Nz+1)+ep(2:Nx+1,2:Ny+1,1:Nz+1))/4;

sigxHx = zeros(Nx+1,Ny,Nz);     sigxHy = zeros(Nx,Ny+1,Nz);     sigxHz = zeros(Nx,Ny,Nz+1);
sigyHx = zeros(Nx+1,Ny,Nz);     sigyHy = zeros(Nx,Ny+1,Nz);     sigyHz = zeros(Nx,Ny,Nz+1);
sigzHx = zeros(Nx+1,Ny,Nz);     sigzHy = zeros(Nx,Ny+1,Nz);     sigzHz = zeros(Nx,Ny,Nz+1);
sigxEx = zeros(Nx,Ny+1,Nz+1);   sigxEy = zeros(Nx+1,Ny,Nz+1);   sigxEz = zeros(Nx+1,Ny+1,Nz);
sigyEx = zeros(Nx,Ny+1,Nz+1);   sigyEy = zeros(Nx+1,Ny,Nz+1);   sigyEz = zeros(Nx+1,Ny+1,Nz);
sigzEx = zeros(Nx,Ny+1,Nz+1);   sigzEy = zeros(Nx+1,Ny,Nz+1);   sigzEz = zeros(Nx+1,Ny+1,Nz);
sigx = zeros(Nx+1,Ny+1,Nz+1);   sigy = zeros(Nx+1,Ny+1,Nz+1);   sigz = zeros(Nx+1,Ny+1,Nz+1);   % For the PML nodes.

m = 2;
R = .000001;
Lx = dx*Np;
Ly = dy*Np;
Lz = dz*Np;
smaxx = -(m+1)/2/n0/Lx*log(abs(R));
smaxy = -(m+1)/2/n0/Ly*log(abs(R));
smaxz = -(m+1)/2/n0/Lz*log(abs(R));

% Bottom
sighx = ones(Ny,1)*smaxz*((dz*fliplr(0:Np-1)+dz/2)/Lz).^m;
sighy = ones(Ny+1,1)*smaxz*((dz*fliplr(0:Np-1)+dz/2)/Lz).^m;
sighz = ones(Ny,1)*smaxz*(dz*fliplr(0:Np)/Lz).^m;
sigex = ones(Ny+1,1)*smaxz*(dz*fliplr(0:Np)/Lz).^m;
sigey = ones(Ny,1)*smaxz*(dz*fliplr(0:Np)/Lz).^m;
sigez = ones(Ny+1,1)*smaxz*((dz*fliplr(0:Np-1)+dz/2)/Lz).^m;
for i = 1:Nx
    sigzHx(i,:,1:Np) = sighx;
    sigzHy(i,:,1:Np) = sighy;
    sigzHz(i,:,1:Np+1) = sighz;
    sigzEx(i,:,1:Np+1) = sigex;
    sigz(i,:,1:Np+1) = sigex;
    sigzEy(i,:,1:Np+1) = sigey;
    sigzEz(i,:,1:Np) = sigez;
end
sigzHx(Nx+1,:,1:Np) = sighx;
sigz(Nx+1,:,1:Np+1) = sigex;
sigzEy(Nx+1,:,1:Np+1) = sigey;
sigzEz(Nx+1,:,1:Np) = sigez;

% Top
sighx = ones(Ny,1)*smaxz*((dz*(0:Np-1)+dz/2)/Lz).^m;
sighy = ones(Ny+1,1)*smaxz*((dz*(0:Np-1)+dz/2)/Lz).^m;
sighz = ones(Ny,1)*smaxz*(dz*(0:Np)/Lz).^m;
sigex = ones(Ny+1,1)*smaxz*(dz*(0:Np)/Lz).^m;
sigey = ones(Ny,1)*smaxz*(dz*(0:Np)/Lz).^m;
sigez = ones(Ny+1,1)*smaxz*((dz*(0:Np-1)+dz/2)/Lz).^m;
for i = 1:Nx
    sigzHx(i,:,Nz-Np+1:Nz) = sighx;
    sigzHy(i,:,Nz-Np+1:Nz) = sighy;
    sigzHz(i,:,Nz-Np+1:Nz+1) = sighz;
    sigzEx(i,:,Nz-Np+1:Nz+1) = sigex;
    sigz(i,:,Nz-Np+1:Nz+1) = sigex;
    sigzEy(i,:,Nz-Np+1:Nz+1) = sigey;
    sigzEz(i,:,Nz-Np+1:Nz) = sigez;
end
sigzHx(Nx+1,:,Nz-Np+1:Nz) = sighx;
sigz(Nx+1,:,Nz-Np+1:Nz+1) = sigex;
sigzEy(Nx+1,:,Nz-Np+1:Nz+1) = sigey;
sigzEz(Nx+1,:,Nz-Np+1:Nz) = sigez;

% Left
sighx = transpose(smaxy*((dy*fliplr(0:Np-1)+dy/2)/Ly).^m)*ones(1,Nz);
sighy = transpose(smaxy*(dy*fliplr(0:Np)/Ly).^m)*ones(1,Nz);
sighz = transpose(smaxy*((dy*fliplr(0:Np-1)+dy/2)/Ly).^m)*ones(1,Nz+1);
sigex = transpose(smaxy*(dy*fliplr(0:Np)/Ly).^m)*ones(1,Nz+1);
sigey = transpose(smaxy*((dy*fliplr(0:Np-1)+dy/2)/Ly).^m)*ones(1,Nz+1);
sigez = transpose(smaxy*(dy*fliplr(0:Np)/Ly).^m)*ones(1,Nz);
for i = 1:Nx
    sigyHx(i,1:Np,:) = sighx;
    sigyHy(i,1:Np+1,:) = sighy;
    sigyHz(i,1:Np,:) = sighz;
    sigyEx(i,1:Np+1,:) = sigex;
    sigy(i,1:Np+1,:) = sigex;
    sigyEy(i,1:Np,:) = sigey;
    sigyEz(i,1:Np+1,:) = sigez;
end
sigyHx(Nx+1,1:Np,:) = sighx;
sigy(Nx+1,1:Np+1,:) = sigex;
sigyEy(Nx+1,1:Np,:) = sigey;
sigyEz(Nx+1,1:Np+1,:) = sigez;

% Right
sighx = transpose(smaxy*((dy*(0:Np-1)+dy/2)/Ly).^m)*ones(1,Nz);          % Nz,Ny,Nx+1
sighy = transpose(smaxy*(dy*(0:Np)/Ly).^m)*ones(1,Nz);            % Nz,Ny+1,Nx
sighz = transpose(smaxy*((dy*(0:Np-1)+dy/2)/Ly).^m)*ones(1,Nz+1);           % Nz+1,Ny,Nx
sigex = transpose(smaxy*(dy*(0:Np)/Ly).^m)*ones(1,Nz+1);            % Nz+1,Ny+1,Nx
sigey = transpose(smaxy*((dy*(0:Np-1)+dy/2)/Ly).^m)*ones(1,Nz+1);            % Nz+1,Ny,Nx+1
sigez = transpose(smaxy*(dy*(0:Np)/Ly).^m)*ones(1,Nz);            % Nz,Ny+1,Nx+1
for i = 1:Nx
    sigyHx(i,Ny-Np+1:Ny,:) = sighx;
    sigyHy(i,Ny-Np+1:Ny+1,:) = sighy;
    sigyHz(i,Ny-Np+1:Ny,:) = sighz;
    sigyEx(i,Ny-Np+1:Ny+1,:) = sigex;
    sigy(i,Ny-Np+1:Ny+1,:) = sigex;
    sigyEy(i,Ny-Np+1:Ny,:) = sigey;
    sigyEz(i,Ny-Np+1:Ny+1,:) = sigez;
end
sigyHx(Nx+1,Ny-Np+1:Ny,:) = sighx;
sigy(Nx+1,Ny-Np+1:Ny+1,:) = sigex;
sigyEy(Nx+1,Ny-Np+1:Ny,:) = sigey;
sigyEz(Nx+1,Ny-Np+1:Ny+1,:) = sigez;

% Back
sighx = transpose(smaxx*(dx*fliplr(0:Np)/Lx).^m)*ones(1,Nz);
sighy = transpose(smaxx*((dx*fliplr(0:Np-1)+dx/2)/Lx).^m)*ones(1,Nz);
sighz = transpose(smaxx*((dx*fliplr(0:Np-1)+dx/2)/Lx).^m)*ones(1,Nz+1);
sigex = transpose(smaxx*((dx*fliplr(0:Np-1)+dx/2)/Lx).^m)*ones(1,Nz+1);
sigey = transpose(smaxx*(dx*fliplr(0:Np)/Lx).^m)*ones(1,Nz+1);
sigez = transpose(smaxx*(dx*fliplr(0:Np)/Lx).^m)*ones(1,Nz);
for i = 1:Ny
    sigxHx(1:Np+1,i,:) = sighx;
    sigxHy(1:Np,i,:) = sighy;
    sigxHz(1:Np,i,:) = sighz;
    sigxEx(1:Np,i,:) = sigex;
    sigxEy(1:Np+1,i,:) = sigey;
    sigx(1:Np+1,i,:) = sigey;
    sigxEz(1:Np+1,i,:) = sigez;
end
sigxHy(1:Np,Ny+1,:) = sighy;
sigxEx(1:Np,Ny+1,:) = sigex;
sigx(1:Np+1,Ny+1,:) = sigey;
sigxEz(1:Np+1,Ny+1,:) = sigez;

% Front
sighx = transpose(smaxx*(dx*(0:Np)/Lx).^m)*ones(1,Nz);
sighy = transpose(smaxx*((dx*(0:Np-1)+dx/2)/Lx).^m)*ones(1,Nz);
sighz = transpose(smaxx*((dx*(0:Np-1)+dx/2)/Lx).^m)*ones(1,Nz+1);
sigex = transpose(smaxx*((dx*(0:Np-1)+dx/2)/Lx).^m)*ones(1,Nz+1);
sigey = transpose(smaxx*(dx*(0:Np)/Lx).^m)*ones(1,Nz+1);
sigez = transpose(smaxx*(dx*(0:Np)/Lx).^m)*ones(1,Nz);
for i = 1:Ny
    sigxHx(Nx-Np+1:Nx+1,i,:) = sighx;
    sigxHy(Nx-Np+1:Nx,i,:) = sighy;
    sigxHz(Nx-Np+1:Nx,i,:) = sighz;
    sigxEx(Nx-Np+1:Nx,i,:) = sigex;
    sigxEy(Nx-Np+1:Nx+1,i,:) = sigey;
    sigx(Nx-Np+1:Nx+1,i,:) = sigey;
    sigxEz(Nx-Np+1:Nx+1,i,:) = sigez;
end
sigxHy(Nx-Np+1:Nx,Ny+1,:) = sighy;
sigxEx(Nx-Np+1:Nx,Ny+1,:) = sigex;
sigx(Nx-Np+1:Nx+1,Ny+1,:) = sigey;
sigxEz(Nx-Np+1:Nx+1,Ny+1,:) = sigez;

clear sighx sighy sighz sigex sigey sigez