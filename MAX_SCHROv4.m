clc, clear, close all
disp(['Starting at: ',datestr(now)]);
% Optimized for speed from v2_3.

% Quantum Constants (Need to be in SI units!)
hbar = 1.054571726e-34; % Planck constant (J*s)
c0 = 2.99792458e8;      % Speed of light in vacuum (m/s)
m = 9.10938291e-31;     % Mass of an electron (kg)
me = 0.023*m/3;           % Effective mass (kg)
q = -1.60217657e-19;    % Charge of an electron (C)
% EM Constants
u0 = 12.566370614e-7;   % Permeability of vacuum (N/A^2)
e0 = 8.854187817e-12;   % Permittivity of vacuum (F/m)
n0 = sqrt(u0/e0);       % Impedance of free space (ohm)

% QD nodes (odd)
Y_qd = 6e-9;  Ny_qd = 28;  dy = Y_qd/(Ny_qd-1);
Nx_qd = 61;  Ny_qd = 61; Nz_qd = 95;  dx = dy;  dz = dy;
X_qd = dx*(Nx_qd-1);       Z_qd = dz*(Nz_qd-1);    Y_qd = dy*(Ny_qd-1);
% EM edges (even)
Nx = 66; Ny = 66; Nz = 100; Np = 10; Ns = 10; Np = Np+Ns;
Xo = Nx*dx; Yo = Ny*dy; Zo = Nz*dz;
X = Xo+Xo/Nx*2*Np; Y = Yo+Yo/Ny*2*Np; Z = Zo+Zo/Nz*2*Np;
Nx = Nx+2*Np; Ny = Ny+2*Np; Nz = Nz+2*Np;

x_qd = 0:dx:X_qd; y_qd = 0:dy:Y_qd; z_qd = 0:dz:Z_qd;
x = dx/2*(1-2*Np):dx:Xo-dx/2*(1-2*Np);
y = dy/2*(1-2*Np):dy:Yo-dy/2*(1-2*Np);
z = dz/2*(1-2*Np):dz:Zo-dz/2*(1-2*Np);

V = zeros(Nx_qd,Ny_qd,Nz_qd);
% 3D harmonic oscillator
lamb = 350e-9;  fo = c0/lamb;
parfor i = 1:Nx_qd
    for j = 1:Ny_qd
        for k = 1:Nz_qd
            % Use the effective mass for M-S.
            V(i,j,k) = .5*me*(2*pi*fo)^2*((x_qd(i)-dx*floor(Nx_qd/2))^2+(y_qd(j)-dy*floor(Ny_qd/2))^2+(z_qd(k)-dz*floor(Nz_qd/2))^2);
        end
    end
end
Vmax = max(max(max(V)));

Np = Np-Ns;     % Actual # of PMLs.
freespace       % Simulation setting applies to PML or non-PML simulations.
% dielectricblock
setupPML
averagednodesPML
Np = Np+Ns;     % Revert it back.

% The quantum dot is right in the middle of the simulation domain.
% Normally, Nx is even, and Nx_qd is odd. In the EM edge coordinates, the
% location of the QD corresponds to:
Nxstart = floor(Nx/2)-floor(Nx_qd/2)+1;  Nxend = floor(Nx/2)-floor(Nx_qd/2)+1+Nx_qd-1;
Nystart = floor(Ny/2)-floor(Ny_qd/2)+1;  Nyend = floor(Ny/2)-floor(Ny_qd/2)+1+Ny_qd-1;
Nzstart = floor(Nz/2)-floor(Nz_qd/2)+1;  Nzend = floor(Nz/2)-floor(Nz_qd/2)+1+Nz_qd-1;

Hx = zeros(Nx+1,Ny,Nz);     Hsyx = Hx;  sumHsyx = Hx;   Hszx = Hx;  sumHszx = Hx;
Hy = zeros(Nx,Ny+1,Nz);     Hsxy = Hy;  sumHsxy = Hy;   Hszy = Hy;  sumHszy = Hy;
Hz = zeros(Nx,Ny,Nz+1);     Hsxz = Hz;  sumHsxz = Hz;   Hsyz = Hz;  sumHsyz = Hz;
Ax = zeros(Nx,Ny+1,Nz+1);   Axm1 = Ax;  Axm3 = Ax;      Axm5 = Ax;  Axm7 = Ax;
Ay = zeros(Nx+1,Ny,Nz+1);   Aym1 = Ay;  Aym3 = Ay;      Aym5 = Ay;  Aym7 = Ay;
Az = zeros(Nx+1,Ny+1,Nz);   Azm1 = Az;  Azm3 = Az;      Azm5 = Az;  Azm7 = Az;
Axm0 = Ax; Aym0 = Ay; Azm0 = Az;
dHzdym1 = Ax;   dHzdym3 = Ax;   dHzdym5 = Ax;   dHydzm1 = Ax;   dHydzm3 = Ax;   dHydzm5 = Ax;
dHxdzm1 = Ay;   dHxdzm3 = Ay;   dHxdzm5 = Ay;   dHzdxm1 = Ay;   dHzdxm3 = Ay;   dHzdxm5 = Ay;
dHydxm1 = Az;   dHydxm3 = Az;   dHydxm5 = Az;   dHxdym1 = Az;   dHxdym3 = Az;   dHxdym5 = Az;
gradfxm1 = Ax;  gradfxm3 = Ax;  gradfxm5 = Ax;  sumgradfx = Ax; gradphix = Ax;  sumgradphix = Ax;
gradfym1 = Ay;  gradfym3 = Ay;  gradfym5 = Ay;  sumgradfy = Ay; gradphiy = Ay;  sumgradphiy = Ay;
gradfzm1 = Az;  gradfzm3 = Az;  gradfzm5 = Az;  sumgradfz = Az; gradphiz = Az;  sumgradphiz = Az;
f = zeros(Nx+1,Ny+1,Nz+1);  fsx = f;    fsy = f;    fsz = f;    sumfsx = f; sumfsy = f; sumfsz = f;
phisx = f;      phisxm0 = f;    phisxm1 = f;
phisy = f;      phisym0 = f;    phisym1 = f;
phisz = f;      phiszm0 = f;    phiszm1 = f;    phim0 = f;

% Parameters for the plane wave source.
% photons = 1;
% Ao = sqrt(photons*hbar*(2*pi*fo)/(X_qd*Y_qd*Z_qd)*c0*2*n0/(2*pi*fo)^2);
Ao = 1e-7;
Ninc = Ny-Np;
Ainc = zeros(1,Ninc+1); Aincm1 = Ainc; Aincm3 = Ainc;
Hincm1 = zeros(1,Ninc);   Hincm3 = Hincm1; dAzdym1 = Hincm1; dAzdym3 = Hincm1;
sigEinc = zeros(1,Ninc+1);  % HS starts at 2, ends at Ninc-Np. Corresponds to Np+2:Ny-Np in Ez. (-)
sigHinc = zeros(1,Ninc);    % HS starts at 1, ends at Ninc-Np. Corresponds to Np+1:Ny-Np in Hx. (+)
sigEinc(Ninc+1-Np:Ninc+1) = smaxy*(dy*(0:Np)/Ly).^m;
sigHinc(Ninc-Np+1:Ninc) = smaxy*((dy*(0:Np-1)+dy/2)/Ly).^m;

% Time steps
dt_qdnew = .9*hbar/(1/(2*me)*(4*hbar^2*(1/dx^2+1/dy^2+1/dz^2)-2*hbar*q*(Ao/dz)+q^2*Ao^2)+Vmax);
dt_qd = .9*hbar/(2*hbar^2/me*(1/dx^2+1/dy^2+1/dz^2)+Vmax);
dt = .9*1/(c0*sqrt(1/dx^2+1/dy^2+1/dz^2));
dt = min([dt_qd dt dt_qdnew]);
Nt = 9094; % 3 wavelengths
T = dt*(Nt-1);

insertHOGM_MS
pdfamp = .99*max(max(max(psi.*conj(psi))));
psihalf = (psim0+psim1)/2;
Jqxm1 = zeros(Nx_qd-1,Ny_qd,Nz_qd); Jqxm3 = Jqxm1; Jqxm5 = Jqxm3;
Jqym1 = zeros(Nx_qd,Ny_qd-1,Nz_qd); Jqym3 = Jqym1; Jqym5 = Jqym3;
Jqzm1 = zeros(Nx_qd,Ny_qd,Nz_qd-1); Jqzm3 = Jqzm1; Jqzm5 = Jqzm3;  rho = zeros(Nx_qd,Ny_qd,Nz_qd);

saveEz = zeros(2*Nt,4);     % There are two phases in this simulation.
saveAz = zeros(2*Nt,4);
saveJz = zeros(2*Nt,1);     % Save the electric current density.

for n = 0:Nt-1
    t = n*dt;
    
    updateMaxWindowedv2     % Update the vector potential up to n+1/2 and scalar potential up to n+1 (but only available up to nth step,
                            % n+1 can be made available).
    updateSchroAndJqv2      % Update the wavefunction up to n+1.
    % Transverse wave
    saveEz(n+1,1) = -(Az(floor(Nx/2),Np-floor(Ns/2),floor(Nz/2))-Azm1(floor(Nx/2),Np-floor(Ns/2),floor(Nz/2)))/dt ...
                      - (phim0(floor(Nx/2),Np-floor(Ns/2),floor(Nz/2)+1)-phim0(floor(Nx/2),Np-floor(Ns/2),floor(Nz/2)))/dz;
    % Longitudinal wave
    saveEz(n+1,2) = -(Az(floor(Nx/2),floor(Ny/2),Nz-Np+floor(Ns/2))-Azm1(floor(Nx/2),floor(Ny/2),Nz-Np+floor(Ns/2)))/dt ...
                      - (phim0(floor(Nx/2),floor(Ny/2),Nz-Np+floor(Ns/2)+1)-phim0(floor(Nx/2),floor(Ny/2),Nz-Np+floor(Ns/2)))/dz;
	% Excitation wave (transverse)
    saveEz(n+1,3) = -(Az(floor(Nx/2),Np+floor(Ns/2),floor(Nz/2))-Azm1(floor(Nx/2),Np+floor(Ns/2),floor(Nz/2)))/dt ...
                      - (phim0(floor(Nx/2),Np+floor(Ns/2),floor(Nz/2)+1)-phim0(floor(Nx/2),Np+floor(Ns/2),floor(Nz/2)))/dz;
	% Excitation wave (longitudinal)
    saveEz(n+1,4) = -(Az(floor(Nx/2),floor(Ny/2),Nz-Np-floor(Ns/2))-Azm1(floor(Nx/2),floor(Ny/2),Nz-Np-floor(Ns/2)))/dt ...
                      - (phim0(floor(Nx/2),floor(Ny/2),Nz-Np-floor(Ns/2)+1)-phim0(floor(Nx/2),floor(Ny/2),Nz-Np-floor(Ns/2)))/dz;
	saveAz(n+1,1) = Az(floor(Nx/2),Np-floor(Ns/2),floor(Nz/2));
    saveAz(n+1,2) = Az(floor(Nx/2),floor(Ny/2),Nz-Np+floor(Ns/2));
    saveAz(n+1,3) = Az(floor(Nx/2),Np+floor(Ns/2),floor(Nz/2));
    saveAz(n+1,4) = Az(floor(Nx/2),floor(Ny/2),Nz-Np-floor(Ns/2));
    saveJz(n+1) = sum(sum(sum(Jqzm1)));
    
    if(mod(n,10) == 0)
        pdf = psi.*conj(psi); % Generate the pdf for plotting.
        subplot(2,2,1);
        imagesc(x,y,Az(:,:,floor(Nz/2)+1)');
        axis xy;
        axis equal;
        colormap('jet');
        caxis([-Ao Ao]);
        subplot(2,2,2);
        pMax(:,:) = Az(floor(Nx/2)+1,:,:);
        imagesc(y,z,pMax');
        axis xy;
        axis equal;
        colormap('jet');
        caxis([-Ao Ao]);
        subplot(2,2,3);
        imagesc(x_qd,y_qd,pdf(:,:,floor(Nz_qd/2))');
        axis xy;
        axis equal;
        colormap('jet');
        caxis([0,pdfamp]);
        subplot(2,2,4);
        pSchro(:,:) = pdf(floor(Nx_qd/2),:,:);
        imagesc(y_qd,z_qd,pSchro');
        axis xy;
        axis equal;
        colormap('jet');
        caxis([0,pdfamp]);
        pause(.0001);
    end
end

Ao = 2e-13;   % Scale the wave down for plotting scattered wave.
for n = Nt:2*Nt-1       % Second phase of the simulation where the excited electron oscillates and emits wave with no further excitation.
    t = n*dt;
    
    updateMaxWindowedv2   % Update Maxwell part up to n+1/2.
    updateSchroAndJqv2    % Update the wavefunction up to n+1.
    % Transverse wave
    saveEz(n+1,1) = -(Az(floor(Nx/2),Np-floor(Ns/2),floor(Nz/2))-Azm1(floor(Nx/2),Np-floor(Ns/2),floor(Nz/2)))/dt ...
                      - (phim0(floor(Nx/2),Np-floor(Ns/2),floor(Nz/2)+1)-phim0(floor(Nx/2),Np-floor(Ns/2),floor(Nz/2)))/dz;
    % Longitudinal wave
    saveEz(n+1,2) = -(Az(floor(Nx/2),floor(Ny/2),Nz-Np+floor(Ns/2))-Azm1(floor(Nx/2),floor(Ny/2),Nz-Np+floor(Ns/2)))/dt ...
                      - (phim0(floor(Nx/2),floor(Ny/2),Nz-Np+floor(Ns/2)+1)-phim0(floor(Nx/2),floor(Ny/2),Nz-Np+floor(Ns/2)))/dz;
	% Excitation wave (transverse)
    saveEz(n+1,3) = -(Az(floor(Nx/2),Np+floor(Ns/2),floor(Nz/2))-Azm1(floor(Nx/2),Np+floor(Ns/2),floor(Nz/2)))/dt ...
                      - (phim0(floor(Nx/2),Np+floor(Ns/2),floor(Nz/2)+1)-phim0(floor(Nx/2),Np+floor(Ns/2),floor(Nz/2)))/dz;
	% Excitation wave (longitudinal)
    saveEz(n+1,4) = -(Az(floor(Nx/2),floor(Ny/2),Nz-Np-floor(Ns/2))-Azm1(floor(Nx/2),floor(Ny/2),Nz-Np-floor(Ns/2)))/dt ...
                      - (phim0(floor(Nx/2),floor(Ny/2),Nz-Np-floor(Ns/2)+1)-phim0(floor(Nx/2),floor(Ny/2),Nz-Np-floor(Ns/2)))/dz;
	saveAz(n+1,1) = Az(floor(Nx/2),Np-floor(Ns/2),floor(Nz/2));
    saveAz(n+1,2) = Az(floor(Nx/2),floor(Ny/2),Nz-Np+floor(Ns/2));
    saveAz(n+1,3) = Az(floor(Nx/2),Np+floor(Ns/2),floor(Nz/2));
    saveAz(n+1,4) = Az(floor(Nx/2),floor(Ny/2),Nz-Np-floor(Ns/2));
    saveJz(n+1) = sum(sum(sum(Jqzm1)));
   
    if(mod(n,10) == 0)
        pdf = psi.*conj(psi); % Generate the pdf for plotting.
        subplot(2,2,1);
        imagesc(x,y,Az(:,:,floor(Nz/2)+1)');
        axis xy;
        axis equal;
        colormap('jet');
        caxis([-Ao Ao]);
        subplot(2,2,2);
        pMax(:,:) = Az(floor(Nx/2)+1,:,:);
        imagesc(y,z,pMax');
        axis xy;
        axis equal;
        colormap('jet');
        caxis([-Ao Ao]);
        subplot(2,2,3);
        imagesc(x_qd,y_qd,pdf(:,:,floor(Nz_qd/2))');
        axis xy;
        axis equal;
        colormap('jet');
        caxis([0,pdfamp]);
        subplot(2,2,4);
        pSchro(:,:) = pdf(floor(Nx_qd/2),:,:);
        imagesc(y_qd,z_qd,pSchro');
        axis xy;
        axis equal;
        colormap('jet');
        caxis([0,pdfamp]);
        pause(.0001);
    end
end

disp(['Finishing at: ',datestr(now)]);