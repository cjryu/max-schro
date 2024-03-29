ep = e0*ones(Nx+1,Ny+1,Nz+1);     % ep and mu defined on nodes.
mu = u0*ones(Nx+1,Ny+1,Nz+1);
epx = (ep(1:Nx,1:Ny+1,1:Nz+1)+ep(2:Nx+1,1:Ny+1,1:Nz+1))/2;      % ep on E grid.
epy = (ep(1:Nx+1,1:Ny,1:Nz+1)+ep(1:Nx+1,2:Ny+1,1:Nz+1))/2;
epz = (ep(1:Nx+1,1:Ny+1,1:Nz)+ep(1:Nx+1,1:Ny+1,2:Nz+1))/2;
mux = (mu(1:Nx+1,1:Ny,1:Nz)+mu(1:Nx+1,2:Ny+1,1:Nz)+mu(1:Nx+1,1:Ny,2:Nz+1)+mu(1:Nx+1,2:Ny+1,2:Nz+1))/4;      % mu on H grid.
muy = (mu(1:Nx,1:Ny+1,1:Nz)+mu(2:Nx+1,1:Ny+1,1:Nz)+mu(1:Nx,1:Ny+1,2:Nz+1)+mu(2:Nx+1,1:Ny+1,2:Nz+1))/4;
muz = (mu(1:Nx,1:Ny,1:Nz+1)+mu(2:Nx+1,1:Ny,1:Nz+1)+mu(1:Nx,2:Ny+1,1:Nz+1)+mu(2:Nx+1,2:Ny+1,1:Nz+1))/4;