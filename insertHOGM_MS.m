% Calculates and inserts the ground mode wavefunctions at different time
% steps for a harmonic oscillator. Also calculates the energy for the
% ground mode.

psi = zeros(Nx_qd,Ny_qd,Nz_qd);  psim0 = psi;  psim1 = psi;
omega = 2*pi*c0/lamb;   omega_g = 1.5*omega;   gndE = hbar*omega_g;
for i = 2:Nx_qd-1
    for j = 2:Ny_qd-1
        for k = 2:Nz_qd-1
            psim0(i,j,k) = (me*omega/(pi*hbar))^(1/4)*exp(-me*omega*(dx*(i-ceil(Nx_qd/2)))^2/(2*hbar))*hermite(0,sqrt(me*omega/hbar)*(dx*(i-ceil(Nx_qd/2))))...
                           *(me*omega/(pi*hbar))^(1/4)*exp(-me*omega*(dy*(j-ceil(Ny_qd/2)))^2/(2*hbar))*hermite(0,sqrt(me*omega/hbar)*(dy*(j-ceil(Ny_qd/2))))...
                           *(me*omega/(pi*hbar))^(1/4)*exp(-me*omega*(dz*(k-ceil(Nz_qd/2)))^2/(2*hbar))*hermite(0,sqrt(me*omega/hbar)*(dz*(k-ceil(Nz_qd/2))));
            psi(i,j,k) = psim0(i,j,k)*exp(-1i*omega_g*dt);
            psim1(i,j,k) = psim0(i,j,k)*exp(-1i*omega_g*(-dt));
        end
    end
end
psi = psi/sqrt(sum(sum(sum(psi.*conj(psi))))*dx*dy*dz);
psim0 = psim0/sqrt(sum(sum(sum(psim0.*conj(psim0))))*dx*dy*dz);
psim1 = psim1/sqrt(sum(sum(sum(psim1.*conj(psim1))))*dx*dy*dz);
