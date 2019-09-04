function [Ux,Uy,Omega,Sij,W,Fxcheck,Fycheck,Lzcheck] = find_rigids_from_GR(GR,Uxinf,Uyinf,Omegainf,Einf)
% [Ux,Uy,Omega,Sij,W,Fxcheck,Fycheck,Lzcheck] = find_rigids_from_GR(GR,Uxinf,Uyinf,Omegainf,Einf)
% find the rigid-body motions Ux,Uy,Omega from object with generalized
% resistance matrix GR, embedded in flow Uxinf, Uyinf, Omegainf,Einf
% also computes the stresslet induced

% create the 3x3 matrix relating velocities to forces, torques - this is
% the submatrix of the grand resistance matrix

W = [GR.A(1,1) GR.A(1,2) GR.Bt(1) ;
     GR.A(2,1) GR.A(2,2) GR.Bt(2) ;
     GR.B(1)   GR.B(2)   GR.C];

% now we need to compute the constant term... 
bx = sum(sum(squeeze(GR.Gt(1,:,:)).*Einf));
by = sum(sum(squeeze(GR.Gt(2,:,:)).*Einf));
bOmega = sum(sum(GR.Ht.*Einf));


% the rigid velocities come from
% W*(U - Uinf) = b, basically
% 
VVinf = [Uxinf ; Uyinf; Omegainf];
btot = [bx; by ; bOmega] + W*VVinf;
VVs = linsolve(W,btot);
Ux = VVs(1);
Uy = VVs(2);
Omega = VVs(3);

% now to find the value of SIJ, we need to solve
% S = G*(Uinf-U) + H*(Omegainf-Omega) + M*Einf
% so this just requires doing the contractions

Sij = NaN*ones(2);
for aa = 1:2
    for bb = 1:2
        Sij(aa,bb) = sum(squeeze(GR.G(aa,bb,:)).*(VVinf(1:2)-VVs(1:2))) + GR.H(aa,bb)*(Omegainf-Omega) + sum(sum(squeeze(GR.M(aa,bb,:,:)).*Einf));
    end
end

% now compute the forces and torques to make sure that they are zero!

Fxcheck = sum((GR.A(1,:).').*(VVinf(1:2)-VVs(1:2)))+GR.Bt(1)*(Omegainf-Omega)+sum(sum(squeeze(GR.Gt(1,:,:)).*Einf));
Fycheck = sum((GR.A(2,:).').*(VVinf(1:2)-VVs(1:2)))+GR.Bt(2)*(Omegainf-Omega)+sum(sum(squeeze(GR.Gt(2,:,:)).*Einf));
Lzcheck = sum(GR.B(:).*(VVinf(1:2)-VVs(1:2))) + GR.C*(Omegainf-Omega) + sum(sum(GR.Ht.*Einf));



end