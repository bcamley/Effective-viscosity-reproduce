% intrinsic viscosity computing script
wd = pwd;
libloc = [wd '/../commoncode'];
addpath(libloc)
lsds = logspace(-2,2,25)
epsilonfactor = 0.5
spacings = 0.05:0.05:0.4
vext = @(x,y) [x,-y]

rerunsimulation = false


henle_levine_interp = @(e)(12/pi^2)*log(1+e)+(3*pi*pi+(3*pi*pi+8*pi-12)*e+4*pi*e.^2)./(pi*pi*(1+e)); % function of epsilon = R/Lsd

%vext = @(x,y) [

if(rerunsimulation)
    [Fx_ext,Fy_ext,Lz_ext,Sij_ext,Fxs,Fys,Lzs,Sijs,xs,ys,ux,uy,fx,fy,Mij,Vxs,Vys,Omegas] = regstokes_circle_refiner_compute_V_Omega(lsds,epsilonfactor,spacings,vext);
else
   load('intrinsic_viscosity_circle_data.mat')
   fprintf('Loading data from file rather than re-running whole simulation; change rerunsimulation to true to run code, will take 10-15 minutes \n')
end

alpha_est = Sij_ext(4,:)/(2*pi);

%% Plotting section

clf
h_do=plot(1./lsds(lsds>0.1),2*ones(size(lsds(lsds>0.1))),'k--','LineWidth',2);
hold on
h_hl = plot(1./lsds,henle_levine_interp(1./lsds),'r','LineWidth',4);
h_hlscale = plot(1./lsds,(2/3)*henle_levine_interp(1./lsds),'b','LineWidth',4);

h_rs = plot(1./lsds,alpha_est,'ko','MarkerSize',14,'LineWidth',3);

set(gca,'FontSize',40);
xlabel('a/L_{sd}'); ylabel('Intrinsic viscosity \alpha')

legend([h_rs h_do h_hl h_hlscale],{'Regularized Stokeslet calculation','Diamant-Oppenheimer (a << L_{sd})','Henle-Levine interpolant','(2/3)xHenle-Levine'})

set(gca,'xscale','log','yscale','log')