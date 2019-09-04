% intrinsic viscosity computing script
wd = pwd;
libloc = [wd '/../commoncode']; % put the folder "commoncode" on the path
addpath(libloc)

rerunsimulation = false

lsds = [0.01 1 100];            % Range of Saffman-Delbruck lengths to be simulated (a = 1 is assumed)
epsilonfactor = 0.5;            % This keeps epsilon to be half the spacing
spacings = 0.05:0.03:0.15;      % spacings 
angles = linspace(0,2*pi,11);   % angles to integrate over
angles =angles(1:end-1);        % so that we don't double-count 0 here
vext = @(x,y) [x,-y];

%henle_levine_interp = @(e)(12/pi^2)*log(1+e)+(3*pi*pi+(3*pi*pi+8*pi-12)*e+4*pi*e.^2)./(pi*pi*(1+e)); % function of epsilon = R/Lsd

nmonomers = [1 2 3 4 5 6];

alphas = NaN*ones(length(nmonomers),length(lsds));

if(rerunsimulation)
clf
for cc = 1:length(nmonomers)
    
    Sij_ext_all = cell(1,length(angles));
    Sijs_all = cell(1,length(angles));
    
    for iii = 1:length(angles)
        cmloc = [0 0];
        shape_function = @(s,R) smooth_pair(s,R,0,nmonomers(cc));
        [Fx_ext,Fy_ext,Lz_ext,Sij_ext,Fxs,Fys,Lzs,Sijs,xs,ys,ux,uy,fx,fy,Mij,Vxs,Vys,Omegas,Vx_ext,Vy_ext,Omega_ext] = regstokes_shape_function_compute_V_Omega(lsds,epsilonfactor,spacings,vext,cmloc,angles(iii),shape_function);
        plot(xs,ys,'ko'); axis equal;
        title(sprintf('theta = %3.3g, Sij=%3.3g %3.3g %3.3g %3.3g',angles(iii),Sij_ext(1),Sij_ext(2),Sij_ext(3),Sij_ext(4)));
        drawnow
        Sij_ext_all{iii} = Sij_ext;
        Sijs_all{iii} = Sijs;
    end
    
    %% Do average and plot
    Sij_ext_rotaverage = NaN*ones(size(Sij_ext));
    Sijs_rotaverage = zeros(size(Sijs));
    
    for iii = 1:length(angles)
        Sijs_rotaverage = Sijs_rotaverage+Sijs_all{iii};
    end
    
    Sijs_rotaverage = Sijs_rotaverage/length(angles);
    
    for i = 1:length(lsds)
        for jjj = 1:4
            p = polyfit(spacings,squeeze(Sijs_rotaverage(jjj,:,i)),1);
            Sij_ext_rotaverage(jjj,i) = p(2);
        end
    end
    
    Aprotein = nmonomers(cc)*pi; %R = 1
    alpha_est = Sij_ext_rotaverage(4,:)/(2*Aprotein);
    
    alphas(cc,:) = alpha_est
end

else % if you don't rerun the simulation, just load pre-existing data
   load('oligomer_data.mat'); 
   fprintf('Loading data from file rather than re-running whole simulation; change rerunsimulation to true to run code, will take long (many hours) time \n')
end

%%
clf
hold on
legendstring = {};
markerlist = {'o','x','d','s'};
if(nmonomers(1)~=1)
    error('Normalization of plot assumes nmonomers(1)=1, but this is not true')
end

for ll = 1:length(lsds)
    plot(nmonomers,alphas(:,ll)/alphas(1,ll),'-','LineWidth',4,'Marker',markerlist{ll},'MarkerSize',12);
    legendstring{ll} = sprintf('a/L_{sd} = %3.3g',1/lsds(ll));
end
set(gca,'FontSize',48)
xlabel('Number of oligomers')
ylabel('\alpha/\alpha_{monomer}')
legend(legendstring)

thf = linspace(0,2*pi,1e2);

%xlim([0 max(nms)])

%for j = 1:length(nms)
%    ax(j) = axes('Position',[0.235+0.1*(j-1) 0.3 0.1 0.5]);
ax = axes('Position',[0.15 0.7 0.5 0.2]);
hold on
for m = 1:nmonomers(end)
    pp = patch(cos(thf)+2*m,sin(thf),'k');
    set(pp,'LineWidth',3,'FaceColor',[0.8 0.8 0.8]);
end
axis equal
xlim([0 (max(nmonomers)+1)*2])
axis off