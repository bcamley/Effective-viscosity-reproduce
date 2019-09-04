%% generate the grand resistance matrix, then average and extrapolate it
wd = pwd;
libloc = [wd '/../commoncode']; % put the folder "commoncode" on the path
addpath(libloc)

lsds = [0.01 1 100];
%lsds = 1
epsilonfactor = 0.5

%spacings = 0.08
spacings = 0.05:0.03:0.15;

angles = 0;  % for the GR method, we don't have to average rotationally
fprintf('mean of cosine of angles = %3.3g, should be about 0 \n',mean(cos(angles)))

nms = 1:6
GRnms = cell(size(nms));
alphas = NaN*ones(length(nms),length(lsds));

for i = 1:length(nms)
    
    nm = nms(i);
    shapefunction = @(s,R) smooth_pair(s,R,0,nm)
    
    GRs = grand_resistance_matrix_Eij(lsds,epsilonfactor,spacings,angles,shapefunction);
    GRs = generate_corrected_M(GRs);  % this gives the "Mcorr" term
    [GRa,GRext] = grand_average_and_extrapolate(GRs,spacings,angles,lsds);
    
    GRnms{i} = GRs;
    Einf = [ [1 0] ; [0 -1]];
    for ll = 1:length(lsds)
        M = GRext{ll}.M + GRext{ll}.Mcorr;
        S11 = (1/4)*(M(1,1,1,1)-M(1,1,2,2)+M(1,2,1,2)+M(1,2,2,1)+M(2,1,1,2)+M(2,1,2,1)-M(2,2,1,1)+M(2,2,2,2));
        
        % This code would allow for computing with angular averaging to make sure
        % this route is consistent
        %         Sij = NaN*ones(2);
        %            for aa = 1:2
        %                for bb = 1:2
        %                    Sij(aa,bb) = sum(sum(squeeze(GRext{ll}.M(aa,bb,:,:)).*Einf))+ sum(sum(squeeze(GRext{ll}.Mcorr(aa,bb,:,:)).*Einf));
        %                end
        %            end
        %            alphas(i,ll) = Sij(1)/(2*pi*nm);
        alphas(i,ll) = S11/(2*pi*nm);
    end
    %     for ll = 1:length(lsds)
    %         alphas(i,ll) = GRext{ll}.Sij(1)/(2*pi*nm);
    %     end
    
end
%% Plot intrinsic viscosity

clf
hold on
legendstring = {};
if(nms(1)~=1)
    warning('Normalization assumes nms(1)=1, but this is not true')
end

for ll = 1:length(lsds)
    plot(nms,alphas(:,ll)/alphas(1,ll),'-o','LineWidth',4);
    legendstring{ll} = sprintf('a/L_{sd} = %3.3g',1/lsds(ll));
end
set(gca,'FontSize',48)
xlabel('Number of oligomers')
ylabel('\alpha/\alpha_{monomer}')
legend(legendstring)

thf = linspace(0,2*pi,1e2);

ax = axes('Position',[0.15 0.7 0.5 0.2]);
hold on
for m = 1:nms(end)
    pp = patch(cos(thf)+2*m,sin(thf),'k');
    set(pp,'LineWidth',3,'FaceColor',[0.8 0.8 0.8]);
end
axis equal
xlim([0 (max(nms)+1)*2])
axis off
%end
