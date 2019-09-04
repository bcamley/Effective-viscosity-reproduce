function [GRa,GRext] = grand_average_and_extrapolate(GRs,spacings,angles,lsds)
% GRs should be cell w/variation along spacings, angles, lsds in that order

figmax = 8;
figcount = 1;

GRa = cell(length(spacings),length(lsds));
GRext = cell(length(lsds),1);

params = {'A','B','Bt','C','H','Ht','G','Gt','M','Sij','Mcorr'}; % the parameters to extrapolate to zero spacing

for l = 1:length(lsds)
    
    for s = 1:length(spacings)
        
        GR = struct;
        for p = 1:length(params)
            %params{p}
            GR.(params{p}) = zeros(size(GRs{1}.(params{p})));
        end
        
        for t = 1:length(angles)
            
            for p = 1:length(params)
                GR.(params{p}) = GR.(params{p}) + GRs{s,t,l}.(params{p})/length(angles);
            end
            
        end
        
        GRa{s,l} = GR;
        
    end
    
    % now we can extrapolate to zero spacing
    % subplot(ceil(sqrt(length(lsds))),ceil(sqrt(length(lsds))),l)
    
    
    GRe = struct;
    for p = 1:length(params)
         GRe.(params{p}) = zeros(size(GRs{1}.(params{p})));
         for j = 1:length(GRe.(params{p})(:))
            varying = zeros(size(spacings));
            for ss = 1:length(spacings)
                varying(ss) = GRa{ss,l}.(params{p})(j);
            end
            ppol = polyfit(spacings,varying,1);
            vfit = polyval(ppol,spacings);
            vresid = varying-vfit;
            SSresid = sum(vresid.^2);
            SStotal = (length(varying)-1)*var(varying);
            rsq = 1-SSresid/SStotal;
            if((rsq<0.4)&&(abs(ppol(2))>5e-2))
                warning('Extrapolation not working well for some variables... see figure')
                if(figcount<figmax)
                    figure
                    plot(spacings,varying,'ko');
                    hold on;
                    sf = linspace(0,max(spacings),1e3);
                    plot(sf,polyval(ppol,sf),'b--');
                    title(sprintf('rsq = %3.3g param = %s no %d lsd = %3.3g',rsq,params{p},j,lsds(l)));
                    figcount = figcount + 1;
                else
                    warning(sprintf('Too many figures being generated... but have error at rsq = %3.3g param = %s no %d lsd = %3.3g',rsq,params{p},j,lsds(l)));
                    
                end
            end
            GRe.(params{p})(j) = ppol(2);
         end
    end
    GRext{l} = GRe;
    
        
end


end