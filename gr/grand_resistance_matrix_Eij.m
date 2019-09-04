function GRs = grand_resistance_matrix_Eij(lsds,epsilonfactor,spacings,angles,shapefunction)
% GRs = grand_resistance_matrix_Eij(lsds,epsilonfactor,spacings,angles,shapefunction)

GRs = cell(length(spacings),length(angles),length(lsds));

eta_m = 1;
R = 1;

for s = 1:length(spacings)
    spacing = spacings(s);
    epsilon = epsilonfactor*spacing;
    
    
    for t = 1:length(angles)
        
        %[xs,ys] = smoother_circle(spacing,R);
        [xs,ys] = shapefunction(spacing,R);
        xs = xs-mean(xs);  % watch out - this part is super sensitive if there is a zero rotation...
        ys = ys-mean(ys);
        [xs,ys] = rotate_pos(xs,ys,angles(t));
        clf
        plot(xs,ys,'ko');
        axis equal
        title(sprintf('Spacing = %3.3g, \\theta = %3.3g',spacing,angles(t)));
        drawnow
        
        for i = 1:length(lsds)
            GR = struct;
            GR.shapefunction = shapefunction;
            GR.A = NaN*ones(2,2);
            GR.B = NaN*ones(2,1);
            GR.Bt = NaN*ones(2,1);
            GR.C = NaN;
            GR.Ht = NaN*ones(2,2);
            GR.G = NaN*ones(2,2,2);
            GR.Gt = NaN*ones(2,2,2);
            GR.H = NaN*ones(2,2);
            GR.M = NaN*ones(2,2,2,2);
            
            Mij = reg_stokeslet_matrix(xs',ys',eta_m,lsds(i),epsilon);
            % First step: set Ux = 1, then read out all the elements
            % (note Ux here is the total velocity - careful with
            % conventions!)
            Ux = 1;
            Uy = 0;
            Omega = 0;
            Eij = [0 0 ; 0 0];
            [Fx,Fy,Lz,Sij] = get_forces(Mij,xs,ys,Ux,Uy,Omega,Eij);
            GR.A(1,1) = Fx;
            GR.A(2,1) = Fy;
            GR.B(1) = Lz;
            GR.G(:,:,1) = Sij;
            
            % Second step: set Uy = 1, read out the appropriate elements
            Ux = 0;
            Uy = 1;
            Omega = 0;
            Eij = [0 0 ; 0 0];
            [Fx,Fy,Lz,Sij] = get_forces(Mij,xs,ys,Ux,Uy,Omega,Eij);
            GR.A(1,2) = Fx;
            GR.A(2,2) = Fy;
            GR.B(2) = Lz;
            GR.G(:,:,2) = Sij;
            
            % Third step: rotation only to read out new elements
            Ux = 0;
            Uy = 0;
            Omega = 1;
            Eij = [0 0 ; 0 0];
            [Fx,Fy,Lz,Sij] = get_forces(Mij,xs,ys,Ux,Uy,Omega,Eij);
            GR.Bt(1) = Fx;
            GR.Bt(2) = Fy;
            GR.C = Lz;
            GR.H = Sij;
            
            % Fourth step: set E = [ [0 1 ; 1 0]] + read out Sij
            % We have a constraint here that, because E is symmetric,
            % Mij21+Mij12 = Sij but we can't identify
            % Mij21-Mij12 -- this is arbitrary; we can choose
            % Mij21 = Mij12
            % Again, be cautious about the sign of Eij here -- this is
            % not the background flow, but setting v_ext = - E dot r
            
            
            Ux = 0;
            Uy = 0;
            Omega = 0;
            Eij = [0 1 ; 1 0];
            %if(aa~=bb)
            [Fx,Fy,Lz,Sij] = get_forces(Mij,xs,ys,Ux,Uy,Omega,Eij);
            GR.M(:,:,2,1) = Sij/2;
            GR.M(:,:,1,2) = Sij/2;
            
            % now, it's a bit harder to figure out the other components of
            % Mijkl... we need to use that Mijkl = Mklij
            % we also constrain Mij11 + Mij22 = 0 (this is allowed
            % because of the arbitrariness in M, which, since it
            % multiplies E, which is trace-free, can have elements added to it...)
            % We could do the same thing with any arbitrary constant, and
            % it would give you the same result.
            %
            % Also, we know from the basic definition that Sij = Mijkl Ekl
            % that when Ekl = [1 0 ; 0 -1] that Sij = Mij11-Mij22
            % therefore adding/subtracting these two equations we get that
            % Mij11 = Sij/2
            % Mij22 = -Sij/2
            
            Ux = 0;
            Uy = 0;
            Omega = 0;
            Eij = [1 0 ; 0 -1];
            [Fx,Fy,Lz,Sij] = get_forces(Mij,xs,ys,Ux,Uy,Omega,Eij);
            arboffset = 0;
            GR.M(:,:,1,1) = Sij/2+arboffset;
            GR.M(:,:,2,2) = -Sij/2+arboffset;
            
            % need to find the others from symmetry relationships
            GR.Ht = GR.H;
            GR.Gt = NaN*ones(size(GR.G));
            for aa = 1:2
                for bb = 1:2
                    for cc = 1:2
                        GR.Gt(cc,aa,bb) = GR.G(aa,bb,cc);
                    end
                end
            end
            
            
            
            % Solve for the velocities and angular velocities of the particle
            % by ensuring that it is force- and torque-free -- this gives the
            % stresslet Sij imposed by
            Uxinf = 0;
            Uyinf = 0;
            Omegainf = 0;
            Einf = [1 0 ; 0 -1]; %
            
            [Ux,Uy,Omega,Sij] = find_rigids_from_GR(GR,Uxinf,Uyinf,Omegainf,Einf);
            GR.Uxinf = Uxinf;
            GR.Uyinf = Uyinf;
            GR.Einf = Einf;
            GR.Ux = Ux;
            GR.Uy = Uy;
            GR.Omega = Omega;
            GR.Sij = Sij;
            
            GRs{s,t,i} = GR;
            
            
        end
        
    end
    
end

end
