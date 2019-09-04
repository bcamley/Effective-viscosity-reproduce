function GRcorr = generate_corrected_M(GRs)
GRcorr = cell(size(GRs));

for s = 1:length(GRs(:))
    GR = GRs{s};
    ABCmatrix = [ [GR.A(1,1) GR.A(1,2) GR.Bt(1)] ; 
                  [GR.A(2,1) GR.A(2,2) GR.Bt(2)] ;
                  [GR.B(1)   GR.B(2)   GR.C    ] ];
    GHtransmatrix = [ [GR.Gt(1,1,1) GR.Gt(1,1,2) GR.Gt(1,2,1) GR.Gt(1,2,2)];
                      [GR.Gt(2,1,1) GR.Gt(2,1,2) GR.Gt(2,2,1) GR.Gt(2,2,2)];
                      [GR.Ht(1,1)   GR.Ht(1,2)   GR.Ht(2,1)   GR.Ht(2,2)  ]];
                        % This multiplies a vector [Exx ; Exy ; Eyx ; Eyy]
                        
    GHmatrix =      [[GR.G(1,1,1) GR.G(1,1,2) GR.H(1,1)] ;
                     [GR.G(1,2,1) GR.G(1,2,2) GR.H(1,2)] ;
                     [GR.G(2,1,1) GR.G(2,1,2) GR.H(2,1)] ;
                     [GR.G(2,2,1) GR.G(2,2,2) GR.H(2,2)] ];
                 
    MMatrix = -(GHmatrix)*((inv(ABCmatrix))*GHtransmatrix);
    
    GR.Mcorr = NaN*ones(2,2,2,2);
    
    Mi = [ [ {[1 1 1 1]}, {[1 1 1 2]}, {[1 1 2 1]}, {[1 1 2 2]}] ;
           [ {[1 2 1 1]}, {[1 2 1 2]}, {[1 2 2 1]}, {[1 2 2 2]}] ;
           [ {[2 1 1 1]}, {[2 1 1 2]}, {[2 1 2 1]}, {[2 1 2 2]}] ;
           [ {[2 2 1 1]}, {[2 2 1 2]}, {[2 2 2 1]}, {[2 2 2 2]}] ];
    
       
   for k = 1:length(Mi(:))
      ii = Mi{k};
      GR.Mcorr(ii(1),ii(2),ii(3),ii(4)) = MMatrix(k);
       
   end
   GRcorr{s} = GR;
end