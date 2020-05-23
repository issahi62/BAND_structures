function C = convmat(A, P, Q, R)


%% Gridsize of the input Array
[Nx, Ny, Nz] = size(A); 
 

%% Handle the number spation dispersion (1D, 2D, 3D); 
if nargin  == 2 
    Q =1; 
    R =1; 
elseif nargin == 3
    R = 1;
end 

NH = P*Q*R;

p= [-floor(P/2): floor(P/2)]; %indices along x
q= [-floor(Q/2): floor(Q/2)]; 
r= [-floor(R/2): floor(R/2)]; 


A  = fftshift(fftn(A))/(Nx*Ny*Nz); 


pO = 1+floor(Nx/2); 
qO = 1+ floor(Ny/2);
rO = 1+floor(Nz/2);


for rrow = 1:R
for qrow = 1:Q
for prow = 1:P
    row = (rrow-1)*Q*P + (qrow-1)*P +prow;
    for rcol = 1:R
    for qcol = 1:Q
    for pcol = 1:P
    col = (rcol-1)*Q*P + (qcol-1)*P +pcol;
    
    pfft = p(prow)-p(pcol);
    qfft = q(qrow)-q(qcol);
    rfft = r(rrow)-r(rcol); 

    C(row,col) = A(pO+pfft, qO+qfft, rO+rfft);
    end
    end
    end
end
end
end 
end 



    
  