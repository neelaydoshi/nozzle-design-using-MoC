%%%%%%%%%%%%%%
% Author : Neelay Doshi
% Project: Using Method of Characteristics for SERN Design
%%%%%%%%%%%%%%

function M = InversePrandtlMeyer(M1, M2, nu, gamma)

nu_1= PrandtlMeyer(M1, gamma);
nu_2= PrandtlMeyer(M2, gamma);

M_arr= [M1, M2];
nu_arr= [nu_1, nu_2];
y_arr= nu_arr - nu;

err= 1;
n= 3;

while err>1e-4
    
    M_n= M_arr(n-1) - y_arr(n-1)*(M_arr(n-1) - M_arr(n-2))/( y_arr(n-1) - y_arr(n-2) );
    nu_n= PrandtlMeyer(M_n, gamma);
    
    M_arr(n)= M_n;
    nu_arr(n)= nu_n;
    y_arr(n)= nu_n - nu;
    
    err= abs( ( nu_arr(n) - nu_arr(n-1) )/nu_arr(n-1) );
    n= n+1;
    
end

M= M_arr(n-1);



