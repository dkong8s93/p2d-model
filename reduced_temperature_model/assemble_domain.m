function [Np,Ns,Nn,Nr,r_half,r_whole]=assemble_domain(h,ha)
    x_0=1e-5; x_p=9e-5; x_s=11.5e-5; x_n=20.3e-5; R=2e-6;
    Np=floor((x_p-x_0)/h); % number of points for section p
    Ns=floor((x_s-x_p)/h); % number of points for section s
    Nn=floor((x_n-x_s)/h); % number of points for section n
    Nr=R/ha; % number of points inside a solid particle
    
    r=0:ha:2e-6;
    r_half=r+ha; r_half=(r_half(1:end-1))';
    r_whole=r+ha/2; r_whole=(r_whole(1:end-1))';
end