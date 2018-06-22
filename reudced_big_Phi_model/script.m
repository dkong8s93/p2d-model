% Housekeeping
clear; clc;
 
% Assemble domain
h=1e-6; % x-spacing
ha=0.25e-6; % r-spacing
[Na,Np,Ns,Nn,Nz,Nr,r_half,r_whole]=assemble_domain(h,ha); % assemble domain 

% Indices for convenience
[id_cs,id_cn,id_ap,id_an,id_sp,id_sn,id_jp,id_jn,id_Phip,id_Phin,id_phip,id_phis,id_phin,id_Ta,id_Tp,id_Ts,id_Tn,id_Tz]=indices(Na,Np,Ns,Nn,Nr);

% Pre-determined parameters
tol=1e-4; % tolerance for Newton's method
k=2; % Timestepping
Tf=720; % End time
numsteps=floor(Tf/k);
Tref=298.15; % Environmental temperature
I=-150; % current density
g=1; % heat exchange coefficient

% Constants and constant matrix A
[d1,d2,d3,d4,d5,d6,d7,d13,d14,d15,d16,d17,d18,d19,d20,d21,d22,d23,d24,d25,d26,d27,d28,d29,d30,d31,d32,d33,d34,d35,d36,d37,d38,d39,d40,d41,d42,d43,d44,d45,d46]=assemble_constants(h,ha,k,I,g,Tref);
A=assemble_A(d2,d5,d20,d22,d23,d25,d30,d33,d38,d40,d41,d43,d45,Na,Np,Ns,Nn,Nz,Nr,r_half,r_whole);
    
% Initial condition
U0=[0;1000*ones(Np,1);0; 0;1000*ones(Ns,1);0; 0;1000*ones(Nn,1);0; 25751*ones(Np*Nr,1); 26128*ones(Nn*Nr,1); zeros(3*Np+Ns+3*Nn+2*3+1+1,1); 0;Tref*ones(Na,1);0; 0;Tref*ones(Np,1);0; 0;Tref*ones(Ns,1);0; 0;Tref*ones(Nn,1);0; 0;Tref*ones(Nz,1);0];

% First guess
[Phip,Phin]=findPhi(25751*ones(Np,1),26128*ones(Nn,1),zeros(Np,1),zeros(Nn,1),Tref*ones(Np,1),Tref*ones(Nn,1),Tref); 
U=[1000*ones(Np+Ns+Nn+2*3,1); 25751*ones(Np*Nr,1); 26128*ones(Nn*Nr,1); 25751*ones(Np,1); 26128*ones(Nn,1); zeros(Np+Nn,1); mean(Phip); mean(Phin); zeros(Np+Ns+Nn+2*3,1); Tref*ones(Na+Np+Ns+Nn+Nz+2*5,1)];

% Initialize plotting
[ap_plt,an_plt,cp_plt,cs_plt,cn_plt,jp_plt,jn_plt,Ta_plt,T_avg,Tz_plt,vol_plt]=initializePlotting(U,h,ha,numsteps);

% Time-stepping
counter=zeros(numsteps,1); % counter for iteration
U_storage=zeros(length(U),numsteps);

[v,Dv]=assemble_vDv(d1,d3,d4,d13,d14,d15,d16,d17,d18,d19,d21,d24,d27,d28,d29,d31,d32,d35,d36,d37,d39,d42,d43,d44,d45,d46,Tref,I,U,Na,Np,Ns,Nn,Nz,Nr,r_half,r_whole);
residual=A*U-U0+v;
J=A+Dv;
U=U-J\residual;
counter(1)=counter(1)+1;

for i=1:1:numsteps
    [v,Dv]=assemble_vDv(d1,d3,d4,d13,d14,d15,d16,d17,d18,d19,d21,d24,d27,d28,d29,d31,d32,d35,d36,d37,d39,d42,d43,d44,d45,d46,Tref,I,U,Na,Np,Ns,Nn,Nz,Nr,r_half,r_whole);
    residual=A*U-U0+v;
    while norm(residual,inf)>tol
        J=A+Dv;
        U=U-J\residual;
        [v,Dv]=assemble_vDv(d1,d3,d4,d13,d14,d15,d16,d17,d18,d19,d21,d24,d27,d28,d29,d31,d32,d35,d36,d37,d39,d42,d43,d44,d45,d46,Tref,I,U,Na,Np,Ns,Nn,Nz,Nr,r_half,r_whole);
        residual=A*U-U0+v;
        counter(i)=counter(i)+1;
    end
    U0=assemble_U0(U,Na,Np,Ns,Nn,Nz,Nr);
    U_storage(:,i)=U;
    i
    
    x=1:1:i;
    vol=U_storage(id_Phip+1,1:i)-U_storage(id_Phin+1,1:i);
    set(vol_plt,'Xdata',x);
    set(vol_plt,'Ydata',vol);
    set(Ta_plt,'Xdata',x);
    set(Ta_plt,'Ydata',U_storage(id_Ta+2,1:i));
    set(Tz_plt,'Xdata',x);
    set(Tz_plt,'Ydata',U_storage(id_Tz+2,1:i));
    set(T_avg,'Xdata',x);
    set(T_avg,'Ydata',mean(U_storage(id_Ta+1:id_Tz+Nz+1,1:i)));
   
%     set(ap_plt,'ZData',U(id_ap+1:id_ap+Np*Nr));
%     set(an_plt,'ZData',U(id_an+1:id_an+Nn*Nr));
    set(ap_plt,'ZData',flip(reshape(U(id_ap+1:id_ap+Np*Nr),Nr,Np)));
    set(an_plt,'ZData',flip(reshape(U(id_an+1:id_an+Nn*Nr),Nr,Nn)));
    set(cp_plt,'Ydata',U(2:Np+1));
    set(cs_plt,'Ydata',U(id_cs+2:id_cs+Ns+1));
    set(cn_plt,'Ydata',U(id_cn+2:id_cn+Nn+1));
    set(jp_plt,'Ydata',U(id_jp+1:id_jp+Np));
    set(jn_plt,'Ydata',U(id_jn+1:id_jn+Nn));
    drawnow
end