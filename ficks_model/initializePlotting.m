function [ap_plt,an_plt,cp_plt,cs_plt,cn_plt,jp_plt,jn_plt,Ta_plt,T_avg,Tz_plt,vol_plt]=initializePlotting(U,h,ha,numsteps)
    % domain
    x_0=1e-5; x_p=9e-5; x_s=11.5e-5; x_n=20.3e-5; R=2e-6;
    [Na,Np,Ns,Nn,Nz,Nr,~,~]=assemble_domain(h,ha); 
    [Xp,Yp]=meshgrid(x_0+h/2:h:x_p,ha/2:ha:R); [Xn,Yn]=meshgrid(x_s+h/2:h:x_n,ha/2:ha:R); 
    
    % indices for convenience
    [id_cs,id_cn,id_ap,id_an,id_sp,id_sn,id_jp,id_jn,id_Phip,id_Phin,id_phip,id_phis,id_phin,id_Ta,id_Tp,id_Ts,id_Tn,id_Tz]=indices(Na,Np,Ns,Nn,Nr);

    % variables
    cp=U(2:Np+1);
    cs=U(id_cs+2:id_cs+Ns+1); 
    cn=U(id_cn+2:id_cn+Nn+1); 
    ap_vec=U(id_ap+1:id_ap+Np*Nr); ap=flip(reshape(ap_vec,Nr,Np));
    an_vec=U(id_an+1:id_an+Nn*Nr); an=flip(reshape(an_vec,Nr,Nn));
    sp=U(id_sp+1:id_sp+Np);
    sn=U(id_sn+1:id_sn+Nn);
    
    Phip=U(id_Phip+1:id_Phip+Np); 
    Phin=U(id_Phin+1:id_Phin+Nn);
    phip=U(id_phip+2:id_phip+Np+1);
    phis=U(id_phis+2:id_phis+Ns+1);
    phin=U(id_phin+2:id_phin+Nn+1);
    
    jp=U(id_jp+1:id_jp+Np);
    jn=U(id_jn+1:id_jn+Nn);
    
    Ta=U(id_Ta+2:id_Ta+Na+1);
    Tp=U(id_Tp+2:id_Tp+Np+1);
    Ts=U(id_Ts+2:id_Ts+Ns+1);
    Tn=U(id_Tn+2:id_Tn+Nn+1);
    Tz=U(id_Tz+2:id_Tz+Nz+1);
    
    % spatial discretization
    xp=linspace(x_0,x_p,Np);
    xs=linspace(x_p,x_s,Ns);
    xn=linspace(x_s,x_n,Nn);

    figure('units','normalized','outerposition',[0 0 1 1])
    % Electrolyte concentration
    subplot(3,4,1);
    cp_plt=scatter(xp,cp,10,'*','y');
    hold on
    cs_plt=scatter(xs,cs,10,'*','m');
    cn_plt=scatter(xn,cn,10,'*','c');
    hold off
    ylabel('Concentration');
    xlim([x_0 x_n]);
    ylim([0 2000]);
    title('Electrolyte Concentration');
    % Ionic flux plot
    subplot(3,4,2);
    jp_plt=scatter(xp,jp,100,'.','y');
    hold on
    jn_plt=scatter(xn,jn,20,'+','c');
    hold off
    title('Ionic Flux')
    legend('cathode','anode','Location','SouthEast');
    ylabel('Ionic Flux');
    xlim([x_0 x_n]);
    % Temperature plot
    subplot(3,4,3);
    Ta_plt=scatter(1,Ta(1),100,'.','b');
    hold on
    Tz_plt=scatter(1,Tz(end),20,'+','g');
    T_avg=scatter(1,mean([Ta;Tp;Ts;Tn;Tz]),10,'*','y');
    hold off
    ylabel('Temperature');
    title('Temperature');   
    xlim([0 numsteps]);
    ylim([298 307]);
    legend('cathode','anode','average','Location','SouthEast');
     % Voltage plot
    subplot(3,4,4);
    vol_plt=scatter(1,Phip(1)-Phin(end),100,'.');
    title('Battery Voltage');
    xlabel('timestep');
    ylabel('Voltage (V)');
    ylim([2.5 4.4]);
    xlim([0 numsteps]); 
    % Solid concentration
    subplot(3,4,[5,6,7,8,9,10,11,12]);
%     an_plt=surf(Xn,Yn,an,'c','MarkerSize',3);
    an_plt=mesh(Xn,Yn,an);
    colorbar
    hold on
%     ap_plt=surf(Xp,Yp,ap,'y','MarkerSize',3);
    ap_plt=mesh(Xp,Yp,ap);
    hold off
    zlim([0 55000]);
    xlabel('x-position');
    ylabel('radius');
    zlabel('Concentration');
    title('Solid concentration');
end