function [ap_plt,an_plt,cp_plt,cs_plt,cn_plt,jp_plt,jn_plt,T_plt,vol_plt]=initializePlotting(U,h,ha,numsteps)
    % domain
    x_0=1e-5; x_p=9e-5; x_s=11.5e-5; x_n=20.3e-5; R=2e-6;
    [Np,Ns,Nn,Nr,~,~]=assemble_domain(h,ha); 
    [Xp,Yp]=meshgrid(x_0+h/2:h:x_p,ha/2:ha:R); [Xn,Yn]=meshgrid(x_s+h/2:h:x_n,ha/2:ha:R); 
    
    % indices for convenience
    [id_cs,id_cn,id_ap,id_an,~,~,id_jp,id_jn,id_Phip,id_Phin,~,~,~,id_T]=indices(Np,Ns,Nn,Nr);

    % variables
    cp=U(2:Np+1);
    cs=U(id_cs+2:id_cs+Ns+1); 
    cn=U(id_cn+2:id_cn+Nn+1); 
    ap_vec=U(id_ap+1:id_ap+Np*Nr); ap=flip(reshape(ap_vec,Nr,Np));
    an_vec=U(id_an+1:id_an+Nn*Nr); an=flip(reshape(an_vec,Nr,Nn));
    
    Phip=U(id_Phip+1:id_Phip+Np); 
    Phin=U(id_Phin+1:id_Phin+Nn);
    
    jp=U(id_jp+1:id_jp+Np);
    jn=U(id_jn+1:id_jn+Nn);
    
    T=U(id_T+1);
    
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
    T_plt=scatter(1,T,100,'.','b');
    ylabel('Temperature');
    title('Temperature');   
    xlim([0 numsteps]);
    ylim([290 330]);
    legend('temperature','Location','SouthEast');
    
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