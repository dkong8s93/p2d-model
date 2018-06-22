% function [ap_plt,an_plt,cp_plt,cs_plt,cn_plt,jp_plt,jn_plt,Phip_plt,Phin_plt,Ta_plt,Tp_plt,Ts_plt,Tn_plt,Tz_plt,vol_plt]=initializePlotting(U,Na,Np,Ns,Nn,Nz,x_le,x_0,x_p,x_s,x_n,x_re,numsteps)
function [ap_plt,an_plt,cp_plt,cs_plt,cn_plt,jp_plt,jn_plt,Phip_plt,Phin_plt,Ta_plt,T_avg,Tz_plt,vol_plt]=initializePlotting(U,Na,Np,Ns,Nn,Nz,x_le,x_0,x_p,x_s,x_n,x_re,numsteps)
    % indices for convenience
    [id_cs,id_cn,id_ap,id_an,~,~,id_jp,id_jn,id_Phip,id_Phin,~,~,~,id_Ta,id_Tp,id_Ts,id_Tn,id_Tz]=indices(Na,Np,Ns,Nn);

    % variables
    cp=U(2:Np+1);
    cs=U(id_cs+2:id_cs+Ns+1); 
    cn=U(id_cn+2:id_cn+Nn+1); 
    ap=U(id_ap+1:id_ap+Np); 
    an=U(id_an+1:id_an+Nn); 
    
    Phip=U(id_Phip+1:id_Phip+Np); 
    Phin=U(id_Phin+1:id_Phin+Nn);
    jp=U(id_jp+1:id_jp+Np);
    jn=U(id_jn+1:id_jn+Nn);
    
    Ta=U(id_Ta+2:id_Ta+Na+1);
    Tp=U(id_Tp+2:id_Tp+Np+1);
    Ts=U(id_Ts+2:id_Ts+Ns+1);
    Tn=U(id_Tn+2:id_Tn+Nn+1);
    Tz=U(id_Tz+2:id_Tz+Nz+1);
    
    % spatial discretization
    xa=linspace(x_le,x_0,Na);
    xp=linspace(x_0,x_p,Np);
    xs=linspace(x_p,x_s,Ns);
    xn=linspace(x_s,x_n,Nn);
    xz=linspace(x_n,x_re,Nz);

    figure('units','normalized','outerposition',[0 0 1 1])
    
    subplot(3,4,[3,4]);
    yyaxis left
    ap_plt=scatter(xp,ap,100,'.','b');
    hold on
    an_plt=scatter(xn,an,20,'+','g');
    hold off
    ylabel('Concentration');
    ylim([0 50000]);
    yyaxis right
    cp_plt=scatter(xp,cp,10,'*','y');
    hold on
    cs_plt=scatter(xs,cs,10,'*','m');
    cn_plt=scatter(xn,cn,10,'*','c');
    hold off
    ylabel('Concentration');
    xlim([x_0 x_n]);
    ylim([0 2000]);
    legend('a_p','a_n','c_p','c_s','c_n');
    title('Electrolyte and Solid Concentration');
    
    subplot(3,4,[7,8]);
    jp_plt=scatter(xp,jp,100,'.','b');
    hold on
    jn_plt=scatter(xn,jn,20,'+','g');
    hold off
    title('Ionic Flux')
    legend('cathode','anode','Location','SouthEast');
    ylabel('Ionic Flux');
    xlim([x_0 x_n]);
    
    subplot(3,4,[11,12]);
    Phip_plt=scatter(xp,Phip,100,'.','b');
    hold on
    Phin_plt=scatter(xn,Phin,20,'+','g');
    hold off
    title('Solid Potential');
    legend('cathode','anode');
    ylabel('Solid Potential');
    xlim([x_0 x_n]);
%     ylim([0 4.4]);
    
%     subplot(3,4,[1,2]);
%     Ta_plt=scatter(xa,Ta,100,'.','b');
%     hold on
%     Tp_plt=scatter(xp,Tp,20,'+','g');
%     Ts_plt=scatter(xs,Ts,10,'*','y');
%     Tn_plt=scatter(xn,Tn,10,'*','m');
%     Tz_plt=scatter(xz,Tz,10,'*','c');
%     hold off
%     title('Temperature');    
%     ylim([295 315]);

    subplot(3,4,[1,2]);
    Ta_plt=scatter(1,Ta(1),100,'.','b');
    hold on
    Tz_plt=scatter(1,Tz(end),20,'+','g');
    T_avg=scatter(1,mean([Ta;Tp;Ts;Tn;Tz]),10,'*','y');
    title('Temperature');   
    xlim([0 numsteps]);
    ylim([298 307]);
    legend('cathode','anode','average','Location','NorthWest');

    subplot(3,4,[5,6,9,10]);
    vol_plt=scatter(1,Phip(1)-Phin(end),100,'.');
    title('Battery Voltage');
    xlabel('timestep');
    ylabel('voltage (V)');
    ylim([2.5 4.4]);
    xlim([0 numsteps]);
end