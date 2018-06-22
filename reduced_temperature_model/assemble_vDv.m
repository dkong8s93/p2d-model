function [v,Dv]=assemble_vDv(d1,d3,d4,d15,d16,d17,d18,d19,d21,d26,d27,d28,d29,d31,d32,d34,d35,d36,d37,d39,d40,d41,d43,d44,d45,d46,h,Tref,I,U,Np,Ns,Nn,Nr,r_half,r_whole)
    dim=Np+2+Ns+2+Nn+2+Np*Nr+Nn*Nr+Np+Nn+Np+Nn+Np+Nn+Np+2+Ns+2+Nn+2+1; % dimension 
    
    % indices for convenience
    [id_cs,id_cn,id_ap,id_an,id_sp,id_sn,id_jp,id_jn,id_Phip,id_Phin,id_phip,id_phis,id_phin,id_T]=indices(Np,Ns,Nn,Nr);
    
    % variables
    cp_g=U(1:Np+2); Cp_g=log(cp_g);
    cs_g=U(id_cs+1:id_cs+Ns+2); Cs_g=log(cs_g);
    cn_g=U(id_cn+1:id_cn+Nn+2); Cn_g=log(cn_g);
    sp=U(id_sp+1:id_sp+Np); 
    sn=U(id_sn+1:id_sn+Nn); 
    jp=U(id_jp+1:id_jp+Np);
    jn=U(id_jn+1:id_jn+Nn);
    Phip=U(id_Phip+1:id_Phip+Np); 
    Phin=U(id_Phin+1:id_Phin+Nn);
    phip_g=U(id_phip+1:id_phip+Np+2);
    phis_g=U(id_phis+1:id_phis+Ns+2);
    phin_g=U(id_phin+1:id_phin+Nn+2);
    T=U(id_T+1);
    
    % Auxiliary functions
    Dp_g=Deffp(cp_g,T); % Diffusion coefficient for the cathode check
    Ds_g=Deffs(cs_g,T); % Diffusion coefficient for the separator check
    Dn_g=Deffn(cn_g,T); % Diffusion coefficient for the anode check
    kp_g=kappap(cp_g,T); % Conductivity for the cathode check
    ks_g=kappas(cs_g,T); % Conductivity for the separator check
    kn_g=kappan(cn_g,T); % Conductivity for the anode check
    DDp_cp=DDeffp_cp(cp_g,T); 
    DDs_cs=DDeffs_cs(cs_g,T);
    DDn_cn=DDeffn_cn(cn_g,T);
    DDp_DTp=DDeffp_Tp(cp_g,T);
    DDs_DTs=DDeffs_Ts(cs_g,T); 
    DDn_DTn=DDeffn_Tn(cn_g,T);
    Dkp_cp=Dkappap_cp(cp_g,T);
    Dks_cs=Dkappas_cs(cs_g,T);
    Dkn_cn=Dkappan_cn(cn_g,T);
    Dkp_Tp=Dkappap_Tp(cp_g,T);
    Dks_Ts=Dkappas_Ts(cs_g,T);
    Dkn_Tn=Dkappan_Tn(cn_g,T);

    % vector v
    v=zeros(dim,1); % check
    % cp_g
    v(2:Np+1)=d1*((Dp_g(2:Np+1)+Dp_g(3:Np+2)).*(cp_g(3:Np+2)-cp_g(2:Np+1))-(Dp_g(1:Np)+Dp_g(2:Np+1)).*(cp_g(2:Np+1)-cp_g(1:Np))); % check
    v(Np+2)=(Dp_g(Np+1)+Dp_g(Np+2))*(cp_g(Np+2)-cp_g(Np+1))-(Ds_g(1)+Ds_g(2))*(cs_g(2)-cs_g(1)); % check
    % cs_g
    v(id_cs+2:id_cs+Ns+1)=d3*((Ds_g(2:Ns+1)+Ds_g(3:Ns+2)).*(cs_g(3:Ns+2)-cs_g(2:Ns+1))-(Ds_g(1:Ns)+Ds_g(2:Ns+1)).*(cs_g(2:Ns+1)-cs_g(1:Ns))); % check
    % cn_g
    v(id_cn+1)=(Ds_g(Ns+1)+Ds_g(Ns+2))*(cs_g(Ns+2)-cs_g(Ns+1))-(Dn_g(1)+Dn_g(2))*(cn_g(2)-cn_g(1)); % check
    v(id_cn+2:id_cn+Nn+1)=d4*((Dn_g(2:Nn+1)+Dn_g(3:Nn+2)).*(cn_g(3:Nn+2)-cn_g(2:Nn+1))-(Dn_g(1:Nn)+Dn_g(2:Nn+1)).*(cn_g(2:Nn+1)-cn_g(1:Nn))); % check
    % ap
    v_ap_A=-d43*((r_half(end)/r_whole(end))^2)*apNL(jp,T,Tref,d44);
    v_ap_B=zeros(Nr,1); v_ap_B(end)=1;
    v(id_ap+1:id_ap+Np*Nr)=kron(v_ap_A,v_ap_B);
    % an
    v_an_A=-d45*((r_half(end)/r_whole(end))^2)*anNL(jn,T,Tref,d46);
    v_an_B=zeros(Nr,1); v_an_B(end)=1;
    v(id_an+1:id_an+Nn*Nr)=kron(v_an_A,v_an_B);
    % sp
    v(id_sp+1:id_sp+Np)=0.5*apNL(jp,T,Tref,d44);
    % sn
    v(id_sn+1:id_sn+Nn)=0.5*anNL(jn,T,Tref,d46);
    % jp
    v(id_jp+1:id_jp+Np)=jpNL(cp_g(2:Np+1),sp,Phip,phip_g(2:Np+1),T,Tref);
    % jn
    v(id_jn+1:id_jn+Nn)=jnNL(cn_g(2:Nn+1),sn,Phin,phin_g(2:Nn+1),T,Tref);
    % Phip
    v(id_Phip+1)=h*I/59;  % check
    % Phin
    v(id_Phin+Nn)=-h*I/48.24; % check
    % phip
    v(id_phip+2:id_phip+Np+1)=d15*((kp_g(2:Np+1)+kp_g(3:Np+2)).*(phip_g(3:Np+2)-phip_g(2:Np+1))-(kp_g(1:Np)+kp_g(2:Np+1)).*(phip_g(2:Np+1)-phip_g(1:Np)))+...
        d16*((kp_g(2:Np+1)+kp_g(3:Np+2)).*(2*T).*(Cp_g(3:Np+2)-Cp_g(2:Np+1))-(kp_g(1:Np)+kp_g(2:Np+1)).*(2*T).*(Cp_g(2:Np+1)-Cp_g(1:Np)));
    v(id_phip+Np+2)=(kp_g(Np+1)+kp_g(Np+2))*(phip_g(Np+2)-phip_g(Np+1))-(ks_g(1)+ks_g(2))*(phis_g(2)-phis_g(1)); % check
    % phis
    v(id_phis+2:id_phis+Ns+1)=(ks_g(2:Ns+1)+ks_g(3:Ns+2)).*(phis_g(3:Ns+2)-phis_g(2:Ns+1))-(ks_g(1:Ns)+ks_g(2:Ns+1)).*(phis_g(2:Ns+1)-phis_g(1:Ns))+...
        d17*((ks_g(2:Ns+1)+ks_g(3:Ns+2)).*(2*T).*(Cs_g(3:Ns+2)-Cs_g(2:Ns+1))-(ks_g(1:Ns)+ks_g(2:Ns+1)).*(2*T).*(Cs_g(2:Ns+1)-Cs_g(1:Ns)));
    % phin
    v(id_phin+1)=(ks_g(Ns+1)+ks_g(Ns+2))*(phis_g(Ns+2)-phis_g(Ns+1))-(kn_g(1)+kn_g(2))*(phin_g(2)-phin_g(1)); % check
    v(id_phin+2:id_phin+Nn+1)=d18*((kn_g(2:Nn+1)+kn_g(3:Nn+2)).*(phin_g(3:Nn+2)-phin_g(2:Nn+1))-(kn_g(1:Nn)+kn_g(2:Nn+1)).*(phin_g(2:Nn+1)-phin_g(1:Nn)))+...
        d19*((kn_g(2:Nn+1)+kn_g(3:Nn+2)).*(2*T).*(Cn_g(3:Nn+2)-Cn_g(2:Nn+1))-(kn_g(1:Nn)+kn_g(2:Nn+1)).*(2*T).*(Cn_g(2:Nn+1)-Cn_g(1:Nn)));
    % T
        v(id_T+1)=d21+...
            d26*((Phip(2)-Phip(1)-h*I/59)^2)+d27*kp_g(2)*((phip_g(3)-phip_g(1))^2)+d28*kp_g(2)*T*(Cp_g(3)-Cp_g(1))*(phip_g(3)-phip_g(1))+Qp(sp(1),jp(1),Phip(1),phip_g(2),T,Tref,d29)+...
            sum(d26*((Phip(3:Np)-Phip(1:Np-2)).^2)+d27*kp_g(3:Np).*((phip_g(4:Np+1)-phip_g(2:Np-1)).^2)+d28*kp_g(3:Np).*T.*(Cp_g(4:Np+1)-Cp_g(2:Np-1)).*(phip_g(4:Np+1)-phip_g(2:Np-1))+Qp(sp(2:Np-1),jp(2:Np-1),Phip(2:Np-1),phip_g(3:Np),T,Tref,d29))+...
            d26*((Phip(Np)-Phip(Np-1))^2)+d27*kp_g(Np+1)*((phip_g(Np+2)-phip_g(Np))^2)+d28*kp_g(Np+1)*T*(Cp_g(Np+2)-Cp_g(Np))*(phip_g(Np+2)-phip_g(Np))+Qp(sp(Np),jp(Np),Phip(Np),phip_g(Np+1),T,Tref,d29)+...
            sum(d31*ks_g(2:Ns+1).*((phis_g(3:Ns+2)-phis_g(1:Ns)).^2)+d32*ks_g(2:Ns+1).*T.*(Cs_g(3:Ns+2)-Cs_g(1:Ns)).*(phis_g(3:Ns+2)-phis_g(1:Ns)))+...
            d34*((Phin(2)-Phin(1))^2)+d35*kn_g(2)*((phin_g(3)-phin_g(1))^2)+d36*kn_g(2)*T*(Cn_g(3)-Cn_g(1))*(phin_g(3)-phin_g(1))+Qn(sn(1),jn(1),Phin(1),phin_g(2),T,Tref,d37)+...
            sum(d34*((Phin(3:Nn)-Phin(1:Nn-2)).^2)+d35*kn_g(3:Nn).*((phin_g(4:Nn+1)-phin_g(2:Nn-1)).^2)+d36*kn_g(3:Nn).*T.*(Cn_g(4:Nn+1)-Cn_g(2:Nn-1)).*(phin_g(4:Nn+1)-phin_g(2:Nn-1))+Qn(sn(2:Nn-1),jn(2:Nn-1),Phin(2:Nn-1),phin_g(3:Nn),T,Tref,d37))+...
            d34*(Phin(Nn)-h*I/48.24-Phin(Nn-1))^2+d35*kn_g(Nn+1)*((phin_g(Nn+2)-phin_g(Nn))^2)+d36*kn_g(Nn+1)*T*(Cn_g(Nn+2)-Cn_g(Nn))*(phin_g(Nn+2)-phin_g(Nn))+Qn(sn(Nn),jn(Nn),Phin(Nn),phin_g(Nn+1),T,Tref,d37)+...
            d39+d40*Tref+d41*Tref;
        
    
    % matrix Dv
    Dv=zeros(dim);
    
    % Dcp_Dcp
    Dv(2:Np+1,1:Np+2)=spdiags([d1*(-DDp_cp(1:Np).*(cp_g(2:Np+1)-cp_g(1:Np))+Dp_g(1:Np)+Dp_g(2:Np+1)),...
        d1*(DDp_cp(2:Np+1).*(cp_g(3:Np+2)-cp_g(2:Np+1))-(Dp_g(2:Np+1)+Dp_g(3:Np+2))-DDp_cp(2:Np+1).*(cp_g(2:Np+1)-cp_g(1:Np))-(Dp_g(1:Np)+Dp_g(2:Np+1))),...
        d1*(DDp_cp(3:Np+2).*(cp_g(3:Np+2)-cp_g(2:Np+1))+(Dp_g(2:Np+1)+Dp_g(3:Np+2)))],0:2,Np,Np+2); % check
    Dv(Np+2,Np+1)=DDp_cp(Np+1)*(cp_g(Np+2)-cp_g(Np+1))-(Dp_g(Np+1)+Dp_g(Np+2));
    Dv(Np+2,Np+2)=DDp_cp(Np+2)*(cp_g(Np+2)-cp_g(Np+1))+Dp_g(Np+1)+Dp_g(Np+2);
    Dv(Np+2,id_cs+1)=-DDs_cs(1)*(cs_g(2)-cs_g(1))+Ds_g(1)+Ds_g(2);
    Dv(Np+2,id_cs+2)=-DDs_cs(2)*(cs_g(2)-cs_g(1))-(Ds_g(1)+Ds_g(2));
    
    % Dcp_DTp
    Dv(2:Np+1,id_T+1)=d1*((DDp_DTp(2:Np+1)+DDp_DTp(3:Np+2)).*(cp_g(3:Np+2)-cp_g(2:Np+1))-(DDp_DTp(1:Np)+DDp_DTp(2:Np+1)).*(cp_g(2:Np+1)-cp_g(1:Np)));
    Dv(Np+2,id_T+1)=(DDp_DTp(Np+1)+DDp_DTp(Np+2))*(cp_g(Np+2)-cp_g(Np+1))-(DDs_DTs(1)+DDs_DTs(2))*(cs_g(2)-cs_g(1));
    
    % Dcs_Dcs
    Dv(id_cs+2:id_cs+Ns+1,id_cs+1:id_cs+Ns+2)=spdiags([d3*(-DDs_cs(1:Ns).*(cs_g(2:Ns+1)-cs_g(1:Ns))+Ds_g(1:Ns)+Ds_g(2:Ns+1)),...
        d3*(DDs_cs(2:Ns+1).*(cs_g(3:Ns+2)-cs_g(2:Ns+1))-(Ds_g(2:Ns+1)+Ds_g(3:Ns+2))-DDs_cs(2:Ns+1).*(cs_g(2:Ns+1)-cs_g(1:Ns))-(Ds_g(1:Ns)+Ds_g(2:Ns+1)))...
        d3*(DDs_cs(3:Ns+2).*(cs_g(3:Ns+2)-cs_g(2:Ns+1))+Ds_g(2:Ns+1)+Ds_g(3:Ns+2))],0:2,Ns,Ns+2);
    
    % Dcs_DTs
    Dv(id_cs+2:id_cs+Ns+1,id_T+1)=d3*((DDs_DTs(2:Ns+1)+DDs_DTs(3:Ns+2)).*(cs_g(3:Ns+2)-cs_g(2:Ns+1))-(DDs_DTs(1:Ns)+DDs_DTs(2:Ns+1)).*(cs_g(2:Ns+1)-cs_g(1:Ns)));
    
    % Dcn_Dcn
    Dv(id_cn+1,id_cs+Ns+1)=DDs_cs(Ns+1)*(cs_g(Ns+2)-cs_g(Ns+1))-(Ds_g(Ns+1)+Ds_g(Ns+2));
    Dv(id_cn+1,id_cs+Ns+2)=DDs_cs(Ns+2)*(cs_g(Ns+2)-cs_g(Ns+1))+Ds_g(Ns+1)+Ds_g(Ns+2);
    Dv(id_cn+1,id_cn+1)=-DDn_cn(1)*(cn_g(2)-cn_g(1))+Dn_g(1)+Dn_g(2);
    Dv(id_cn+1,id_cn+2)=-DDn_cn(2)*(cn_g(2)-cn_g(1))-(Dn_g(1)+Dn_g(2));
    Dv(id_cn+2:id_cn+Nn+1,id_cn+1:id_cn+Nn+2)=spdiags([d4*(-DDn_cn(1:Nn).*(cn_g(2:Nn+1)-cn_g(1:Nn))+Dn_g(1:Nn)+Dn_g(2:Nn+1)),...
        d4*(DDn_cn(2:Nn+1).*(cn_g(3:Nn+2)-cn_g(2:Nn+1))-(Dn_g(2:Nn+1)+Dn_g(3:Nn+2))-DDn_cn(2:Nn+1).*(cn_g(2:Nn+1)-cn_g(1:Nn))-(Dn_g(1:Nn)+Dn_g(2:Nn+1))),...
        d4*(DDn_cn(3:Nn+2).*(cn_g(3:Nn+2)-cn_g(2:Nn+1))+Dn_g(2:Nn+1)+Dn_g(3:Nn+2))],0:2,Nn,Nn+2);
    
    % Dcn_DTn
    Dv(id_cn+1:id_T+1)=(DDs_DTs(Ns+1)+DDs_DTs(Ns+2))*(cs_g(Ns+2)-cs_g(Ns+1))-(DDn_DTn(1)+DDn_DTn(2))*(cn_g(2)-cn_g(1));
    Dv(id_cn+2:id_cn+Nn+1,id_T+1)=d4*((DDn_DTn(2:Nn+1)+DDn_DTn(3:Nn+2)).*(cn_g(3:Nn+2)-cn_g(2:Nn+1))-(DDn_DTn(1:Nn)+DDn_DTn(2:Nn+1)).*(cn_g(2:Nn+1)-cn_g(1:Nn)));
    
    % Dap_Djp
    Dv(id_ap+1:id_ap+Np*Nr,id_jp+1:id_jp+Np)=sparse(Nr:Nr:Np*Nr,1:1:Np,-d43*((r_half(end)/r_whole(end))^2)*DapNL_jp(T,Tref,d44),Np*Nr,Np);
    
    % Dap_DTp
    Dv(id_ap+1:id_ap+Np*Nr,id_T+1)=kron(-d43*((r_half(end)/r_whole(end))^2)*DapNL_Tp(jp,T,Tref,d44),[zeros(Nr-1,1);1]);
    
    % Dan_Djn
    Dv(id_an+1:id_an+Nn*Nr,id_jn+1:id_jn+Nn)=sparse(Nr:Nr:Nn*Nr,1:1:Nn,-d45*((r_half(end)/r_whole(end))^2)*DanNL_jn(T,Tref,d46),Nn*Nr,Nn);
    
    % Dan_DTn
    Dv(id_an+1:id_an+Nn*Nr,id_T+1)=kron(-d45*((r_half(end)/r_whole(end))^2)*DanNL_Tn(jn,T,Tref,d46),[zeros(Nr-1,1);1]);
    
    % Dsp_Djp
    Dv(id_sp+1:id_sp+Np,id_jp+1:id_jp+Np)=spdiags(0.5*DapNL_jp(T*ones(Np,1),Tref,d44),0,Np,Np);
    
    % Dsp_DTp
    Dv(id_sp+1:id_sp+Np,id_T+1)=0.5*DapNL_Tp(jp,T,Tref,d44);
    
    % Dsn_Djn
    Dv(id_sn+1:id_sn+Nn,id_jn+1:id_jn+Nn)=spdiags(0.5*DanNL_jn(T*ones(Nn,1),Tref,d46),0,Nn,Nn);
    
    % Dsn_dTn
    Dv(id_sn+1:id_sn+Nn,id_T+1)=0.5*DanNL_Tn(jn,T,Tref,d46);
    
    % Djp_Dcp
    Dv(id_jp+1:id_jp+Np,2:Np+1)=spdiags(DjpNL_cp(cp_g(2:Np+1),sp,Phip,phip_g(2:Np+1),T,Tref),0,Np,Np);
    
    % Djp_Dsp
    Dv(id_jp+1:id_jp+Np,id_sp+1:id_sp+Np)=spdiags(DjpNL_sp(cp_g(2:Np+1),sp,Phip,phip_g(2:Np+1),T,Tref),0,Np,Np);
    
    % Djp_DPhip
    Dv(id_jp+1:id_jp+Np,id_Phip+1:id_Phip+Np)=spdiags(DjpNL_Phip(cp_g(2:Np+1),sp,Phip,phip_g(2:Np+1),T,Tref),0,Np,Np);
    
    % Djp_Dphip
    Dv(id_jp+1:id_jp+Np,id_phip+2:id_phip+Np+1)=spdiags(DjpNL_sphip(cp_g(2:Np+1),sp,Phip,phip_g(2:Np+1),T,Tref),0,Np,Np);
    
    % Djp_DTp
    Dv(id_jp+1:id_jp+Np,id_T+1)=DjpNL_Tp(cp_g(2:Np+1),sp,Phip,phip_g(2:Np+1),T,Tref);
    
    % Djn_Dcn
    Dv(id_jn+1:id_jn+Nn,id_cn+2:id_cn+Nn+1)=spdiags(DjnNL_cn(cn_g(2:Nn+1),sn,Phin,phin_g(2:Nn+1),T,Tref),0,Nn,Nn);
    
    % Djn_Dsn
    Dv(id_jn+1:id_jn+Nn,id_sn+1:id_sn+Nn)=spdiags(DjnNL_sn(cn_g(2:Nn+1),sn,Phin,phin_g(2:Nn+1),T,Tref),0,Nn,Nn);
    
    % Djn_DPhin
    Dv(id_jn+1:id_jn+Nn,id_Phin+1:id_Phin+Nn)=spdiags(DjnNL_Phin(cn_g(2:Nn+1),sn,Phin,phin_g(2:Nn+1),T,Tref),0,Nn,Nn);
    
    % Djn_Dphin
    Dv(id_jn+1:id_jn+Nn,id_phin+2:id_phin+Nn+1)=spdiags(DjnNL_sphin(cn_g(2:Nn+1),sn,Phin,phin_g(2:Nn+1),T,Tref),0,Nn,Nn);
    
    % Djn_DTn
    Dv(id_jn+1:id_jn+Nn,id_T+1)=DjnNL_Tn(cn_g(2:Nn+1),sn,Phin,phin_g(2:Nn+1),T,Tref);
    
    % Dphip_Dcp
    Dv(id_phip+2:id_phip+Np+1,1:Np+2)=spdiags([d15*(-Dkp_cp(1:Np).*(phip_g(2:Np+1)-phip_g(1:Np)))+d16*(-Dkp_cp(1:Np).*(2*T).*(Cp_g(2:Np+1)-Cp_g(1:Np))+1./cp_g(1:Np).*(kp_g(1:Np)+kp_g(2:Np+1)).*(2*T)),...
        d15*(Dkp_cp(2:Np+1).*(phip_g(3:Np+2)-phip_g(2:Np+1))-Dkp_cp(2:Np+1).*(phip_g(2:Np+1)-phip_g(1:Np)))+d16*(Dkp_cp(2:Np+1).*(2*T).*(Cp_g(3:Np+2)-Cp_g(2:Np+1))-1./cp_g(2:Np+1).*(kp_g(2:Np+1)+kp_g(3:Np+2)).*(2*T)-Dkp_cp(2:Np+1).*(2*T).*(Cp_g(2:Np+1)-Cp_g(1:Np))-1./cp_g(2:Np+1).*(kp_g(1:Np)+kp_g(2:Np+1)).*(2*T)),...
        d15*(Dkp_cp(3:Np+2).*(phip_g(3:Np+2)-phip_g(2:Np+1)))+d16*(Dkp_cp(3:Np+2).*(2*T).*(Cp_g(3:Np+2)-Cp_g(2:Np+1))+1./cp_g(3:Np+2).*(kp_g(2:Np+1)+kp_g(3:Np+2)).*(2*T))],0:2,Np,Np+2);
    Dv(id_phip+Np+2,Np+1)=Dkp_cp(Np+1)*(phip_g(Np+2)-phip_g(Np+1));
    Dv(id_phip+Np+2,Np+2)=Dkp_cp(Np+2)*(phip_g(Np+2)-phip_g(Np+1));
    Dv(id_phip+Np+2,id_cs+1)=-Dks_cs(1)*(phis_g(2)-phis_g(1));
    Dv(id_phip+Np+2,id_cs+2)=-Dks_cs(2)*(phis_g(2)-phis_g(1));
    
    % Dphis_Dcs
    Dv(id_phis+2:id_phis+Ns+1,id_cs+1:id_cs+Ns+2)=spdiags([-Dks_cs(1:Ns).*(phis_g(2:Ns+1)-phis_g(1:Ns))+d17*(-Dks_cs(1:Ns).*(2*T).*(Cs_g(2:Ns+1)-Cs_g(1:Ns))+1./cs_g(1:Ns).*(ks_g(1:Ns)+ks_g(2:Ns+1)).*(2*T)),...
        Dks_cs(2:Ns+1).*(phis_g(3:Ns+2)-phis_g(2:Ns+1))-Dks_cs(2:Ns+1).*(phis_g(2:Ns+1)-phis_g(1:Ns))+d17*(Dks_cs(2:Ns+1).*(2*T).*(Cs_g(3:Ns+2)-Cs_g(2:Ns+1))-1./cs_g(2:Ns+1).*(ks_g(2:Ns+1)+ks_g(3:Ns+2)).*(2*T)-Dks_cs(2:Ns+1).*(2*T).*(Cs_g(2:Ns+1)-Cs_g(1:Ns))-1./cs_g(2:Ns+1).*(ks_g(1:Ns)+ks_g(2:Ns+1)).*(2*T)),...
        Dks_cs(3:Ns+2).*(phis_g(3:Ns+2)-phis_g(2:Ns+1))+d17*(Dks_cs(3:Ns+2).*(2*T).*(Cs_g(3:Ns+2)-Cs_g(2:Ns+1))+1./cs_g(3:Ns+2).*(ks_g(2:Ns+1)+ks_g(3:Ns+2)).*(2*T))],0:2,Ns,Ns+2);
    
    % Dphin_Dcn
    Dv(id_phin+1,id_cs+Ns+1)=Dks_cs(Ns+1)*(phis_g(Ns+2)-phis_g(Ns+1)); 
    Dv(id_phin+1,id_cs+Ns+2)=Dks_cs(Ns+2)*(phis_g(Ns+2)-phis_g(Ns+1));
    Dv(id_phin+1,id_cn+1)=-Dkn_cn(1)*(phin_g(2)-phin_g(1));
    Dv(id_phin+1,id_cn+2)=-Dkn_cn(2)*(phin_g(2)-phin_g(1));
    Dv(id_phin+2:id_phin+Nn+1,id_cn+1:id_cn+Nn+2)=spdiags([d18*(-Dkn_cn(1:Nn).*(phin_g(2:Nn+1)-phin_g(1:Nn)))+d19*(-Dkn_cn(1:Nn).*(2*T).*(Cn_g(2:Nn+1)-Cn_g(1:Nn))+1./cn_g(1:Nn).*(kn_g(1:Nn)+kn_g(2:Nn+1)).*(2*T)),...
        d18*(Dkn_cn(2:Nn+1).*(phin_g(3:Nn+2)-phin_g(2:Nn+1))-Dkn_cn(2:Nn+1).*(phin_g(2:Nn+1)-phin_g(1:Nn)))+d19*(Dkn_cn(2:Nn+1).*(2*T).*(Cn_g(3:Nn+2)-Cn_g(2:Nn+1))-1./cn_g(2:Nn+1).*(kn_g(2:Nn+1)+kn_g(3:Nn+2)).*(2*T)-Dkn_cn(2:Nn+1).*(2*T).*(Cn_g(2:Nn+1)-Cn_g(1:Nn))-1./cn_g(2:Nn+1).*(kn_g(1:Nn)+kn_g(2:Nn+1)).*(2*T)),...
        d18*(Dkn_cn(3:Nn+2).*(phin_g(3:Nn+2)-phin_g(2:Nn+1)))+d19*(Dkn_cn(3:Nn+2).*(2*T).*(Cn_g(3:Nn+2)-Cn_g(2:Nn+1))+1./cn_g(3:Nn+2).*(kn_g(2:Nn+1)+kn_g(3:Nn+2)).*(2*T))],0:2,Nn,Nn+2);
    
    % Dphip_Dphip 
    Dv(id_phip+2:id_phip+Np+1,id_phip+1:id_phip+Np+2)=spdiags([d15*(kp_g(1:Np)+kp_g(2:Np+1)),...
        d15*(-(kp_g(2:Np+1)+kp_g(3:Np+2))-(kp_g(1:Np)+kp_g(2:Np+1))),...
        d15*(kp_g(2:Np+1)+kp_g(3:Np+2))],0:2,Np,Np+2);
    Dv(id_phip+Np+2,id_phip+Np+1)=-(kp_g(Np+1)+kp_g(Np+2));
    Dv(id_phip+Np+2,id_phip+Np+2)=kp_g(Np+1)+kp_g(Np+2);
    Dv(id_phip+Np+2,id_phis+1)=ks_g(1)+ks_g(2);
    Dv(id_phip+Np+2,id_phis+2)=-(ks_g(1)+ks_g(2));
    
    % Dphis_Dphis
    Dv(id_phis+2:id_phis+Ns+1,id_phis+1:id_phis+Ns+2)=spdiags([ks_g(1:Ns)+ks_g(2:Ns+1),...
        -(ks_g(2:Ns+1)+ks_g(3:Ns+2))-(ks_g(1:Ns)+ks_g(2:Ns+1)),...
        ks_g(2:Ns+1)+ks_g(3:Ns+2)],0:2,Ns,Ns+2);
    
    % Dphin_Dphin
    Dv(id_phin+1,id_phis+Ns+1)=-(ks_g(Ns+1)+ks_g(Ns+2));
    Dv(id_phin+1,id_phis+Ns+2)=ks_g(Ns+1)+ks_g(Ns+2);
    Dv(id_phin+1,id_phin+1)=kn_g(1)+kn_g(2);
    Dv(id_phin+1,id_phin+2)=-(kn_g(1)+kn_g(2));
    Dv(id_phin+2:id_phin+Nn+1,id_phin+1:id_phin+Nn+2)=spdiags([d18*(kn_g(1:Nn)+kn_g(2:Nn+1)),...
        d18*(-(kn_g(2:Nn+1)+kn_g(3:Nn+2))-(kn_g(1:Nn)+kn_g(2:Nn+1))),...
        d18*(kn_g(2:Nn+1)+kn_g(3:Nn+2))],0:2,Nn,Nn+2);
    
    % Dphip_DTp
    Dv(id_phip+2:id_phip+Np+1,id_T+1)=d15*((Dkp_Tp(2:Np+1)+Dkp_Tp(3:Np+2)).*(phip_g(3:Np+2)-phip_g(2:Np+1))-(Dkp_Tp(1:Np)+Dkp_Tp(2:Np+1)).*(phip_g(2:Np+1)-phip_g(1:Np)))+...
        d16*((Cp_g(3:Np+2)-Cp_g(2:Np+1)).*(2*(kp_g(2:Np+1)+kp_g(3:Np+2))+2*T*(Dkp_Tp(2:Np+1)+Dkp_Tp(3:Np+2)))-(Cp_g(2:Np+1)-Cp_g(1:Np)).*(2*(kp_g(1:Np)+kp_g(2:Np+1))+2*T*(Dkp_Tp(1:Np)+Dkp_Tp(2:Np+1))));
    Dv(id_phip+Np+2,id_T+1)=(Dkp_Tp(Np+1)+Dkp_Tp(Np+2))*(phip_g(Np+2)-phip_g(Np+1))-(Dks_Ts(1)+Dks_Ts(2))*(phis_g(2)-phis_g(1));

    % Dphis_DTs
    Dv(id_phis+2:id_phis+Ns+1,id_T+1)=(Dks_Ts(2:Ns+1)+Dks_Ts(3:Ns+2)).*(phis_g(3:Ns+2)-phis_g(2:Ns+1))-(Dks_Ts(1:Ns)+Dks_Ts(2:Ns+1)).*(phis_g(2:Ns+1)-phis_g(1:Ns))+...
        d17*((Cs_g(3:Ns+2)-Cs_g(2:Ns+1)).*(2*(ks_g(2:Ns+1)+ks_g(3:Ns+2))+2*T*(Dks_Ts(2:Ns+1)+Dks_Ts(3:Ns+2)))-(Cs_g(3:Ns+2)-Cs_g(2:Ns+1)).*(2*(ks_g(1:Ns)+ks_g(2:Ns+1))+2*T*(Dks_Ts(1:Ns)+Dks_Ts(2:Ns+1))));

    % Dphin_DTn
    Dv(id_phin+1,id_T+1)=(Dks_Ts(Ns+1)+Dks_Ts(Ns+2))*(phis_g(Ns+2)-phis_g(Ns+1))-(Dkn_Tn(1)+Dkn_Tn(2))*(phin_g(2)-phin_g(1));
    Dv(id_phin+2:id_phin+Nn+1,id_T+1)=d18*((Dkn_Tn(2:Nn+1)+Dkn_Tn(3:Nn+2)).*(phin_g(3:Nn+2)-phin_g(2:Nn+1))-(Dkn_Tn(1:Nn)+Dkn_Tn(2:Nn+1)).*(phin_g(2:Nn+1)-phin_g(1:Nn)))+...
        d19*((Cn_g(3:Nn+2)-Cn_g(2:Nn+1)).*(2*(kn_g(2:Nn+1)+kn_g(3:Nn+2))+2*T*(Dkn_Tn(2:Nn+1)+Dkn_Tn(3:Nn+2)))-(Cn_g(2:Nn+1)-Cn_g(1:Nn)).*(2*(kn_g(1:Nn)+kn_g(2:Nn+1))+2*T*(Dkn_Tn(1:Nn)+Dkn_Tn(2:Nn+1))));
    
    % DTp_Dcp
    Dv(id_T+1,1:Np+2)=sum(spdiags([-d28./cp_g(1:Np).*kp_g(2:Np+1).*T.*(phip_g(3:Np+2)-phip_g(1:Np)),...
        d27*Dkp_cp(2:Np+1).*((phip_g(3:Np+2)-phip_g(1:Np)).^2)+d28*Dkp_cp(2:Np+1).*T.*(Cp_g(3:Np+2)-Cp_g(1:Np)).*(phip_g(3:Np+2)-phip_g(1:Np)),...
        d28./cp_g(3:Np+2).*kp_g(2:Np+1).*T.*(phip_g(3:Np+2)-phip_g(1:Np))],0:2,Np,Np+2));
    
    % DTs_Dcs
    Dv(id_T+1,id_cs+1:id_cs+Ns+2)=sum(spdiags([-d32./cs_g(1:Ns).*ks_g(2:Ns+1).*T.*(phis_g(3:Ns+2)-phis_g(1:Ns)),...
        d31*Dks_cs(2:Ns+1).*((phis_g(3:Ns+2)-phis_g(1:Ns)).^2)+d32*Dks_cs(2:Ns+1).*T.*(Cs_g(3:Ns+2)-Cs_g(1:Ns)).*(phis_g(3:Ns+2)-phis_g(1:Ns)),...
        d32./cs_g(3:Ns+2).*ks_g(2:Ns+1).*T.*(phis_g(3:Ns+2)-phis_g(1:Ns))],0:2,Ns,Ns+2));
    
    % DTn_Dcn
    Dv(id_T+1,id_cn+1:id_cn+Nn+2)=sum(spdiags([-d36./cn_g(1:Nn).*kn_g(2:Nn+1).*T.*(phin_g(3:Nn+2)-phin_g(1:Nn)),...
        d35*Dkn_cn(2:Nn+1).*((phin_g(3:Nn+2)-phin_g(1:Nn)).^2)+d36*Dkn_cn(2:Nn+1).*T.*(Cn_g(3:Nn+2)-Cn_g(1:Nn)).*(phin_g(3:Nn+2)-phin_g(1:Nn)),...
        d36./cn_g(3:Nn+2).*kn_g(2:Nn+1).*T.*(phin_g(3:Nn+2)-phin_g(1:Nn))],0:2,Nn,Nn+2));
    
    % DTp_Dsp
    Dv(id_T+1,id_sp+1:id_sp+Np)=sum(spdiags(DQp_sp(sp,jp,T,Tref,d29),0,Np,Np));
    
    % DTn_Dsn
    Dv(id_T+1,id_sn+1:id_sn+Nn)=sum(spdiags(DQn_sn(sn,jn,T,Tref,d37),0,Nn,Nn));
    
    % DTp_Djp
    Dv(id_T+1,id_jp+1:id_jp+Np)=sum(spdiags(DQp_jp(sp,Phip,phip_g(2:Np+1),T,Tref,d29),0,Np,Np));
    
    % DTn_Djn
    Dv(id_T+1,id_jn+1:id_jn+Nn)=sum(spdiags(DQn_jn(sn,Phin,phin_g(2:Nn+1),T,Tref,d37),0,Nn,Nn));
    
    % DTp_DPhip
    DTp_DPhip=zeros(Np);
    bottom_TpPhip=[-2*d26*(Phip(3:Np)-Phip(1:Np-2));-2*d26*(Phip(Np)-Phip(Np-1))];
    top_TpPhip=[0;2*d26*(Phip(2)-Phip(1)-h*I/59);2*d26*(Phip(3:Np)-Phip(1:Np-2))];
    DTp_DPhip(1:Np,1:Np)=spdiags(bottom_TpPhip,-1,Np,Np)+spdiags(DQp_Phip(jp,d29),0,Np,Np)+spdiags(top_TpPhip,1,Np,Np);
    DTp_DPhip(1,1)=DTp_DPhip(1,1)-2*d26*(Phip(2)-Phip(1)-h*I/59);
    DTp_DPhip(Np,Np)=DTp_DPhip(Np,Np)+2*d26*(Phip(Np)-Phip(Np-2));
    Dv(id_T+1,id_Phip+1:id_Phip+Np)=sum(DTp_DPhip);
    
    % DTn_DPhin
    DTn_DPhin=zeros(Nn);
    bottom_TnPhin=[-2*d34*(Phin(3:Nn)-Phin(1:Nn-2));-2*d34*(Phin(Nn)-h*I/48.24-Phin(Nn-1))];
    top_TnPhin=[0;2*d34*(Phin(2)-Phin(1));2*d34*(Phin(3:Nn)-Phin(1:Nn-2))];
    DTn_DPhin(1:Nn,1:Nn)=spdiags(bottom_TnPhin,-1,Nn,Nn)+spdiags(DQn_Tn(sn,jn,d37),0,Nn,Nn)+spdiags(top_TnPhin,1,Nn,Nn);
    DTn_DPhin(1,1)=DTn_DPhin(1,1)-2*d34*(Phin(2)-Phin(1));
    DTn_DPhin(Nn,Nn)=DTn_DPhin(Nn,Nn)-h*I/48.24-Phin(Nn-1);
    Dv(id_T+1,id_Phin+1:id_Phin+Nn)=sum(DTn_DPhin);
    
    % DTp_Dphip
    Dv(id_T+1,id_phip+1:id_phip+Np+2)=sum(spdiags([-2*d27*kp_g(2:Np+1).*(phip_g(3:Np+2)-phip_g(1:Np))-d28*kp_g(2:Np+1).*T.*(Cp_g(3:Np+2)-Cp_g(1:Np)),...
        DQp_sphip(jp,d29),2*d27*kp_g(2:Np+1).*(phip_g(3:Np+2)-phip_g(1:Np))+d28*kp_g(2:Np+1).*T.*(Cp_g(3:Np+2)-Cp_g(1:Np))],0:2,Np,Np+2));
    
    % DTs_Dphis
    Dv(id_T+1,id_phis+1:id_phis+Ns+2)=sum(spdiags([-2*d31*ks_g(2:Ns+1).*(phis_g(3:Ns+2)-phis_g(1:Ns))-d32*ks_g(2:Ns+1).*T.*(Cs_g(3:Ns+2)-Cs_g(1:Ns)),...
        zeros(Ns,1),2*d31*ks_g(2:Ns+1).*(phis_g(3:Ns+2)-phis_g(1:Ns))+d32*ks_g(2:Ns+1).*T.*(Cs_g(3:Ns+2)-Cs_g(1:Ns))],0:2,Ns,Ns+2));
    
    % DTn_Dphin
    Dv(id_T+1,id_phin+1:id_phin+Nn+2)=sum(spdiags([-2*d35*kn_g(2:Nn+1).*(phin_g(3:Nn+2)-phin_g(1:Nn))-d36*kn_g(2:Nn+1).*T.*(Cn_g(3:Nn+2)-Cn_g(1:Nn)),...
        DQn_sphin(jn,d37),2*d35*kn_g(2:Nn+1).*(phin_g(3:Nn+2)-phin_g(1:Nn))+d36*kn_g(2:Nn+1).*T.*(Cn_g(3:Nn+2)-Cn_g(1:Nn))],0:2,Nn,Nn+2));
    
    % DT_DT
    Dv(id_T+1,id_T+1)=sum(d27*Dkp_Tp(2:Np+1).*((phip_g(3:Np+2)-phip_g(1:Np)).^2)+d28*(Cp_g(3:Np+2)-Cp_g(1:Np)).*(phip_g(3:Np+2)-phip_g(1:Np)).*(Dkp_Tp(2:Np+1).*T+kp_g(2:Np+1))+DQp_Tp(sp,jp,d29))+...
        sum(d31*Dks_Ts(2:Ns+1).*((phis_g(3:Ns+2)-phis_g(1:Ns)).^2)+d32*(Cs_g(3:Ns+2)-Cs_g(1:Ns)).*(phis_g(3:Ns+2)-phis_g(1:Ns)).*(Dks_Ts(2:Ns+1).*T+ks_g(2:Ns+1)))+...
        sum(d35*Dkn_Tn(2:Nn+1).*((phin_g(3:Nn+2)-phin_g(1:Nn)).^2)+d36*(Cn_g(3:Nn+2)-Cn_g(1:Nn)).*(phin_g(3:Nn+2)-phin_g(1:Nn)).*(Dkn_Tn(2:Nn+1).*T+kn_g(2:Nn+1))+DQn_Tn(sn,jn,d37));
end