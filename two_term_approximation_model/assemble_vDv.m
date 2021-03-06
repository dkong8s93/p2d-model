function [v,Dv]=assemble_vDv(d1,d3,d4,d15,d16,d17,d18,d19,d21,d24,d26,d27,d28,d29,d31,d32,d34,d35,d36,d37,d39,d42,h,Tref,I,U,Na,Np,Ns,Nn,Nz)
    dim=Na+7*Np+3*Ns+7*Nn+Nz+2*11; 
    
    % indices for convenience
    [id_cs,id_cn,~,~,id_sp,id_sn,id_jp,id_jn,id_Phip,id_Phin,id_phip,id_phis,id_phin,id_Ta,id_Tp,id_Ts,id_Tn,id_Tz]=indices(Na,Np,Ns,Nn);
    
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
    Ta_g=U(id_Ta+1:id_Ta+Na+2);
    Tp_g=U(id_Tp+1:id_Tp+Np+2);
    Ts_g=U(id_Ts+1:id_Ts+Ns+2);
    Tn_g=U(id_Tn+1:id_Tn+Nn+2);
    Tz_g=U(id_Tz+1:id_Tz+Nz+2);
    
    % Auxiliary functions
    Dp_g=Deffp(cp_g,Tp_g); % Diffusion coefficient for the cathode check
    Ds_g=Deffs(cs_g,Ts_g); % Diffusion coefficient for the separator check
    Dn_g=Deffn(cn_g,Tn_g); % Diffusion coefficient for the anode check
    kp_g=kappap(cp_g,Tp_g); % Conductivity for the cathode check
    ks_g=kappas(cs_g,Ts_g); % Conductivity for the separator check
    kn_g=kappan(cn_g,Tn_g); % Conductivity for the anode check
    DDp_cp=DDeffp_cp(cp_g,Tp_g); 
    DDs_cs=DDeffs_cs(cs_g,Ts_g);
    DDn_cn=DDeffn_cn(cn_g,Tn_g);
    DDp_DTp=DDeffp_Tp(cp_g,Tp_g);
    DDs_DTs=DDeffs_Ts(cs_g,Ts_g); 
    DDn_DTn=DDeffn_Tn(cn_g,Tn_g);
    Dkp_cp=Dkappap_cp(cp_g,Tp_g);
    Dks_cs=Dkappas_cs(cs_g,Ts_g);
    Dkn_cn=Dkappan_cn(cn_g,Tn_g);
    Dkp_Tp=Dkappap_Tp(cp_g,Tp_g);
    Dks_Ts=Dkappas_Ts(cs_g,Ts_g);
    Dkn_Tn=Dkappan_Tn(cn_g,Tn_g);

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
    % sp
    v(id_sp+1:id_sp+Np)=spNL(Tp_g(2:Np+1),Tref).*jp;
    % sn
    v(id_sn+1:id_sn+Nn)=snNL(Tn_g(2:Nn+1),Tref).*jn;
    % jp
    v(id_jp+1:id_jp+Np)=jpNL(cp_g(2:Np+1),sp,Phip,phip_g(2:Np+1),Tp_g(2:Np+1),Tref);
    % jn
    v(id_jn+1:id_jn+Nn)=jnNL(cn_g(2:Nn+1),sn,Phin,phin_g(2:Nn+1),Tn_g(2:Nn+1),Tref);
    % Phip
    v(id_Phip+1)=h*I/59;  % check
    % Phin
    v(id_Phin+Nn)=-h*I/48.24; % check
    % phip
    v(id_phip+2:id_phip+Np+1)=d15*((kp_g(2:Np+1)+kp_g(3:Np+2)).*(phip_g(3:Np+2)-phip_g(2:Np+1))-(kp_g(1:Np)+kp_g(2:Np+1)).*(phip_g(2:Np+1)-phip_g(1:Np)))+...
        d16*((kp_g(2:Np+1)+kp_g(3:Np+2)).*(Tp_g(2:Np+1)+Tp_g(3:Np+2)).*(Cp_g(3:Np+2)-Cp_g(2:Np+1))-(kp_g(1:Np)+kp_g(2:Np+1)).*(Tp_g(1:Np)+Tp_g(2:Np+1)).*(Cp_g(2:Np+1)-Cp_g(1:Np)));
    v(id_phip+Np+2)=(kp_g(Np+1)+kp_g(Np+2))*(phip_g(Np+2)-phip_g(Np+1))-(ks_g(1)+ks_g(2))*(phis_g(2)-phis_g(1)); % check
    % phis
    v(id_phis+2:id_phis+Ns+1)=(ks_g(2:Ns+1)+ks_g(3:Ns+2)).*(phis_g(3:Ns+2)-phis_g(2:Ns+1))-(ks_g(1:Ns)+ks_g(2:Ns+1)).*(phis_g(2:Ns+1)-phis_g(1:Ns))+...
        d17*((ks_g(2:Ns+1)+ks_g(3:Ns+2)).*(Ts_g(2:Ns+1)+Ts_g(3:Ns+2)).*(Cs_g(3:Ns+2)-Cs_g(2:Ns+1))-(ks_g(1:Ns)+ks_g(2:Ns+1)).*(Ts_g(1:Ns)+Ts_g(2:Ns+1)).*(Cs_g(2:Ns+1)-Cs_g(1:Ns)));
    % phin
    v(id_phin+1)=(ks_g(Ns+1)+ks_g(Ns+2))*(phis_g(Ns+2)-phis_g(Ns+1))-(kn_g(1)+kn_g(2))*(phin_g(2)-phin_g(1)); % check
    v(id_phin+2:id_phin+Nn+1)=d18*((kn_g(2:Nn+1)+kn_g(3:Nn+2)).*(phin_g(3:Nn+2)-phin_g(2:Nn+1))-(kn_g(1:Nn)+kn_g(2:Nn+1)).*(phin_g(2:Nn+1)-phin_g(1:Nn)))+...
        d19*((kn_g(2:Nn+1)+kn_g(3:Nn+2)).*(Tn_g(2:Nn+1)+Tn_g(3:Nn+2)).*(Cn_g(3:Nn+2)-Cn_g(2:Nn+1))-(kn_g(1:Nn)+kn_g(2:Nn+1)).*(Tn_g(1:Nn)+Tn_g(2:Nn+1)).*(Cn_g(2:Nn+1)-Cn_g(1:Nn)));
    % Ta_g
    v(id_Ta+1)=d24;
    v(id_Ta+2:id_Ta+Na+1)=d21;
    % Tp_g
    v(id_Tp+2)=d26*((Phip(2)-Phip(1)-h*I/59)^2)+d27*kp_g(2)*((phip_g(3)-phip_g(1))^2)+d28*kp_g(2)*Tp_g(2)*(Cp_g(3)-Cp_g(1))*(phip_g(3)-phip_g(1))+Qp(sp(1),jp(1),Phip(1),phip_g(2),Tp_g(2),Tref,d29);
    v(id_Tp+3:id_Tp+Np)=d26*((Phip(3:Np)-Phip(1:Np-2)).^2)+d27*kp_g(3:Np).*((phip_g(4:Np+1)-phip_g(2:Np-1)).^2)+d28*kp_g(3:Np).*Tp_g(3:Np).*(Cp_g(4:Np+1)-Cp_g(2:Np-1)).*(phip_g(4:Np+1)-phip_g(2:Np-1))+Qp(sp(2:Np-1),jp(2:Np-1),Phip(2:Np-1),phip_g(3:Np),Tp_g(3:Np),Tref,d29);
    v(id_Tp+Np+1)=d26*((Phip(Np)-Phip(Np-1))^2)+d27*kp_g(Np+1)*((phip_g(Np+2)-phip_g(Np))^2)+d28*kp_g(Np+1)*Tp_g(Np+1)*(Cp_g(Np+2)-Cp_g(Np))*(phip_g(Np+2)-phip_g(Np))+Qp(sp(Np),jp(Np),Phip(Np),phip_g(Np+1),Tp_g(Np+1),Tref,d29);
    % Ts_g
    v(id_Ts+2:id_Ts+Ns+1)=d31*ks_g(2:Ns+1).*((phis_g(3:Ns+2)-phis_g(1:Ns)).^2)+d32*ks_g(2:Ns+1).*Ts_g(2:Ns+1).*(Cs_g(3:Ns+2)-Cs_g(1:Ns)).*(phis_g(3:Ns+2)-phis_g(1:Ns));
    % Tn_g
    v(id_Tn+2)=d34*((Phin(2)-Phin(1))^2)+d35*kn_g(2)*((phin_g(3)-phin_g(1))^2)+d36*kn_g(2)*Tn_g(2)*(Cn_g(3)-Cn_g(1))*(phin_g(3)-phin_g(1))+Qn(sn(1),jn(1),Phin(1),phin_g(2),Tn_g(2),Tref,d37);
    v(id_Tn+3:id_Tn+Nn)=d34*((Phin(3:Nn)-Phin(1:Nn-2)).^2)+d35*kn_g(3:Nn).*((phin_g(4:Nn+1)-phin_g(2:Nn-1)).^2)+d36*kn_g(3:Nn).*Tn_g(3:Nn).*(Cn_g(4:Nn+1)-Cn_g(2:Nn-1)).*(phin_g(4:Nn+1)-phin_g(2:Nn-1))+Qn(sn(2:Nn-1),jn(2:Nn-1),Phin(2:Nn-1),phin_g(3:Nn),Tn_g(3:Nn),Tref,d37);
    v(id_Tn+Nn+1)=d34*(Phin(Nn)-h*I/48.24-Phin(Nn-1))^2+d35*kn_g(Nn+1)*((phin_g(Nn+2)-phin_g(Nn))^2)+d36*kn_g(Nn+1)*Tn_g(Nn+1)*(Cn_g(Nn+2)-Cn_g(Nn))*(phin_g(Nn+2)-phin_g(Nn))+Qn(sn(Nn),jn(Nn),Phin(Nn),phin_g(Nn+1),Tn_g(Nn+1),Tref,d37);
    % Tz_g
    v(id_Tz+2:id_Tz+Nz+1)=d39; 
    v(id_Tz+Nz+2)=d42;
    
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
    Dv(2:Np+1,id_Tp+1:id_Tp+Np+2)=spdiags([d1*(-DDp_DTp(1:Np).*(cp_g(2:Np+1)-cp_g(1:Np))),d1*(DDp_DTp(2:Np+1).*(cp_g(3:Np+2)-cp_g(2:Np+1))-DDp_DTp(2:Np+1).*(cp_g(2:Np+1)-cp_g(1:Np))),d1*(DDp_DTp(3:Np+2).*(cp_g(3:Np+2)-cp_g(2:Np+1)))],0:2,Np,Np+2);
    Dv(Np+2,id_Tp+Np+1)=DDp_DTp(Np+1)*(cp_g(Np+2)-cp_g(Np+1));
    Dv(Np+2,id_Tp+Np+2)=DDp_DTp(Np+2)*(cp_g(Np+2)-cp_g(Np+1));
    Dv(Np+2,id_Ts+1)=-DDs_DTs(1)*(cs_g(2)-cs_g(1));
    Dv(Np+2,id_Ts+2)=-DDs_DTs(2)*(cs_g(2)-cs_g(1));
    % Dcs_Dcs
    Dv(id_cs+2:id_cs+Ns+1,id_cs+1:id_cs+Ns+2)=spdiags([d3*(-DDs_cs(1:Ns).*(cs_g(2:Ns+1)-cs_g(1:Ns))+Ds_g(1:Ns)+Ds_g(2:Ns+1)),...
        d3*(DDs_cs(2:Ns+1).*(cs_g(3:Ns+2)-cs_g(2:Ns+1))-(Ds_g(2:Ns+1)+Ds_g(3:Ns+2))-DDs_cs(2:Ns+1).*(cs_g(2:Ns+1)-cs_g(1:Ns))-(Ds_g(1:Ns)+Ds_g(2:Ns+1)))...
        d3*(DDs_cs(3:Ns+2).*(cs_g(3:Ns+2)-cs_g(2:Ns+1))+Ds_g(2:Ns+1)+Ds_g(3:Ns+2))],0:2,Ns,Ns+2);
    % Dcs_DTs
    Dv(id_cs+2:id_cs+Ns+1,id_Ts+1:id_Ts+Ns+2)=spdiags([d3*(-DDs_DTs(1:Ns).*(cs_g(2:Ns+1)-cs_g(1:Ns))),...
        d3*(DDs_DTs(2:Ns+1).*(cs_g(3:Ns+2)-cs_g(2:Ns+1))),d3*(DDs_DTs(3:Ns+2).*(cs_g(3:Ns+2)-cs_g(2:Ns+1)))],0:2,Ns,Ns+2);
    % Dcn_Dcn
    Dv(id_cn+1,id_cs+Ns+1)=DDs_cs(Ns+1)*(cs_g(Ns+2)-cs_g(Ns+1))-(Ds_g(Ns+1)+Ds_g(Ns+2));
    Dv(id_cn+1,id_cs+Ns+2)=DDs_cs(Ns+2)*(cs_g(Ns+2)-cs_g(Ns+1))+Ds_g(Ns+1)+Ds_g(Ns+2);
    Dv(id_cn+1,id_cn+1)=-DDn_cn(1)*(cn_g(2)-cn_g(1))+Dn_g(1)+Dn_g(2);
    Dv(id_cn+1,id_cn+2)=-DDn_cn(2)*(cn_g(2)-cn_g(1))-(Dn_g(1)+Dn_g(2));
    Dv(id_cn+2:id_cn+Nn+1,id_cn+1:id_cn+Nn+2)=spdiags([d4*(-DDn_cn(1:Nn).*(cn_g(2:Nn+1)-cn_g(1:Nn))+Dn_g(1:Nn)+Dn_g(2:Nn+1)),...
        d4*(DDn_cn(2:Nn+1).*(cn_g(3:Nn+2)-cn_g(2:Nn+1))-(Dn_g(2:Nn+1)+Dn_g(3:Nn+2))-DDn_cn(2:Nn+1).*(cn_g(2:Nn+1)-cn_g(1:Nn))-(Dn_g(1:Nn)+Dn_g(2:Nn+1))),...
        d4*(DDn_cn(3:Nn+2).*(cn_g(3:Nn+2)-cn_g(2:Nn+1))+Dn_g(2:Nn+1)+Dn_g(3:Nn+2))],0:2,Nn,Nn+2);
    % Dcn_DTn
    Dv(id_cn+1,id_Ts+Ns+1)=DDs_DTs(Ns+1)*(cs_g(Ns+2)-cs_g(Ns+1));
    Dv(id_cn+1,id_Ts+Ns+2)=DDs_DTs(Ns+2)*(cs_g(Ns+2)-cs_g(Ns+1));
    Dv(id_cn+1,id_Tn+1)=-DDn_DTn(1)*(cn_g(2)-cn_g(1));
    Dv(id_cn+1,id_Tn+2)=-DDn_DTn(2)*(cn_g(2)-cn_g(1));
    Dv(id_cn+2:id_cn+Nn+1,id_Tn+1:id_Tn+Nn+2)=spdiags([d4*(-DDn_DTn(1:Nn).*(cn_g(2:Nn+1)-cn_g(1:Nn))),...
        d4*(DDn_DTn(2:Nn+1).*(cn_g(3:Nn+2)-cn_g(2:Nn+1))-DDn_DTn(2:Nn+1).*(cn_g(2:Nn+1)-cn_g(1:Nn))),...
        d4*(DDn_DTn(3:Nn+2).*(cn_g(3:Nn+2)-cn_g(2:Nn+1)))],0:2,Nn,Nn+2);
    % Dsp_Djp
    Dv(id_sp+1:id_sp+Np,id_jp+1:id_jp+Np)=spdiags(spNL(Tp_g(2:Np+1),Tref),0,Np,Np);
    % Dsp_DTp
    Dv(id_sp+1:id_sp+Np,id_Tp+2:id_Tp+Np+1)=spdiags(DspNL_Tp(Tp_g(2:Np+1),Tref).*jp,0,Np,Np);
    % Dsn_Djn
    Dv(id_sn+1:id_sn+Nn,id_jn+1:id_jn+Nn)=spdiags(snNL(Tn_g(2:Nn+1),Tref),0,Nn,Nn);
    % Dsn_dTn
    Dv(id_sn+1:id_sn+Nn,id_Tn+2:id_Tn+Nn+1)=spdiags(DsnNL_Tn(Tn_g(2:Nn+1),Tref).*jn,0,Nn,Nn);
    % Djp_Dcp
    Dv(id_jp+1:id_jp+Np,2:Np+1)=spdiags(DjpNL_cp(cp_g(2:Np+1),sp,Phip,phip_g(2:Np+1),Tp_g(2:Np+1),Tref),0,Np,Np);
    % Djp_Dsp
    Dv(id_jp+1:id_jp+Np,id_sp+1:id_sp+Np)=spdiags(DjpNL_sp(cp_g(2:Np+1),sp,Phip,phip_g(2:Np+1),Tp_g(2:Np+1),Tref),0,Np,Np);
    % Djp_DPhip
    Dv(id_jp+1:id_jp+Np,id_Phip+1:id_Phip+Np)=spdiags(DjpNL_Phip(cp_g(2:Np+1),sp,Phip,phip_g(2:Np+1),Tp_g(2:Np+1),Tref),0,Np,Np);
    % Djp_Dphip
    Dv(id_jp+1:id_jp+Np,id_phip+2:id_phip+Np+1)=spdiags(DjpNL_sphip(cp_g(2:Np+1),sp,Phip,phip_g(2:Np+1),Tp_g(2:Np+1),Tref),0,Np,Np);
    % Djp_DTp
    Dv(id_jp+1:id_jp+Np,id_Tp+2:id_Tp+Np+1)=spdiags(DjpNL_Tp(cp_g(2:Np+1),sp,Phip,phip_g(2:Np+1),Tp_g(2:Np+1),Tref),0,Np,Np);
    % Djn_Dcn
    Dv(id_jn+1:id_jn+Nn,id_cn+2:id_cn+Nn+1)=spdiags(DjnNL_cn(cn_g(2:Nn+1),sn,Phin,phin_g(2:Nn+1),Tn_g(2:Nn+1),Tref),0,Nn,Nn);
    % Djn_Dsn
    Dv(id_jn+1:id_jn+Nn,id_sn+1:id_sn+Nn)=spdiags(DjnNL_sn(cn_g(2:Nn+1),sn,Phin,phin_g(2:Nn+1),Tn_g(2:Nn+1),Tref),0,Nn,Nn);
    % Djn_DPhin
    Dv(id_jn+1:id_jn+Nn,id_Phin+1:id_Phin+Nn)=spdiags(DjnNL_Phin(cn_g(2:Nn+1),sn,Phin,phin_g(2:Nn+1),Tn_g(2:Nn+1),Tref),0,Nn,Nn);
    % Djn_Dphin
    Dv(id_jn+1:id_jn+Nn,id_phin+2:id_phin+Nn+1)=spdiags(DjnNL_sphin(cn_g(2:Nn+1),sn,Phin,phin_g(2:Nn+1),Tn_g(2:Nn+1),Tref),0,Nn,Nn);
    % Djn_DTn
    Dv(id_jn+1:id_jn+Nn,id_Tn+2:id_Tn+Nn+1)=spdiags(DjnNL_Tn(cn_g(2:Nn+1),sn,Phin,phin_g(2:Nn+1),Tn_g(2:Nn+1),Tref),0,Nn,Nn);
    % Dphip_Dcp
    Dv(id_phip+2:id_phip+Np+1,1:Np+2)=spdiags([d15*(-Dkp_cp(1:Np).*(phip_g(2:Np+1)-phip_g(1:Np)))+d16*(-Dkp_cp(1:Np).*(Tp_g(1:Np)+Tp_g(2:Np+1)).*(Cp_g(2:Np+1)-Cp_g(1:Np))+1./cp_g(1:Np).*(kp_g(1:Np)+kp_g(2:Np+1)).*(Tp_g(1:Np)+Tp_g(2:Np+1))),...
        d15*(Dkp_cp(2:Np+1).*(phip_g(3:Np+2)-phip_g(2:Np+1))-Dkp_cp(2:Np+1).*(phip_g(2:Np+1)-phip_g(1:Np)))+d16*(Dkp_cp(2:Np+1).*(Tp_g(2:Np+1)+Tp_g(3:Np+2)).*(Cp_g(3:Np+2)-Cp_g(2:Np+1))-1./cp_g(2:Np+1).*(kp_g(2:Np+1)+kp_g(3:Np+2)).*(Tp_g(2:Np+1)+Tp_g(3:Np+2))-Dkp_cp(2:Np+1).*(Tp_g(1:Np)+Tp_g(2:Np+1)).*(Cp_g(2:Np+1)-Cp_g(1:Np))-1./cp_g(2:Np+1).*(kp_g(1:Np)+kp_g(2:Np+1)).*(Tp_g(1:Np)+Tp_g(2:Np+1))),...
        d15*(Dkp_cp(3:Np+2).*(phip_g(3:Np+2)-phip_g(2:Np+1)))+d16*(Dkp_cp(3:Np+2).*(Tp_g(2:Np+1)+Tp_g(3:Np+2)).*(Cp_g(3:Np+2)-Cp_g(2:Np+1))+1./cp_g(3:Np+2).*(kp_g(2:Np+1)+kp_g(3:Np+2)).*(Tp_g(2:Np+1)+Tp_g(3:Np+2)))],0:2,Np,Np+2);
    Dv(id_phip+Np+2,Np+1)=Dkp_cp(Np+1)*(phip_g(Np+2)-phip_g(Np+1));
    Dv(id_phip+Np+2,Np+2)=Dkp_cp(Np+2)*(phip_g(Np+2)-phip_g(Np+1));
    Dv(id_phip+Np+2,id_cs+1)=-Dks_cs(1)*(phis_g(2)-phis_g(1));
    Dv(id_phip+Np+2,id_cs+2)=-Dks_cs(2)*(phis_g(2)-phis_g(1));
    % Dphis_Dcs
    Dv(id_phis+2:id_phis+Ns+1,id_cs+1:id_cs+Ns+2)=spdiags([-Dks_cs(1:Ns).*(phis_g(2:Ns+1)-phis_g(1:Ns))+d17*(-Dks_cs(1:Ns).*(Ts_g(1:Ns)+Ts_g(2:Ns+1)).*(Cs_g(2:Ns+1)-Cs_g(1:Ns))+1./cs_g(1:Ns).*(ks_g(1:Ns)+ks_g(2:Ns+1)).*(Ts_g(1:Ns)+Ts_g(2:Ns+1))),...
        Dks_cs(2:Ns+1).*(phis_g(3:Ns+2)-phis_g(2:Ns+1))-Dks_cs(2:Ns+1).*(phis_g(2:Ns+1)-phis_g(1:Ns))+d17*(Dks_cs(2:Ns+1).*(Ts_g(2:Ns+1)+Ts_g(3:Ns+2)).*(Cs_g(3:Ns+2)-Cs_g(2:Ns+1))-1./cs_g(2:Ns+1).*(ks_g(2:Ns+1)+ks_g(3:Ns+2)).*(Ts_g(2:Ns+1)+Ts_g(3:Ns+2))-Dks_cs(2:Ns+1).*(Ts_g(1:Ns)+Ts_g(2:Ns+1)).*(Cs_g(2:Ns+1)-Cs_g(1:Ns))-1./cs_g(2:Ns+1).*(ks_g(1:Ns)+ks_g(2:Ns+1)).*(Ts_g(1:Ns)+Ts_g(2:Ns+1))),...
        Dks_cs(3:Ns+2).*(phis_g(3:Ns+2)-phis_g(2:Ns+1))+d17*(Dks_cs(3:Ns+2).*(Ts_g(2:Ns+1)+Ts_g(3:Ns+2)).*(Cs_g(3:Ns+2)-Cs_g(2:Ns+1))+1./cs_g(3:Ns+2).*(ks_g(2:Ns+1)+ks_g(3:Ns+2)).*(Ts_g(2:Ns+1)+Ts_g(3:Ns+2)))],0:2,Ns,Ns+2);
    % Dphin_Dcn
    Dv(id_phin+1,id_cs+Ns+1)=Dks_cs(Ns+1)*(phis_g(Ns+2)-phis_g(Ns+1)); 
    Dv(id_phin+1,id_cs+Ns+2)=Dks_cs(Ns+2)*(phis_g(Ns+2)-phis_g(Ns+1));
    Dv(id_phin+1,id_cn+1)=-Dkn_cn(1)*(phin_g(2)-phin_g(1));
    Dv(id_phin+1,id_cn+2)=-Dkn_cn(2)*(phin_g(2)-phin_g(1));
    Dv(id_phin+2:id_phin+Nn+1,id_cn+1:id_cn+Nn+2)=spdiags([d18*(-Dkn_cn(1:Nn).*(phin_g(2:Nn+1)-phin_g(1:Nn)))+d19*(-Dkn_cn(1:Nn).*(Tn_g(1:Nn)+Tn_g(2:Nn+1)).*(Cn_g(2:Nn+1)-Cn_g(1:Nn))+1./cn_g(1:Nn).*(kn_g(1:Nn)+kn_g(2:Nn+1)).*(Tn_g(1:Nn)+Tn_g(2:Nn+1))),...
        d18*(Dkn_cn(2:Nn+1).*(phin_g(3:Nn+2)-phin_g(2:Nn+1))-Dkn_cn(2:Nn+1).*(phin_g(2:Nn+1)-phin_g(1:Nn)))+d19*(Dkn_cn(2:Nn+1).*(Tn_g(2:Nn+1)+Tn_g(3:Nn+2)).*(Cn_g(3:Nn+2)-Cn_g(2:Nn+1))-1./cn_g(2:Nn+1).*(kn_g(2:Nn+1)+kn_g(3:Nn+2)).*(Tn_g(2:Nn+1)+Tn_g(3:Nn+2))-Dkn_cn(2:Nn+1).*(Tn_g(1:Nn)+Tn_g(2:Nn+1)).*(Cn_g(2:Nn+1)-Cn_g(1:Nn))-1./cn_g(2:Nn+1).*(kn_g(1:Nn)+kn_g(2:Nn+1)).*(Tn_g(1:Nn)+Tn_g(2:Nn+1))),...
        d18*(Dkn_cn(3:Nn+2).*(phin_g(3:Nn+2)-phin_g(2:Nn+1)))+d19*(Dkn_cn(3:Nn+2).*(Tn_g(2:Nn+1)+Tn_g(3:Nn+2)).*(Cn_g(3:Nn+2)-Cn_g(2:Nn+1))+1./cn_g(3:Nn+2).*(kn_g(2:Nn+1)+kn_g(3:Nn+2)).*(Tn_g(2:Nn+1)+Tn_g(3:Nn+2)))],0:2,Nn,Nn+2);
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
    Dv(id_phip+2:id_phip+Np+1,id_Tp+1:id_Tp+Np+2)=spdiags([d15*(-Dkp_Tp(1:Np).*(phip_g(2:Np+1)-phip_g(1:Np)))+d16*(-Dkp_Tp(1:Np).*(Tp_g(1:Np)+Tp_g(2:Np+1)).*(Cp_g(2:Np+1)-Cp_g(1:Np))-(kp_g(1:Np)+kp_g(2:Np+1)).*(Cp_g(2:Np+1)-Cp_g(1:Np))),...
        d15*(Dkp_Tp(2:Np+1).*(phip_g(3:Np+2)-phip_g(2:Np+1))-Dkp_Tp(2:Np+1).*(phip_g(2:Np+1)-phip_g(1:Np)))+d16*(Dkp_Tp(2:Np+1).*(Tp_g(2:Np+1)+Tp_g(3:Np+2)).*(Cp_g(3:Np+2)-Cp_g(2:Np+1))+(kp_g(2:Np+1)+kp_g(3:Np+2)).*(Cp_g(3:Np+2)-Cp_g(2:Np+1))-Dkp_Tp(2:Np+1).*(Tp_g(1:Np)+Tp_g(2:Np+1)).*(Cp_g(2:Np+1)-Cp_g(1:Np))-(kp_g(1:Np)+kp_g(2:Np+1)).*(Cp_g(2:Np+1)-Cp_g(1:Np))),...
        d15*(Dkp_Tp(3:Np+2).*(phip_g(3:Np+2)-phip_g(2:Np+1)))+d16*(Dkp_Tp(3:Np+2).*(Tp_g(2:Np+1)+Tp_g(3:Np+2)).*(Cp_g(3:Np+2)-Cp_g(2:Np+1))+(kp_g(2:Np+1)+kp_g(3:Np+2)).*(Cp_g(3:Np+2)-Cp_g(2:Np+1)))],0:2,Np,Np+2);
    Dv(id_phip+Np+2,id_Tp+Np+1)=Dkp_Tp(Np+1)*(phip_g(Np+2)-phip_g(Np+1));
    Dv(id_phip+Np+2,id_Tp+Np+2)=Dkp_Tp(Np+2)*(phip_g(Np+2)-phip_g(Np+1));
    Dv(id_phip+Np+2,id_Ts+1)=-Dks_Ts(1)*(phis_g(2)-phis_g(1));
    Dv(id_phip+Np+2,id_Ts+2)=-Dks_Ts(2)*(phis_g(2)-phis_g(1));
    % Dphis_DTs
    Dv(id_phis+2:id_phis+Ns+1,id_Ts+1:id_Ts+Ns+2)=spdiags([-Dks_Ts(1:Ns).*(phis_g(2:Ns+1)-phis_g(1:Ns))+d17*(Dks_Ts(1:Ns).*(Ts_g(1:Ns)+Ts_g(2:Ns+1)).*(Cs_g(2:Ns+1)-Cs_g(1:Ns))-(ks_g(1:Ns)+ks_g(2:Ns+1)).*(Cs_g(2:Ns+1)-Cs_g(1:Ns))),...
        Dks_Ts(2:Ns+1).*(phis_g(3:Ns+2)-phis_g(2:Ns+1))-Dks_Ts(2:Ns+1).*(phis_g(2:Ns+1)-phis_g(1:Ns))+d17*(Dks_Ts(2:Ns+1).*(Ts_g(2:Ns+1)+Ts_g(3:Ns+2)).*(Cs_g(3:Ns+2)-Cs_g(2:Ns+1))+(ks_g(2:Ns+1)+ks_g(3:Ns+2)).*(Cs_g(3:Ns+2)-Cs_g(2:Ns+1))-Dks_Ts(2:Ns+1).*(Ts_g(1:Ns)+Ts_g(2:Ns+1)).*(Cs_g(2:Ns+1)-Cs_g(1:Ns))-(ks_g(1:Ns)+ks_g(2:Ns+1)).*(Cs_g(2:Ns+1)-Cs_g(1:Ns))),...
        Dks_Ts(3:Ns+2).*(phis_g(3:Ns+2)-phis_g(2:Ns+1))+d17*(Dks_Ts(3:Ns+2).*(Ts_g(2:Ns+1)+Ts_g(3:Ns+2)).*(Cs_g(3:Ns+2)-Cs_g(2:Ns+1))+(ks_g(2:Ns+1)+ks_g(3:Ns+2)).*(Cs_g(3:Ns+2)-Cs_g(2:Ns+1)))],0:2,Ns,Ns+2);
    % Dphin_DTn
    Dv(id_phin+1,id_Ts+Ns+1)=Dks_Ts(Ns+1)*(phis_g(Ns+2)-phis_g(Ns+1));
    Dv(id_phin+1,id_Ts+Ns+2)=Dks_Ts(Ns+2)*(phis_g(Ns+2)-phis_g(Ns+1));
    Dv(id_phin+1,id_Tn+1)=-Dkn_Tn(1)*(phin_g(2)-phin_g(1));
    Dv(id_phin+1,id_Tn+2)=-Dkn_Tn(2)*(phin_g(2)-phin_g(1));
    Dv(id_phin+2:id_phin+Nn+1,id_Tn+1:id_Tn+Nn+2)=spdiags([d18*(-Dkn_Tn(1:Nn).*(phin_g(2:Nn+1)-phin_g(1:Nn)))+d19*(-Dkn_Tn(1:Nn).*(Tn_g(1:Nn)+Tn_g(2:Nn+1)).*(Cn_g(2:Nn+1)-Cn_g(1:Nn))-(kn_g(1:Nn)+kn_g(2:Nn+1)).*(Cn_g(2:Nn+1)-Cn_g(1:Nn))),...
        d18*(Dkn_Tn(2:Nn+1).*(phin_g(3:Nn+2)-phin_g(2:Nn+1))-Dkn_Tn(2:Nn+1).*(phin_g(2:Nn+1)-phin_g(1:Nn)))+d19*(Dkn_Tn(2:Nn+1).*(Tn_g(2:Nn+1)+Tn_g(3:Nn+2)).*(Cn_g(3:Nn+2)-Cn_g(2:Nn+1))+(kn_g(2:Nn+1)+kn_g(3:Nn+2)).*(Cn_g(3:Nn+2)-Cn_g(2:Nn+1))-Dkn_Tn(2:Nn+1).*(Tn_g(1:Nn)+Tn_g(2:Nn+1)).*(Cn_g(2:Nn+1)-Cn_g(1:Nn))-(kn_g(1:Nn)+kn_g(2:Nn+1)).*(Cn_g(2:Nn+1)-Cn_g(1:Nn))),...
        d18*Dkn_Tn(3:Nn+2).*(phin_g(3:Nn+2)-phin_g(2:Nn+1))+d19*(Dkn_Tn(3:Nn+2).*(Tn_g(2:Nn+1)+Tn_g(3:Nn+2)).*(Cn_g(3:Nn+2)-Cn_g(2:Nn+1))+(kn_g(2:Nn+1)+kn_g(3:Nn+2)).*(Cn_g(3:Nn+2)-Cn_g(2:Nn+1)))],0:2,Nn,Nn+2);
    % DTp_Dcp
    Dv(id_Tp+2:id_Tp+Np+1,1:Np+2)=spdiags([-d28./cp_g(1:Np).*kp_g(2:Np+1).*Tp_g(2:Np+1).*(phip_g(3:Np+2)-phip_g(1:Np)),...
        d27*Dkp_cp(2:Np+1).*((phip_g(3:Np+2)-phip_g(1:Np)).^2)+d28*Dkp_cp(2:Np+1).*Tp_g(2:Np+1).*(Cp_g(3:Np+2)-Cp_g(1:Np)).*(phip_g(3:Np+2)-phip_g(1:Np)),...
        d28./cp_g(3:Np+2).*kp_g(2:Np+1).*Tp_g(2:Np+1).*(phip_g(3:Np+2)-phip_g(1:Np))],0:2,Np,Np+2);
    % DTs_Dcs
    Dv(id_Ts+2:id_Ts+Ns+1,id_cs+1:id_cs+Ns+2)=spdiags([-d32./cs_g(1:Ns).*ks_g(2:Ns+1).*Ts_g(2:Ns+1).*(phis_g(3:Ns+2)-phis_g(1:Ns)),...
        d31*Dks_cs(2:Ns+1).*((phis_g(3:Ns+2)-phis_g(1:Ns)).^2)+d32*Dks_cs(2:Ns+1).*Ts_g(2:Ns+1).*(Cs_g(3:Ns+2)-Cs_g(1:Ns)).*(phis_g(3:Ns+2)-phis_g(1:Ns)),...
        d32./cs_g(3:Ns+2).*ks_g(2:Ns+1).*Ts_g(2:Ns+1).*(phis_g(3:Ns+2)-phis_g(1:Ns))],0:2,Ns,Ns+2);
    % DTn_Dcn
    Dv(id_Tn+2:id_Tn+Nn+1,id_cn+1:id_cn+Nn+2)=spdiags([-d36./cn_g(1:Nn).*kn_g(2:Nn+1).*Tn_g(2:Nn+1).*(phin_g(3:Nn+2)-phin_g(1:Nn)),...
        d35*Dkn_cn(2:Nn+1).*((phin_g(3:Nn+2)-phin_g(1:Nn)).^2)+d36*Dkn_cn(2:Nn+1).*Tn_g(2:Nn+1).*(Cn_g(3:Nn+2)-Cn_g(1:Nn)).*(phin_g(3:Nn+2)-phin_g(1:Nn)),...
        d36./cn_g(3:Nn+2).*kn_g(2:Nn+1).*Tn_g(2:Nn+1).*(phin_g(3:Nn+2)-phin_g(1:Nn))],0:2,Nn,Nn+2);
    % DTp_Dsp
    Dv(id_Tp+2:id_Tp+Np+1,id_sp+1:id_sp+Np)=spdiags(DQp_sp(sp,jp,Tp_g(2:Np+1),Tref,d29),0,Np,Np);
    % DTn_Dsn
    Dv(id_Tn+2:id_Tn+Nn+1,id_sn+1:id_sn+Nn)=spdiags(DQn_sn(sn,jn,Tn_g(2:Nn+1),Tref,d37),0,Nn,Nn);
    % DTp_Djp
    Dv(id_Tp+2:id_Tp+Np+1,id_jp+1:id_jp+Np)=spdiags(DQp_jp(sp,Phip,phip_g(2:Np+1),Tp_g(2:Np+1),Tref,d29),0,Np,Np);
    % DTn_Djn
    Dv(id_Tn+2:id_Tn+Nn+1,id_jn+1:id_jn+Nn)=spdiags(DQn_jn(sn,Phin,phin_g(2:Nn+1),Tn_g(2:Nn+1),Tref,d37),0,Nn,Nn);
    % DTp_DPhip
    bottom_TpPhip=[-2*d26*(Phip(3:Np)-Phip(1:Np-2));-2*d26*(Phip(Np)-Phip(Np-1))];
    top_TpPhip=[0;2*d26*(Phip(2)-Phip(1)-h*I/59);2*d26*(Phip(3:Np)-Phip(1:Np-2))];
    Dv(id_Tp+2:id_Tp+Np+1,id_Phip+1:id_Phip+Np)=spdiags(bottom_TpPhip,-1,Np,Np)+spdiags(DQp_Phip(jp,d29),0,Np,Np)+spdiags(top_TpPhip,1,Np,Np);
    Dv(id_Tp+2,id_Phip+1)=Dv(id_Tp+2,id_Phip+1)-2*d26*(Phip(2)-Phip(1)-h*I/59);
    Dv(id_Tp+Np+1,id_Phip+Np)=Dv(id_Tp+Np+1,id_Phip+Np)+2*d26*(Phip(Np)-Phip(Np-2));
    % DTn_DPhin
    bottom_TnPhin=[-2*d34*(Phin(3:Nn)-Phin(1:Nn-2));-2*d34*(Phin(Nn)-h*I/48.24-Phin(Nn-1))];
    top_TnPhin=[0;2*d34*(Phin(2)-Phin(1));2*d34*(Phin(3:Nn)-Phin(1:Nn-2))];
    Dv(id_Tn+2:id_Tn+Nn+1,id_Phin+1:id_Phin+Nn)=spdiags(bottom_TnPhin,-1,Nn,Nn)+spdiags(DQn_Tn(sn,jn,d37),0,Nn,Nn)+spdiags(top_TnPhin,1,Nn,Nn);
    Dv(id_Tn+2,id_Phin+1)=Dv(id_Tn+2,id_Phin+1)-2*d34*(Phin(2)-Phin(1));
    Dv(id_Tn+Nn+1,id_Phin+Nn)=Dv(id_Tn+Nn+1,id_Phin+Nn)+2*d34*(Phin(Nn)-h*I/48.24-Phin(Nn-1));
    % DTp_Dphip
    Dv(id_Tp+2:id_Tp+Np+1,id_phip+1:id_phip+Np+2)=spdiags([-2*d27*kp_g(2:Np+1).*(phip_g(3:Np+2)-phip_g(1:Np))-d28*kp_g(2:Np+1).*Tp_g(2:Np+1).*(Cp_g(3:Np+2)-Cp_g(1:Np)),...
        DQp_sphip(jp,d29),2*d27*kp_g(2:Np+1).*(phip_g(3:Np+2)-phip_g(1:Np))+d28*kp_g(2:Np+1).*Tp_g(2:Np+1).*(Cp_g(3:Np+2)-Cp_g(1:Np))],0:2,Np,Np+2);
    % DTs_Dphis
    Dv(id_Ts+2:id_Ts+Ns+1,id_phis+1:id_phis+Ns+2)=spdiags([-2*d31*ks_g(2:Ns+1).*(phis_g(3:Ns+2)-phis_g(1:Ns))-d32*ks_g(2:Ns+1).*Ts_g(2:Ns+1).*(Cs_g(3:Ns+2)-Cs_g(1:Ns)),...
        zeros(Ns,1),2*d31*ks_g(2:Ns+1).*(phis_g(3:Ns+2)-phis_g(1:Ns))+d32*ks_g(2:Ns+1).*Ts_g(2:Ns+1).*(Cs_g(3:Ns+2)-Cs_g(1:Ns))],0:2,Ns,Ns+2);
    % DTn_Dphin
    Dv(id_Tn+2:id_Tn+Nn+1,id_phin+1:id_phin+Nn+2)=spdiags([-2*d35*kn_g(2:Nn+1).*(phin_g(3:Nn+2)-phin_g(1:Nn))-d36*kn_g(2:Nn+1).*Tn_g(2:Nn+1).*(Cn_g(3:Nn+2)-Cn_g(1:Nn)),...
        DQn_sphin(jn,d37),2*d35*kn_g(2:Nn+1).*(phin_g(3:Nn+2)-phin_g(1:Nn))+d36*kn_g(2:Nn+1).*Tn_g(2:Nn+1).*(Cn_g(3:Nn+2)-Cn_g(1:Nn))],0:2,Nn,Nn+2);
    % DTp_DTp
    Dv(id_Tp+2:id_Tp+Np+1,id_Tp+2:id_Tp+Np+1)=spdiags(d27*Dkp_Tp(2:Np+1).*((phip_g(3:Np+2)-phip_g(1:Np)).^2)+d28*(Cp_g(3:Np+2)-Cp_g(1:Np)).*(phip_g(3:Np+2)-phip_g(1:Np)).*(Dkp_Tp(2:Np+1).*Tp_g(2:Np+1)+kp_g(2:Np+1))+DQp_Tp(sp,jp,d29),0,Np,Np);
    % DTs_DTs
    Dv(id_Ts+2:id_Ts+Ns+1,id_Ts+2:id_Ts+Ns+1)=spdiags(d31*Dks_Ts(2:Ns+1).*((phis_g(3:Ns+2)-phis_g(1:Ns)).^2)+d32*(Cs_g(3:Ns+2)-Cs_g(1:Ns)).*(phis_g(3:Ns+2)-phis_g(1:Ns)).*(Dks_Ts(2:Ns+1).*Ts_g(2:Ns+1)+ks_g(2:Ns+1)),0,Ns,Ns);
    % DTn_DTn
    Dv(id_Tn+2:id_Tn+Nn+1,id_Tn+2:id_Tn+Nn+1)=spdiags(d35*Dkn_Tn(2:Nn+1).*((phin_g(3:Nn+2)-phin_g(1:Nn)).^2)+d36*(Cn_g(3:Nn+2)-Cn_g(1:Nn)).*(phin_g(3:Nn+2)-phin_g(1:Nn)).*(Dkn_Tn(2:Nn+1).*Tn_g(2:Nn+1)+kn_g(2:Nn+1))+DQn_Tn(sn,jn,d37),0,Nn,Nn);
end