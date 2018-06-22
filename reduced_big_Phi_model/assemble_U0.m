function U0=assemble_U0(U,Na,Np,Ns,Nn,Nz,Nr)
    % indices for convenience
    [id_cs,id_cn,id_ap,id_an,id_sp,id_sn,id_jp,id_jn,id_Phip,id_Phin,id_phip,id_phis,id_phin,id_Ta,id_Tp,id_Ts,id_Tn,id_Tz]=indices(Na,Np,Ns,Nn,Nr);
    
    % variables
    cp_g=U(1:Np+2); 
    cs_g=U(id_cs+1:id_cs+Ns+2); 
    cn_g=U(id_cn+1:id_cn+Nn+2); 
    ap=U(id_ap+1:id_ap+Np*Nr); 
    an=U(id_an+1:id_an+Nn*Nr); 
    Ta_g=U(id_Ta+1:id_Ta+Na+2);
    Tp_g=U(id_Tp+1:id_Tp+Np+2);
    Ts_g=U(id_Ts+1:id_Ts+Ns+2);
    Tn_g=U(id_Tn+1:id_Tn+Nn+2);
    Tz_g=U(id_Tz+1:id_Tz+Nz+2);
    
    % constructing U0
    U0=zeros(length(U),1);
    U0(2:Np+1)=cp_g(2:Np+1);
    U0(id_cs+2:id_cs+Ns+1)=cs_g(2:Ns+1);
    U0(id_cn+2:id_cn+Nn+1)=cn_g(2:Nn+1);
    U0(id_ap+1:id_ap+Np*Nr)=ap;
    U0(id_an+1:id_an+Nn*Nr)=an;
    U0(id_Ta+2:id_Ta+Na+1)=Ta_g(2:Na+1);
    U0(id_Tp+2:id_Tp+Np+1)=Tp_g(2:Np+1);
    U0(id_Ts+2:id_Ts+Ns+1)=Ts_g(2:Ns+1);
    U0(id_Tn+2:id_Tn+Nn+1)=Tn_g(2:Nn+1);
    U0(id_Tz+2:id_Tz+Nz+1)=Tz_g(2:Nz+1);
end