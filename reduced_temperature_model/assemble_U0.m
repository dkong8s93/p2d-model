function U0=assemble_U0(U,Np,Ns,Nn,Nr)
    % indices for convenience
    [id_cs,id_cn,id_ap,id_an,id_sp,id_sn,id_jp,id_jn,id_Phip,id_Phin,id_phip,id_phis,id_phin,id_T]=indices(Np,Ns,Nn,Nr);
    
    % variables
    cp_g=U(1:Np+2); 
    cs_g=U(id_cs+1:id_cs+Ns+2); 
    cn_g=U(id_cn+1:id_cn+Nn+2); 
    ap=U(id_ap+1:id_ap+Np*Nr); 
    an=U(id_an+1:id_an+Nn*Nr); 
    T=U(id_T+1);
    
    % constructing U0
    U0=zeros(length(U),1);
    U0(2:Np+1)=cp_g(2:Np+1);
    U0(id_cs+2:id_cs+Ns+1)=cs_g(2:Ns+1);
    U0(id_cn+2:id_cn+Nn+1)=cn_g(2:Nn+1);
    U0(id_ap+1:id_ap+Np*Nr)=ap;
    U0(id_an+1:id_an+Nn*Nr)=an;
    U0(id_T+1)=T;
end