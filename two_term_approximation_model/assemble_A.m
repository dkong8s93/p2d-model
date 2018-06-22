function A=assemble_A(d2,d5,d6,d7,d13,d14,d20,d22,d23,d25,d30,d33,d38,d40,d41,Na,Np,Ns,Nn,Nz)
    dim=Na+7*Np+3*Ns+7*Nn+Nz+2*11; % dimension of matrix A
    
    A=zeros(dim);
    
    % indices for convenience
    [id_cs,id_cn,id_ap,id_an,id_sp,id_sn,id_jp,id_jn,id_Phip,id_Phin,id_phip,id_phis,id_phin,id_Ta,id_Tp,id_Ts,id_Tn,id_Tz]=indices(Na,Np,Ns,Nn);
    
    % A0101
    A0101=eye(Np+2); % check
    A0101(1,1)=-1; A0101(1,2)=1; A0101(end,end)=0; % check
    A(1:Np+2,1:Np+2)=A0101; % check
    % A0108
    A(2:Np+1,id_jp+1:id_jp+Np)=d2*eye(Np); % check
    % A0201
    A(id_cs+1,id_cs-1)=1; A(id_cs+1,id_cs)=1; % check
    % A0202
    A0202=eye(Ns+2); % check
    A0202(1,1)=-1; A0202(1,2)=-1; A0202(end,end-1)=1; % check
    A(id_cs+1:id_cs+Ns+2,id_cs+1:id_cs+Ns+2)=A0202; % check
    % A0203
    A(id_cn,id_cn+1)=-1; A(id_cn,id_cn+2)=-1; % check
    % A0303
    A0303=eye(Nn+2); % check
    A0303(1,1)=0; A0303(end,end-1)=-1; % check
    A(id_cn+1:id_cn+Nn+2,id_cn+1:id_cn+Nn+2)=A0303; % check
    % A0309
    A(id_cn+2:id_cn+Nn+1,id_jn+1:id_jn+Nn)=d5*eye(Nn); % check
    % A0404
    A(id_ap+1:id_ap+Np,id_ap+1:id_ap+Np)=eye(Np); % check
    % A0408
    A(id_ap+1:id_ap+Np,id_jp+1:id_jp+Np)=d6*eye(Np); % check
    % A0505
    A(id_an+1:id_an+Nn,id_an+1:id_an+Nn)=eye(Nn); % check
    % A0509
    A(id_an+1:id_an+Nn,id_jn+1:id_jn+Nn)=d7*eye(Nn); % check
    % A0604
    A(id_sp+1:id_sp+Np,id_ap+1:id_ap+Np)=eye(Np); % check
    % A0606;
    A(id_sp+1:id_sp+Np,id_sp+1:id_sp+Np)=-eye(Np); % check
    % A0705
    A(id_sn+1:id_sn+Nn,id_an+1:id_an+Nn)=eye(Nn); % check
    % A0707
    A(id_sn+1:id_sn+Nn,id_sn+1:id_sn+Nn)=-eye(Nn); % check
    % A0808
    A(id_jp+1:id_jp+Np,id_jp+1:id_jp+Np)=eye(Np); % check
    % A0909
    A(id_jn+1:id_jn+Nn,id_jn+1:id_jn+Nn)=eye(Nn); % check
    % A1008
    A(id_Phip+1:id_Phip+Np,id_jp+1:id_jp+Np)=d13*eye(Np); % check
    % A1010
    A1010=spdiags([ones(Np,1),-2*ones(Np,1),ones(Np,1)],-1:1,Np,Np); % check
    A1010(1,1)=-1; A1010(end,end)=-1; % check
    A(id_Phip+1:id_Phip+Np,id_Phip+1:id_Phip+Np)=A1010; % check
    % A1109
    A(id_Phin+1:id_Phin+Nn,id_jn+1:id_jn+Nn)=d14*eye(Nn); % check
    % A1111
    A1111=spdiags([ones(Nn,1),-2*ones(Nn,1),ones(Nn,1)],-1:1,Nn,Nn); % check
    A1111(1,1)=-1; A1111(end,end)=-1; % check
    A(id_Phin+1:id_Phin+Nn,id_Phin+1:id_Phin+Nn)=A1111; % check
    % A1208
    A(id_phip+2:id_phip+Np+1,id_jp+1:id_jp+Np)=eye(Np); % check
    % A1212
    A(id_phip+1,id_phip+1)=-1; A(id_phip+1,id_phip+2)=1; % check
    % A1312
    A(id_phis+1,id_phis-1)=1; A(id_phis+1,id_phis)=1; % check
    % A1313
    A(id_phis+1,id_phis+1)=-1; A(id_phis+1,id_phis+2)=-1; % check
    A(id_phin,id_phin-1)=1; A(id_phin,id_phin)=1; % check
    % A1314
    A(id_phin,id_phin+1)=-1; A(id_phin,id_phin+2)=-1; % check
    % A1409
    A(id_phin+2:id_phin+Nn+1,id_jn+1:id_jn+Nn)=eye(Nn);
    % A1414
    A(id_phin+Nn+2,id_phin+Nn+1)=1; A(id_phin+Nn+2,id_phin+Nn+2)=1;
    % A1515
    A1515=spdiags([d20*ones(Na+2,1),(1-2*d20)*ones(Na+2,1),d20*ones(Na+2,1)],-1:1,(Na+2),(Na+2));
    A1515(1,1)=d22; A1515(1,2)=d23; A1515(Na+2,Na+1)=1; A1515(Na+2,Na+2)=1;
    A(id_Ta+1:id_Ta+Na+2,id_Ta+1:id_Ta+Na+2)=A1515;
    % A1516
    A(id_Ta+Na+2,id_Tp+1)=-1; A(id_Ta+Na+2,id_Tp+2)=-1;
    % A1615
    A(id_Tp+1,id_Ta+Na+1)=-237; A(id_Tp+1,id_Ta+Na+2)=237;
    % A1616
    A1616=spdiags([d25*ones(Np+2,1),(1-2*d25)*ones(Np+2,1),d25*ones(Np+2,1)],-1:1,(Np+2),(Np+2));
    A1616(1,1)=2.1; A1616(1,2)=-2.1; A1616(Np+2,Np+1)=-2.1; A1616(Np+2,Np+2)=2.1;
    A(id_Tp+1:id_Tp+Np+2,id_Tp+1:id_Tp+Np+2)=A1616; 
    % A1617
    A(id_Tp+Np+2,id_Ts+1)=0.16; A(id_Tp+Np+2,id_Ts+2)=-0.16;
    % A1716
    A(id_Ts+1,id_Tp+Np+1)=1; A(id_Ts+1,id_Tp+Np+2)=1; 
    % A1717
    A1717=spdiags([d30*ones(Ns+2,1),(1-2*d30)*ones(Ns+2,1),d30*ones(Ns+2,1)],-1:1,(Ns+2),(Ns+2));
    A1717(1,1)=-1; A1717(1,2)=-1; A1717(Ns+2,Ns+1)=1; A1717(Ns+2,Ns+2)=1;
    A(id_Ts+1:id_Ts+Ns+2,id_Ts+1:id_Ts+Ns+2)=A1717;
    % A1718
    A(id_Ts+Ns+2,id_Tn+1)=-1; A(id_Ts+Ns+2,id_Tn+2)=-1;
    % A1817
    A(id_Tn+1,id_Ts+Ns+1)=-0.16; A(id_Tn+1,id_Ts+Ns+2)=0.16;
    % A1818
    A1818=spdiags([d33*ones(Nn+2,1),(1-2*d33)*ones(Nn+2,1),d33*ones(Nn+2,1)],-1:1,(Nn+2),(Nn+2));
    A1818(1,1)=1.7; A1818(1,2)=-1.7; A1818(Nn+2,Nn+1)=-1.7; A1818(Nn+2,Nn+2)=1.7;
    A(id_Tn+1:id_Tn+Nn+2,id_Tn+1:id_Tn+Nn+2)=A1818;
    % A1819
    A(id_Tn+Nn+2,id_Tz+1)=401; A(id_Tn+Nn+2,id_Tz+2)=-401;
    % A1918
    A(id_Tz+1,id_Tn+Nn+1)=1; A(id_Tz+1,id_Tn+Nn+2)=1; 
    % A1919
    A1919=spdiags([d38*ones(Nz+2,1),(1-2*d38)*ones(Nz+2,1),d38*ones(Nz+2,1)],-1:1,(Nz+2),(Nz+2));
    A1919(1,1)=-1; A1919(1,2)=-1; A1919(Nz+2,Nz+1)=d40; A1919(Nz+2,Nz+2)=d41;
    A(id_Tz+1:id_Tz+Nz+2,id_Tz+1:id_Tz+Nz+2)=A1919;
end