function output=jnNL(cn,sn,Phin,phin,Tn,Tref)
    tn=sn./30555;
    G=-(0.7222+0.1387.*tn+0.029.*sqrt(tn)-0.0172./tn+0.0019./tn.^(1.5)+0.2808.*exp(0.9-15.*tn)-0.7984.*exp(0.4465.*tn-0.4108))-(Tn-Tref).*...
        0.001.*(0.005268056+3.299265709.*tn-91.79325798.*tn.^2+1004.911008.*tn.^3-5812.278127.*tn.^4+19329.7549.*tn.^5-...
        37147.8947.*tn.^6+38379.18127.*tn.^7-16515.05308.*tn.^8)./(1-48.09287227.*tn+1017.234804.*tn.^2-10481.80419.*tn.^3+59431.3.*tn.^4-...
        195881.6488.*tn.^5+374577.3152.*tn.^6-385821.1607.*tn.^7+165705.8597.*tn.^8);
    output=-2.*5.031.*(10.^(-11)).*exp((-5000./8.314472).*(1./Tn-1./Tref)).*sqrt(cn.*(30555-sn).*sn).*sinh(0.5.*96485./(8.314472.*Tn).*(Phin-phin+G));
end