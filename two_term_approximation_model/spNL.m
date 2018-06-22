function output=spNL(Tp,Tref)
    output=-(2e-6)./(5.*(1e-14).*exp((-5000./8.314472).*(1./Tp-1./Tref)));
end