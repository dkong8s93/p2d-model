function output=DapNL_jp(Tp,Tref,d44)
    output=d44./exp((-5000./8.314472).*(1./Tp-1./Tref));
end