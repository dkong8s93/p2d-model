function output=apNL(jp,Tp,Tref,d44)
    output=d44*jp./exp((-5000./8.314472).*(1./Tp-1./Tref));
end