function output=snNL(Tn,Tref)
    output=-(2e-6)./(5.*(3.9e-14).*exp((-5000./8.314472).*(1./Tn-1./Tref))); 
end