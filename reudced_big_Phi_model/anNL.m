function output=anNL(jn,Tn,Tref,d46)
    output=d46*jn./exp((-5000./8.314472).*(1./Tn-1./Tref));
end