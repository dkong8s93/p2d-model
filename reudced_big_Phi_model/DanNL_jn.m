function output=DanNL_jn(Tn,Tref,d46)
    output=d46./exp((-5000./8.314472).*(1./Tn-1./Tref));
end