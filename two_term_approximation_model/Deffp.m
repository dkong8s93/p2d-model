function Dp=Deffp(cp,Tp)
    Dp=(0.385.^4).*(10.^(-4)).*10.^(-4.43-54./(Tp-229-(5e-3).*cp)-(0.22e-3).*cp);
end