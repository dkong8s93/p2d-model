function Ds=Deffs(cs,Ts)
    Ds=(0.724.^4).*(10.^(-4)).*10.^(-4.43-54./(Ts-229-(5e-3).*cs)-(0.22e-3).*cs);
end