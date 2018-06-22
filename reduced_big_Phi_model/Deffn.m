function Dn=Deffn(cn,Tn)
    Dn=(0.485.^4).*(10.^(-4)).*10.^(-4.43-54./(Tn-229-(5e-3).*cn)-(0.22e-3).*cn);
end