function Dks_DTs=Dkappas_Ts(cs,Ts)
    Dks_DTs=-(8109497979584201.*cs.*((164176022256015.*cs)./9223372036854775808 + (2142216552357123.*cs.^2)./2417851639229258349412352 - 2.*Ts.*((2115620184325601.*cs)./75557863725914323419136 - 5135573550120739./73786976294838206464) - 37./500).*((3080606260309495.*cs)./4611686018427387904 + Ts.^2.*((2115620184325601.*cs)./75557863725914323419136 - 5135573550120739./73786976294838206464) + (4665698085075209.*cs.^2)./9444732965739290427392 - Ts.*((2142216552357123.*cs.^2)./2417851639229258349412352 + (164176022256015.*cs)./9223372036854775808 - 37./500) - 21./2))./147573952589676412928;
end