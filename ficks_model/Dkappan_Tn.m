function Dkn_DTn=Dkappan_Tn(cn,Tn)
    Dkn_DTn=-(3266153979273035.*cn.*((164176022256015.*cn)./9223372036854775808 + (2142216552357123.*cn.^2)./2417851639229258349412352 - 2.*Tn.*((2115620184325601.*cn)./75557863725914323419136 - 5135573550120739./73786976294838206464) - 37./500).*((3080606260309495.*cn)./4611686018427387904 + Tn.^2.*((2115620184325601.*cn)./75557863725914323419136 - 5135573550120739./73786976294838206464) + (4665698085075209.*cn.^2)./9444732965739290427392 - Tn.*((2142216552357123.*cn.^2)./2417851639229258349412352 + (164176022256015.*cn)./9223372036854775808 - 37./500) - 21./2))./295147905179352825856;
end