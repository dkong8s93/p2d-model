function kn=kappan(cn,Tn)
    kn=(0.485.^4).*(10.^(-4)).*cn.*(-10.5+(0.668e-3).*cn+(0.494e-6).*(cn.^2)+(0.074-(1.78e-5).*cn-(8.86e-10).*(cn.^2)).*Tn+((-6.96e-5)+(2.8e-8).*cn).*(Tn.^2)).^2;
end