function DDn_DTn=DDeffn_Tn(cn,Tn)
    DDn_DTn=(88186157440371945.*10.^(54./(cn./200 - Tn + 229) - (11.*cn)./50000 - 443./100).*log(10))./(295147905179352825856.*(cn./200 - Tn + 229).^2);
end