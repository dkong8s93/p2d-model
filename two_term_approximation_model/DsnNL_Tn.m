function output=DsnNL_Tn(Tn,Tref)
    output=(6244890665979544178027278379745.*exp(1322407036521381./(2199023255552.*Tn) - 1322407036521381./(2199023255552.*Tref)))./(1012497887414291267584.*Tn.^2);
end