function output=DjnNL_Phin(cn,sn,Phin,phin,Tn,Tref)
    output=-(1467079944916146455.*cosh((27158113127927644160.*(phin - Phin + (1387.*sn)./305550000 + (351.*exp(9./10 - sn./2037))./1250 - (499.*exp((893.*sn)./61110000 - 1027./2500))./625 + (29.*(sn./30555).^(1./2))./1000 + 19./(10000.*(sn./30555).^(3./2)) - 262773./(500.*sn) + ((Tn./1000 - Tref./1000).*(- (4539623223699461.*sn.^8)./208832525047565003860943627541738317414400000000 + (659349313607623.*sn.^7)./427165204237369096426410627437690880000000 - (1021113554251211.*sn.^6)./22368330118795305327516184058265600000 + (110694220180261.*sn.^5)./152514114702962590189250150400000 - (639066738450461.*sn.^4)./95836066185464956034482176000 + (1104911338176101.*sn.^3)./31365101026170825080832000 - (1291875257605787.*sn.^2)./13139364854687827230720 + (103184526511449.*sn)./955607545932677120 + 759206881234141./144115188075855872))./((1138721997067303.*sn.^8)./5220813126189125096523590688543457935360000000 - (1104726178207507.*sn.^7)./71194200706228182737735104572948480000000 + (47668068699481.*sn.^6)./103557083883311598738500852121600000 - (6730442203860461.*sn.^5)./915084688217775541135500902400000 + (453787537526693.*sn.^4)./6655282373990621946839040000 - (5762432793488099.*sn.^3)./15682550513085412540416000 + (8947691961411523.*sn.^2)./8212103034179892019200 - (161154048835081.*sn)./102386522778501120 + 1) + 3611./5000))./(4680631625122803.*Tn)).*exp(1322407036521381./(2199023255552.*Tref) - 1322407036521381./(2199023255552.*Tn)).*(-cn.*sn.*(sn - 30555)).^(1./2))./(2512894969315721358606336.*Tn);
end