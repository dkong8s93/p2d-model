function output=DQp_sp(sp,jp,Tp,Tref,d29)
    output=-d29.*jp.*(((38125257829877025.*sp.^9)./1166577203843018372260999966340558972783512283121166013431808 - (1016984484018389.*sp.^7)./13716356771358538095117022418028474899071850512384 + (9048778386456969.*sp.^5)./82572231429879073461363762881544306294784 - (3528280036975051.*sp.^3)./15533856577596066411225559334912 + (88699.*sp)./1328907458000)./((2399.*sp.^10)./3315611831575945337422694951614400260903927014400 - (73083.*sp.^8)./49899815244711273754062940199322747136000 + (37311.*sp.^6)./18774751749760769930927327296000 - (19883.*sp.^4)./1765995031928021764000 + (18993.*sp.^2)./2657814916000 - 1) + ((Tp./1000 - Tref./1000).*((16524920140432371.*sp.^2)./1234175360628340924474917388288 - (6145389974508567.*sp)./5984867132658673063559168 + 4181023998819723./232178575189458550784))./((6865186082835357.*sp.^4)./15906669135458372005094972758949888 - (5530596059357619.*sp.^3)./77135960039271307779682336768 + (6460617402680619.*sp.^2)./1496216783164668265889792 - (3187129838682343.*sp)./29022321898682318848 + 1) - (((2399.*sp.^9)./331561183157594533742269495161440026090392701440 - (73083.*sp.^7)./6237476905588909219257867524915343392000 + (111933.*sp.^5)./9387375874880384965463663648000 - (19883.*sp.^3)./441498757982005441000 + (18993.*sp)./1328907458000).*((7625051565975405.*sp.^10)./2333154407686036744521999932681117945567024566242332026863616 - (1016984484018389.*sp.^8)./109730854170868304760936179344227799192574804099072 + (3016259462152323.*sp.^6)./165144462859758146922727525763088612589568 - (3528280036975051.*sp.^4)./62135426310384265644902237339648 + (88699.*sp.^2)./2657814916000 - 582./125))./((2399.*sp.^10)./3315611831575945337422694951614400260903927014400 - (73083.*sp.^8)./49899815244711273754062940199322747136000 + (37311.*sp.^6)./18774751749760769930927327296000 - (19883.*sp.^4)./1765995031928021764000 + (18993.*sp.^2)./2657814916000 - 1).^2 - (Tp.*((16524920140432371.*sp.^2)./1234175360628340924474917388288000 - (6145389974508567.*sp)./5984867132658673063559168000 + 4181023998819723./232178575189458550784000))./((6865186082835357.*sp.^4)./15906669135458372005094972758949888 - (5530596059357619.*sp.^3)./77135960039271307779682336768 + (6460617402680619.*sp.^2)./1496216783164668265889792 - (3187129838682343.*sp)./29022321898682318848 + 1) - ((Tp./1000 - Tref./1000).*((6865186082835357.*sp.^3)./3976667283864593001273743189737472 - (16591788178072857.*sp.^2)./77135960039271307779682336768 + (6460617402680619.*sp)./748108391582334132944896 - 3187129838682343./29022321898682318848).*((5508306713477457.*sp.^3)./1234175360628340924474917388288 - (6145389974508567.*sp.^2)./11969734265317346127118336 + (4181023998819723.*sp)./232178575189458550784 - 3594251507571897./18014398509481984))./((6865186082835357.*sp.^4)./15906669135458372005094972758949888 - (5530596059357619.*sp.^3)./77135960039271307779682336768 + (6460617402680619.*sp.^2)./1496216783164668265889792 - (3187129838682343.*sp)./29022321898682318848 + 1).^2 + (Tp.*((6865186082835357.*sp.^3)./3976667283864593001273743189737472 - (16591788178072857.*sp.^2)./77135960039271307779682336768 + (6460617402680619.*sp)./748108391582334132944896 - 3187129838682343./29022321898682318848).*((5508306713477457.*sp.^3)./1234175360628340924474917388288000 - (6145389974508567.*sp.^2)./11969734265317346127118336000 + (4181023998819723.*sp)./232178575189458550784000 - 3594251507571897./18014398509481984000))./((6865186082835357.*sp.^4)./15906669135458372005094972758949888 - (5530596059357619.*sp.^3)./77135960039271307779682336768 + (6460617402680619.*sp.^2)./1496216783164668265889792 - (3187129838682343.*sp)./29022321898682318848 + 1).^2);
end