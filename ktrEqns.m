% Mathematical model of KTR system, as published by Regot et al. 2014 (DOI: 10.1016/j.cell.2014.04.039).
% 2014-07-21, JJH.

function dy = ktrEqns(t, y, p, tKin, kinC, kinN)

%% parameters
% kv -> ratio of cytosolic volume to nuclear volume

% kiu -> import of unphosphorylated
% keu -> export of unphosphorylated
% kip -> import of phosphorylated
% kep -> export of phosphorylated

% kcat -> kcat of kinase and reporter
% Km -> Km of kinase and reporter
% kdc -> dephosphorylation in cytosol
% kdn -> dephosphorylation in nucleus
% Kmd -> dephosphorylation Km

% tKin -> time vector for kinC and kinN
% kinC -> time-dependent concentration of active kinase in cytosol
% kinN -> time-dependent concentration of active kinase in nucleus

%% state variables
% y(1) unphosphorylated reporter in cytosol
% y(2) unphosphorylated reporter in nucleus
% y(3) phosphorylated reporter in cytosol
% y(4) phosphorylated reporter in nucleus

%% equations
v2struct(p);

kincNow = interp1(tKin, kinC, t);
kinnNow = interp1(tKin, kinN, t);

dy = zeros(4, 1);

dy(1) = - kincNow*kcat*y(1)/(y(1)+Km) + kdc*y(3)/(y(3)+Kmd) - kiu*y(1) + keu*y(2);

dy(2) = - kinnNow*kcat*y(2)/(y(2)+Km) + kdn*y(4)/(y(4)+Kmd) + kv*kiu*y(1) - kv*keu*y(2);

dy(3) = kincNow*kcat*y(1)/(y(1)+Km) - kdc*y(3)/(y(3)+Kmd) - kip*y(3) + kep*y(4);

dy(4) = kinnNow*kcat*y(2)/(y(2)+Km) - kdn*y(4)/(y(4)+Kmd) + kv*kip*y(3) - kv*kep*y(4);
