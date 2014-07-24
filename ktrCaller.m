% Mathematical model of KTR system, as published by Regot et al. 2014 (DOI: 10.1016/j.cell.2014.04.039).
% Parameters and initial conditions for JNK KTR in 3T3 cells, along with an example simulation.
% 2014-07-21, JJH.

%% time-independent parameters
reporterTotal = 0.4; % (uM) total concentration of reporter if it were all in cytosol
kv = 4; % ratio of cytosolic volume to nuclear volume

kiu = 0.44; % (1/min) import of unphosphorylated reporter
keu = 0.18; % (1/min) export of unphosphorylated reporter
kip = 0.16; % (1/min) import of phosphorylated reporter
kep = 0.2; % (1/min) export of phosphorylated reporter

kcat = 20; % (1/min) kcat of kinase and reporter
Km = 3; % (uM) Km of kinase and reporter

kdc = 0.03; % (uM/min) dephosphorylation Vmax in cytosol
kdn = kdc; % (uM/min) dephosphorylation Vmax in nucleus
Kmd = 0.1; % (uM) dephosphorylation Km

params = v2struct({'reporterTotal', 'kv', 'kiu', 'keu', 'kip', 'kep', 'kcat', 'Km', 'kdc', 'kdn', 'Kmd', 'fieldnames'});

%% time-dependent parameters
% tKin (min) time vector for kinC and kinN, should have the same limits as tSpan
% kinC (uM) time-dependent concentration of active kinase in cytosol
% kinN (uM) time-dependent concentration of active kinase in nucleus

%% initial conditions
% y(1) unphosphorylated reporter in cytosol
% y(2) unphosphorylated reporter in nucleus
% y(3) phosphorylated reporter in cytosol
% y(4) phosphorylated reporter in nucleus
yStart = [reporterTotal; 0; 0; 0];

%% get model to steady state (set the basal active kinase concentrations)
tStart = [-100, 0];
tKin0 = tStart;
kinC0 = 0.01 * [1, 1];
kinN0 = kinC0;

options = odeset('RelTol', 1e-4);
[t0, y0] = ode15s(@ktrEqns, tStart, yStart, options, params, tKin0, kinC0, kinN0);

%% add input
tSpan = [0, 30];
tKin1 = tSpan;
kinC1 = 0.1 * [1, 1];
kinN1 = kinC1;

[tStim, yStim] = ode15s(@ktrEqns, tSpan, y0(end,:), options, params, tKin1, kinC1, kinN1);

%% add inhibitor
tSpan = [tSpan(2), tSpan(2)+50];
tKin2 = tSpan;
kinC2 = [0, 0];
kinN2 = kinC2;

[tInhib, yInhib] = ode15s(@ktrEqns, tSpan, yStim(end,:), options, params, tKin2, kinC2, kinN2);

%% construct outputs
t = [t0(1:end-1); tStim; tInhib(2:end)]';
y = [y0(1:end-1,:); yStim; yInhib(2:end,:)]';

tKin = [tKin0, tKin1, tKin2];
kinC = [kinC0, kinC1, kinC2];
kinN = [kinN0, kinN1, kinN2];

%% plotting
figure
plot(tKin, [kinC; kinN])
xlim([-10, 80])
ylabel('active kinase (\muM)')
legend({'kinC', 'kinN'})

figure
plot(t, y(1:4,:))
xlim([-10, 80])
ylabel('conc (\muM)')
legend({'uc', 'un', 'pc', 'pn'})

figure
plot(t, (y(1,:) + y(3,:)) ./ (y(2,:) + y(4,:)))
xlim([-10, 80])
ylim([0.2, 1.2])
ylabel('cytosolic / nuclear')
