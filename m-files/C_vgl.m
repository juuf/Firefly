%%
c_phi_ext = table2array(importfile('./phi_ext.txt'));

figure
plot(c_phi_ext); hold on
plot(phi_ext,"--")
%%
figure
c_phi = table2array(importfile('./phi.txt'));
plot(c_phi);hold on
plot(phase,"--");

%%
figure
c_omegas = table2array(importfile('./omegas.txt'));
plot(c_omegas);hold on
plot(omegas,"--");
%%
figure(3)
clf
yyaxis left
plot(c_phi,'b-'); hold on; plot(phi_ext,'r-'); xlabel('Time'); ylabel('Phase'); set(gca,'YColor','k')
yyaxis right
plot(c_omegas,'g-'); ylabel('Frequency'); set(gca,'YColor','g')
