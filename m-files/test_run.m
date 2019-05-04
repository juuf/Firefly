%% one node synchronises to external signal
[phaseOut,omegaOut,phi_ext] = fireflySimulation3_Sich_fkt_('nodes',1,'fire',0,'simlength',8e4);
% plot(phi_ext); hold on; plot(test)
figure(1)
clf
yyaxis left
plot(phaseOut,'b-'); hold on; plot(phi_ext,'r-'); xlabel('Time'); ylabel('Phase'); set(gca,'YColor','k')
yyaxis right
plot(omegaOut,'g-'); ylabel('Frequency'); set(gca,'YColor','g')
%%
[phaseOut,omegaOut] = fireflySimulation2('nodes',2);
figure(2)
clf
yyaxis left
plot(phaseOut(:,1),'b-'); hold on; plot(phaseOut(:,2),'r-'); xlabel('Time'); ylabel('Phase'); set(gca,'YColor','k')
yyaxis right
plot(omegaOut,'g-'); ylabel('Frequency'); set(gca,'YColor','g')