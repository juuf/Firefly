%% loop outside
clear 
clear global

global pas_

pas_.in = 0;
pas_.interactions = zeros(1,11);
pas_.phi = zeros(1,1);
pas_.C = 1;
pas_.E = 1;
pas_.domega = zeros(1,2); 
pas_.fire = 0; 
rng(1234)

simlength = 8e4;
sr =1e3;

phi_ext = zeros(simlength,1);
omegas = zeros(simlength,1);
phase = zeros(simlength,1);
temp2 = rand(2,1);
pas_.phi = temp2(1); phase(1) = pas_.phi;
phi_ext(1) = temp2(2);
temp = 2*2.^((rand(2,1)-0.5)*2);
pas_.omega = temp(1); omegas(1)=pas_.omega;
omega_ext= temp(2);
f = zeros(1);



for i_=2:simlength
    % external oscillator with phase phi_ext and constant frequency omega_ext
    if phi_ext(i_-1) < 1
        phi_ext(i_) = phi_ext(i_-1)+omega_ext/sr; % tooth saw -> increase (phase=amplitude)
        if phi_ext(i_)>=1
            phi_ext(i_)=1;
        end
    else
        f(1)=mod(f(1)+1,2);
        phi_ext(i_)=0;
    end
    pas_.phi_ext = phi_ext(i_);
    pas_.ii = i_;
    [phaseOut,omegaOut] = fireflySimulation4('nodes',1,'f',f,'samplerate',sr);
    omegas(i_)=omegaOut;
    phase(i_)=phaseOut;
end
%%

% plot(phi_ext); hold on; plot(test)
figure(2)
clf
yyaxis left
plot(phase,'b-'); hold on; plot(phi_ext,'r-'); xlabel('Time'); ylabel('Phase'); set(gca,'YColor','k')
yyaxis right
plot(omegas,'g-'); ylabel('Frequency'); set(gca,'YColor','g')
% %%



%%
% [phaseOut,omegaOut,phi_ext] = fireflySimulation3_Sich_fkt('nodes',1,'fire',0,'simlength',12e4);
% % plot(phi_ext); hold on; plot(test)
% figure(1)
% clf
% yyaxis left
% plot(phaseOut,'b-'); hold on; plot(phi_ext,'r-'); xlabel('Time'); ylabel('Phase'); set(gca,'YColor','k')
% yyaxis right
% plot(omegaOut,'g-'); ylabel('Frequency'); set(gca,'YColor','g')
% %%
% 
% % [phaseOut,omegaOut] = fireflySimulation3_Sich_fkt('nodes',1,'fire',0);
% [phaseOut,omegaOut] = fireflySimulation2('nodes',6);




