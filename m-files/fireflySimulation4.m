function [phaseOut,omegaOut] = fireflySimulation4(varargin)
% loop poutside 

global p pas_

p.simulationlength = 6e4;
p.filterlength = 8;
p.nodes = 5;       % number of nodes
p.samplerate = 1e3;
% p.plotlength = 10000;
p.tref = 5; %ms
p.doublingthreshold = 20;
% p.animspeed = 0.2;
% p.animate = 0;
p.alpha = 0.4;
p.beta = 0.4;
p.filtertype = 'median';
% p.pas_.fire=0;

for h_ = 2 : 2 : length(varargin)
    switch lower(varargin{h_-1})
        case {'alpha','alp'}
            p.alpha = varargin{h_};
        case {'beta','bet'}
            p.beta = varargin{h_};
        case {'simulationlength','simlength'}
            p.simulationlength = varargin{h_};
        case {'filterlength'}
            p.filterlength = varargin{h_};
        case {'filterype'}
            p.filtertype = varargin{h_};
        case {'nodes'}
            p.nodes = varargin{h_};
        case {'samplerate'}
            p.samplerate = varargin{h_};
        case {'tref'}
            p.tref = varargin{h_};
        case {'doublingtreshold'}
            p.doublingtreshold = varargin{h_};
        case {'pas_.fire','f'}
            p.f = varargin{h_}; 
        otherwise
            error('Unknown value given with option ''% s'': ''% s''!',varargin{h_-1},any2str(varargin{h_}));
    end
end

% it is also possible to change phase coupling function further below
% SimulationLength = p.simulationlength;
FilterLength = p.filterlength;
% n = p.nodes;
sr = p.samplerate;
tref = p.tref;
DoublingThreshold = p.doublingthreshold;
Alpha = p.alpha;
Beta = p.beta;
filterType = str2func(p.filtertype);
f = p.f;

%set random seed
% rng(1234)

% for i = 2:SimulationLength
    i = 2;
        j=1;
        if pas_.phi < 1
            phi = pas_.phi+pas_.omega(j)/sr; % tooth saw -> increase (phase=amplitude)
            if phi>=1
                phi=1;
            end
        else 
            % fprintf('%d %d %f %d \n',i,j,phi,f(j))
            pas_.in=pas_.in+1;
            pas_.interactions(pas_.in,[1:3,11]) = [2 i j pas_.fire];
%             f(j)=mod(f(j)+1,2); % pas_.fire only on every other peak
            pas_.fire = mod(pas_.fire+1,2);
            phi = 0;
        end
%     end

    phiref = tref*pas_.omega*5e-4;
    
    if pas_.ii==12635
        
    end

%     for j=randperm(n) % go through all oscillators in random order.
%        fprintf('time %d - pas_.phi(%d,:)=[%f %f], node %d\n',i,i, pas_.phi(i,:),j) % (for n=2 only)
         if pas_.phi_ext >= 1        % if external node is at maximum
           % fprintf('pas_.phis %d %d %d \n',i,j,f(j))
           % f(j)=1;             % fire on every peak
            if f(j) == 0
                pas_.in=pas_.in+1;
                pas_.interactions(pas_.in,1:11) = [-1 i j 0 0 0 0 0 0 0 0];
            end
            if f(j) == 1 && pas_.fire == 1      % if node j pas_.fires and firing has already occured at this time.
                pas_.in=pas_.in+1;
                pas_.interactions(pas_.in,1:3) = [0 i j];
            end

            if f(j) == 1 && pas_.fire == 0       % if node j pas_.fires and firing has not already occured at this time.
%                 pas_.fire = 1; % set flag for pas_.fire event at time i
%                 phi = 1;   % j restrict self to max 1
                pas_.phi_ext = 1;
%                for k = setdiff(randperm(n),j,'stable') % k = all nodes but j
                    k = 1;
                    if phi > phiref % if the phase of node k is above phiref
                       % PHASE COUPLING FUNCTION
                        phi = phi+Alpha*-sin(phi*2*pi)*abs(sin(phi*2*pi));
                        if phi<0 % trim phase to [0 1] range
                            phi=0;
                        elseif phi>=1
                            phi=1;
                        end

                        if phi < phiref || phi > 1-phiref
                            pas_.E(end+1) = 0;
                            %E(i:SimulationLength,k) = 0; %Error is 0 if firing happens within phiref
                        else
                            pas_.E(end+1) = (1-cos(2*pi*phi))/2; %else, error is calculated like this
                            %E(i:SimulationLength,k) = (1-cos(2*pi*phi))/2; %else, error is calculated like this
                        end
%                         pas_.Ei{k}(length(pas_.Ei{k})+1) = i;


                        %this line makes oscillator frequencies converge tighter at the risk of instability. if disabled, phase coupling will ensure sync anyway
                        %if phi==1 && pas_.phi(i-1,k) < 0.95,pas_.domega(1,1)=pas_.domega(1,1)*1.01;end

                     if length(pas_.E) < FilterLength
                        pas_.C(end+1) = filterType(pas_.E(1:length(pas_.E)));
                     else
                        pas_.C(end+1) = filterType(pas_.E(length(pas_.E)-FilterLength+1:length(pas_.E)));
                     end
%
%                     % FREQUENCY COUPLING FUNCTION
%                      freqAdjConfidence = 1; %set to 0 to remove confidence from frequency update function
%
%                      if freqAdjConfidence
%                         pas_.domega(1,1) = ((pas_.omega * 2^(adapting(k)*Beta*sin(2*pi*phi)*pas_.C(i,k)))+prod(pas_.domega(1,1:2)))/(pas_.domega(1,2)+1);
%                      else
%                         pas_.domega(1,1) = ((pas_.omega * 2^(adapting(k)*Beta*sin(2*pi*phi)*0.2))+prod(pas_.domega(1,1:2)))/(pas_.domega(1,2)+1); %removing confidence measure
%                      end
%                         %USE FILTERED SIN(E) INSTEAD OF CONFIDENCE pas_.domega(1,1) = ((pas_.omega * 2^(adapting(k)*Beta*sin(2*pi*phi)*(filterType(E(length(E(:,k))-FilterLength+1:length(E(:,k)),k)))))+prod(pas_.domega(1,1:2)))/(pas_.domega(1,2)+1); % using sin(E) instead of

                        % FREQUENCY COUPLING FUNCTION

                        pas_.domega(1,1) = ((pas_.omega * 2^(Beta*-sin(2*pi*phi)*pas_.C(end)))+prod(pas_.domega(1,1:2)))/(pas_.domega(1,2)+1);
%                        pas_.domega(1,1) = ((pas_.omega * 2^(Beta*-sin(2*pi*phi)*(filterType(E(length(E(:,k))-FilterLength+1:length(E(:,k)),k)))))+prod(pas_.domega(1,1:2)))/(pas_.domega(1,2)+1);
                        pas_.domega(1,2) = pas_.domega(1,2)+1;
                        if pas_.domega(1,2) > DoublingThreshold
                            pas_.omega = 2*pas_.omega;
                            pas_.domega(1,:) = [0 0];
                            disp(['Someone''s stuck. Doubling frequency of fly ' num2str(k) ' Time ' num2str(i/sr) 's'])
                        end
                        pas_.in=pas_.in+1;
                        pas_.interactions(pas_.in,1:9) = [1 i j k 1 pas_.C(end) phi-pas_.phi Beta*sin(2*pi*phi) pas_.domega(1,1)];
                    else
                        pas_.in=pas_.in+1;
                        pas_.interactions(pas_.in,1:6) = [1 i j k 0 pas_.C(end)]; %col 5 == 0 means pas_.fire within tref
                    end
%                 end
            end

            if pas_.domega(j,2) > 0 % if there has been freq adjustments
                pas_.omega(j,1) = pas_.domega(j,1);    % Reachback pas_.firefly algorithm
                pas_.domega(j,:) = [0 0];
%                 if pas_.omega(j) == 0
%                         disp('a node has died... pas_.omega = zero')
%                 end
                pas_.in=pas_.in+1;
                pas_.interactions(pas_.in,1:10) = [-2 i j j 1 pas_.C(end) -1 0 0 pas_.omega(j)];
            end
%            [j;E(:,j)]

            if length(pas_.E) < FilterLength
                allzeroes = sum(pas_.E(1:length(pas_.E)));
                %allzeroes = filterType(pas_.E{j}(1:length(pas_.E{j})));
            else
                allzeroes = sum(pas_.E(length(pas_.E)-FilterLength+1:length(pas_.E)));
                %allzeroes = filterType(pas_.E{j}(length(pas_.E{j})-FilterLength+1:length(pas_.E{j})));
            end

                %tmp(length(tmp)+1)=allzeroes;
% 
%             if allzeroes == 0 && length(pas_.interactions(pas_.interactions(:,3)==j & pas_.interactions(:,11)==1,1)) > 8
%                 % length(pas_.interactions(pas_.interactions(:,3)==j,1))
%                 % disp('selfmean adjust')
%                 sorted.Omegas = sort(pas_.interactions(pas_.interactions(:,3)==j & pas_.interactions(:,11)==1,2),'descend')';
%                 
%                 pas_.omega(j) = -2000/mean(diff(sorted.Omegas(1:FilterLength)));
%                 pas_.in=pas_.in+1;
%                 pas_.interactions(pas_.in,1:10) = [-3 i j j 1 pas_.C(end) -1 0 0 pas_.omega(j)];
%             end
        end

%     end
for ii=1:size(pas_.interactions,1)
    if pas_.interactions(ii,10)==-Inf
    %    t =1; 
%     return
    end
end
%     for j=1:n
        if phi >= 1
            phi = 1;
        end
%             pas_.omegas(i,j) = pas_.omega(j);
%     end
% end

pas_.phi = phi;

if nargout > 0
    phaseOut = pas_.phi;
    if nargout > 1
        omegaOut = pas_.omega;
    end
end
% 
% figure(1)
% clf
% 
% % subplot 311
% % plot(C)
% % subplot 312
% % plot(pas_.E{1})
% 
% 
% % subplot(3,1,3)
% 
% yyaxis left
% plot(phaseOut,'b-'); hold on; plot(phi_ext,'r-'); xlabel('Time'); ylabel('Phase'); set(gca,'YColor','k')
% yyaxis right
% plot(omegaOut,'g-'); ylabel('Frequency'); set(gca,'YColor','g')
% xlim([2 SimulationLength])
% drawnow


end
