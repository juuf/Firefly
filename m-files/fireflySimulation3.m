function [phaseOut,omegaOut,phi_ext] = fireflySimulation3(varargin)
% Simulate decentralised pulse-coupled oscillators synchronising both phase
% and frequency to Harmonic Synchrony. 
% Developed by Kristian Nymoen, University of Oslo, 2013.
% Edit: Adapted for external input by Julian Fuhrer
%
% Usage:
% - [phaseOut] = fireflySimulation('argument',argumentValue)
% - [phaseOut omegaOut] = fireflySimulation('argument',argumentValue)
%
% Arguments: 
% 'alpha': Phase coupling constant (default 0.4)
% 'beta': Frequency coupling constant (default 0.4)
% 'SimulationLength': Number of samples to run the experiment (default 60000)
% 'FilterLength': Number of samples used in the running filter (default 8)
% 'FilterType': FilterType of the running filter: 'median' (default) or 'mean'
% 'Nodes': Number of oscillators (default 5)
% 'SampleRate': (default 1000)
% 'tref': refractory interval in ms (default 1)
% 'DoublingThreshold': Number of repeated negative phase jumps without 
%                      firing before frequency is doubled (default 20)
% 
% Example:
% [phaseOut,omegaOut] = fireflySimulation('alpha',0.2,'beta',0.1,'Nodes',4)
%   Starts an experiment with 4 oscillators, alpha value set to 0.2 and 
%   beta value set to 0.1. The output variable phaseOut contains the phi 
%   values of all involved oscillators.

global p pas_



p.simulationlength = 6e4; 
p.filterlength = 8;        
p.nodes = 5;       % number of nodes
p.samplerate = 1000;  
p.plotlength = 10000; 
p.tref = 5; %ms
p.doublingthreshold = 20;
% p.animspeed = 0.2;
% p.animate = 0;
p.alpha = 0.4;
p.beta = 0.4;
p.filtertype = 'median';
p.f=0;
   
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
       case {'fire','f'}
            p.f = varargin{h_};     
        otherwise
            error('Unknown value given with option ''% s'': ''% s''!',varargin{h_-1},any2str(varargin{h_}));
    end
end
    
% it is also possible to change phase coupling function further below
SimulationLength = p.simulationlength;
FilterLength = p.filterlength;
n = p.nodes;
sr = p.samplerate;
tref = p.tref;
DoublingThreshold = p.doublingthreshold;
Alpha = p.alpha;
Beta = p.beta;
f = p.f;
filterType = str2func(p.filtertype);


 

phi = zeros(SimulationLength,n);           
% phi(1,1:n) = rand(n,1);    % randomised initial phase
% phi(1,1:n) = [0 0.4 0.7];
temp2 = rand(2,1);
phi(1,1:n) = temp2(2); 

omegas = zeros(SimulationLength,n);
% omega = 2*2.^((rand(1,n)-0.5)*2);   % randomised initial frequency
% omega = [1.0635 2.6719 1.7153];
temp = 2*2.^((rand(2,1)-0.5)*2);
% omega = [2.6719];
omega = temp(1);

E = ones(SimulationLength,n);
% f = zeros(1,n);

domega = zeros(n,2); % domega variable for RFA

in = 0;
% pas_.interactions = zeros(1,11);
fire = p.fire; % input indicating if other node fired or not, i.e. p.fire is 1 or 0


phi_ext = zeros(SimulationLength,n);
% phi(1,1:n) = 0;
phi_ext(1) = temp2(2);
% omega_ext = 1.0635;%*ones(SimulationLength,n);
omega_ext= temp(2);

j=1; k=1; i=2; % can be replaced

% for i = 2:SimulationLength
    % external oscillator with phase phi_ext and constant frequency omega_ext
%     if phi_ext(i-1) < 1
%         phi_ext(i) = phi_ext(i-1)+omega_ext/sr; % tooth saw -> increase (phase=amplitude)
%         if phi_ext(i)>=1
%             phi_ext(i)=1;
%             f(1)=1;
% %             test(i) = f(1);
% %             f(1)=mod(f(1)+1,2);
%         else
%             f(1)=0;
% %             test(i) = f(1); 
%         end
%     end

    
    % increase the phase of each oscillator, else reset to zero (i.e. leave phi(i,j) at initial value = 0)
%     for j=1:n
        if phi(i-1,j) < 1
            phi(i,j) = phi(i-1,j)+omega(j)/sr; % tooth saw -> increase (phase=amplitude)
            if phi(i)>=1
                phi(i)=1;
            end
        else % -> phi at time i >=1 -> fire! 
            % fprintf('%d %d %f %d \n',i,j,phi(i,j),f(j))
            pas_.in=pas_.in+1;
            pas_.interactions(pas_.in,[1:3,11]) = [2 i j fire(j)];
%             f(j)=mod(f(j)+1,2); % fire only on every other peak
            fire = mod(fire+1,2);
%             fire =1;
        end
%     end
    
    phiref = tref*omega*5e-4;
    
%     for j=randperm(n) % go through all oscillators in random order.
%        fprintf('time %d - phi(%d,:)=[%f %f], node %d\n',i,i, phi(i,:),j) % (for n=2 only)
         if phi_ext(i) >= 1        % if external node is at maximum
           % fprintf('phis %d %d %d \n',i,j,f(j))
           % f(j)=1;             % fire on every peak
            if f(j) == 0
                pas_.in=pas_.in+1;
                pas_.interactions(pas_.in,1:11) = [-1 i j 0 0 0 0 0 0 0 0];
            end
            if f(j) == 1 && fire == 1      % if node j fires and firing has already occured at this time. 
                pas_.in=pas_.in+1;
                pas_.interactions(pas_.in,1:3) = [0 i j];
            end
           
            if f(j) == 0 && fire == 1       % if node j fires and firing has not already occured at this time.
%                 fire = 1; % set flag for fire event at time i 
%                 phi(i,j) = 1;   % j restrict self to max 1    
                phi_ext(i) = 1;
%                for k = setdiff(randperm(n),j,'stable') % k = all nodes but j                  
                    if phi(i,k) > phiref(k) % if the phase of node k is above phiref
                       % PHASE COUPLING FUNCTION
                        phi(i,k) = phi(i,k)+Alpha*-sin(phi(i,k)*2*pi)*abs(-sin(phi(i,k)*2*pi));                 
                        if phi(i,k)<0 % trim phase to [0 1] range
                            phi(i,k)=0;  
                        elseif phi(i,k)>=1
                            phi(i,k)=1;
                        end                       
                        if phi(i,k) > 1-phiref(k)
                            E(i:SimulationLength,k) = 0; % Error is 0 if firing happens within phiref
                        else
                            E(i:SimulationLength,k) = (1-cos(2*pi*phi(i,k)))/2; % else, error is calculated like this
                        end 
                        % FREQUENCY COUPLING FUNCTION
                        domega(k,1) = ((omega(k) * 2^(Beta*-sin(2*pi*phi(i,k))*(filterType(E(length(E(:,k))-FilterLength+1:length(E(:,k)),k)))))+prod(domega(k,1:2)))/(domega(k,2)+1);
                        domega(k,2) = domega(k,2)+1;
                        if domega(k,2) > DoublingThreshold
                            omega(k) = 2*omega(k);
                            domega(k,:) = [0 0];
                            disp(['Someone''s stuck. Doubling frequency of fly ' num2str(k) ' Time ' num2str(i/sr) 's'])
                        end
                        pas_.in=pas_.in+1;
                        pas_.interactions(pas_.in,1:9) = [1 i j k 1 filterType(E(length(E(:,k))-FilterLength+1:length(E(:,k)),k)) phi(i,k)-phi(i-1,k) Beta*sin(2*pi*phi(i,k)) domega(k,1)]; 
                    else
                        pas_.in=pas_.in+1;
                        pas_.interactions(pas_.in,1:6) = [1 i j k 0 filterType(E(length(E(:,k))-FilterLength+1:length(E(:,k)),k))]; %col 5 == 0 means fire within tref
                    end
%                 end
            end
      
            if domega(j,2) > 0 % if there has been freq adjustments
                omega(j) = domega(j,1);    % Reachback firefly algorithm
                domega(j,:) = [0 0];
%                 if omega(j) == 0
%                         disp('a node has died... omega = zero')
%                 end
                pas_.in=pas_.in+1;
                pas_.interactions(pas_.in,1:10) = [-2 i j j 1 filterType(E(length(E(:,j))-FilterLength+1:length(E(:,j)),j)) -1 0 0 omega(j)]; 
            end
%            [j;E(:,j)]
            if sum(E(length(E(:,j))-FilterLength+1:length(E(:,j)),j)) == 0 && length(pas_.interactions(pas_.interactions(:,3)==j & pas_.interactions(:,11)==1,1)) > 8
                % length(pas_.interactions(pas_.interactions(:,3)==j,1))
                % disp('selfmean adjust')
                sortedOmegas = sort(pas_.interactions(pas_.interactions(:,3)==j & pas_.interactions(:,11)==1,2),'descend')';
                omega(j) = -2000/mean(diff(sortedOmegas(1:FilterLength)));
                pas_.in=pas_.in+1;
                pas_.interactions(pas_.in,1:10) = [-3 i j j 1 filterType(E(length(E(:,j))-FilterLength+1:length(E(:,j)),j)) -1 0 0 omega(j)]; 
            end
        end
%     end

    for j=1:n
        if phi(i,j) >= 1
            phi(i,j) = 1;
        end
            omegas(i,j) = omega(j);
    end
% end



if nargout > 0
    phaseOut = phi;
    if nargout > 1
        omegaOut = omegas;
        phi_ext;
    end 
end


end