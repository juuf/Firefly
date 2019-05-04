
function [phaseOut,omegaOut] = fireflySimulation(varargin)
%
% Simulate decentralised pulse-coupled oscillators synchronising both phase
% and frequency to Harmonic Synchrony. 
% Developed by Kristian Nymoen, University of Oslo, 2013.
% 
% Usage:
% - fireflySimulation
% - fireflySimulation('argument',argumentValue,'argument2',argumentValue2,...)
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
% 'PlotLength': Number of samples to plot (default 10000)
% 'tref': refractory interval in ms (default 1)
% 'DoublingThreshold': Number of repeated negative phase jumps without firing before frequency is doubled (default 20)
% 'PlotSpectrogram': 0 (default) or 1
% 'RandomSeed': Specify a random seed in order to analyse different parameters
% 'Animate': 0 (default) or 1
% 
% Examples: 
% fireflySimulation
%   - Starts an experiment with default values (no output...)
%
% fireflySimulation('Animate',1)
%   - Starts an experiment with default values (showing animation)
%
% fireflySimulation('Nodes',6)
%   - Starts an experiment with 6 oscillators
%
% fireflySimulation('alpha',0.2,'beta',0.1,'Nodes',4,'PlotSpectrogram',1)
%   - Starts an experiment with 4 oscillators, alpha value set to 0.2 and 
%   beta value set to 0.1. Plotting spectrogram for all nodes
%
% omegaOut = fireflySimulation('RandomSeed',534);
%   - Starts an experiment with default settings using the random seed 534.
%   The experiment results will be the same every time. The output variable
%   w contains a matrix of the frequency values over time
%
% [phaseOut,omegaOut] = fireflySimulation;
%   - The output variable phaseOut contains the phi values of all involved
%   oscillators.
%




   p = inputParser;
   p.KeepUnmatched = true;
   %default settings
    defaultSimulationLength = 60000; 
    defaultFilterLength = 8;        
    defaultNodes = 5;       % number of nodes
    defaultSampleRate = 1000;  
    defaultPlotLength = 10000; 
    defaulttref = 5; %ms
    defaultDoublingThreshold = 20;
    defaultSeed = 0;
    defaultAnimSpeed = 0.2;
    defaultListenProb = 0.98;
    defaultActivateThreshold = 0.05;
    
    defaultAlpha = 0.4;
    defaultBeta = 0.4;
    defaultFilterType = 'median';
  
%it is also possible to change phase coupling function further below

%initialisation
addOptional(p,'SimulationLength',defaultSimulationLength,@isscalar);
addOptional(p,'FilterLength',defaultFilterLength,@isscalar);
addOptional(p,'Nodes',defaultNodes,@isscalar);
addOptional(p,'SampleRate',defaultSampleRate,@isscalar);
addOptional(p,'PlotLength',defaultPlotLength,@isscalar);
addOptional(p,'tref',defaulttref,@isscalar);
addOptional(p,'DoublingThreshold',defaultDoublingThreshold,@isscalar);
addOptional(p,'alpha',defaultAlpha,@isscalar);
addOptional(p,'beta',defaultBeta,@isscalar);
addOptional(p,'PlotSpectrogram',0,@isscalar);
addOptional(p,'RandomSeed',defaultSeed,@isscalar);
addOptional(p,'AnimSpeed',defaultAnimSpeed,@isscalar);
addOptional(p,'FilterType',defaultFilterType);
addOptional(p,'Animate',0);
addOptional(p,'ListenProb',defaultListenProb,@isscalar);
addOptional(p,'ActivateThreshold',defaultActivateThreshold,@isscalar);

parse(p,varargin{:});

SimulationLength = p.Results.SimulationLength;
FilterLength = p.Results.FilterLength;
n = p.Results.Nodes;
sr = p.Results.SampleRate;
PlotLength = p.Results.PlotLength;
tref = p.Results.tref;
DoublingThreshold = p.Results.DoublingThreshold;
Alpha = p.Results.alpha;
Beta = p.Results.beta;
plotSpectrogram = p.Results.PlotSpectrogram;
randomSeed = p.Results.RandomSeed;
filterType = str2func(p.Results.FilterType);
animate = p.Results.Animate;
animSpeed = p.Results.AnimSpeed*200;
listenProb = p.Results.ListenProb;
activateThreshold = p.Results.ActivateThreshold;

if randomSeed ~= 0
	rng(randomSeed);
end

omegas = zeros(SimulationLength,n);

adapting = ones(1,n); 

tmp = [];

phi = zeros(SimulationLength,n);           
phi(1,1:n) = rand(n,1);    %randomised initial phase
%phi(1,1:n) = [0 0.4 0.7]; 

omega = 2*2.^((rand(1,n)-0.5)*2);   %randomised initial frequency
%omega = ones(1,n);
%omega = [1.0635 2.6719 1.7153];

E{1,n} = []; for i=1:n,E{1,i}=1;end
Ei{1,n} = []; for i=1:n,Ei{1,i}=1;end
f = zeros(1,n);
C = ones(SimulationLength,n);

domega = zeros(n,2); %domega variable for RFA

in = 0;


for i = 2:SimulationLength
    
    fire = 0; %clear flag for fire event at time i
    
    %increase the phase of each oscillator, else reset to zero (i.e. leave phi(i,j) at initial value = 0)
    for j=1:n
        if phi(i-1,j) < 1
            phi(i,j) = phi(i-1,j)+omega(j)/sr;
        else
            %fprintf('%d %d %f %d \n',i,j,phi(i,j),f(j))
            in=in+1;interactions(in,[1:3,11]) = [2 i j f(j)];
            f(j)=mod(f(j)+1,2); % fire only on every other peak
            phi(i,j) = 0;
        end
    end
    
    phiref = tref*omega*0.0005;
    
    for j=randperm(n) % go through all oscillators in random order.
       % fprintf('time %d - %f %f %f node %d\n',i, phi(i,:),j)

        if phi(i,j) >= 1        % if node j is at maximum
           % fprintf('phis %d %d %d \n',i,j,f(j))
           %f(j)=1;             % fire on every peak
            if f(j) == 0
                in=in+1;interactions(in,1:11) = [-1 i j 0 0 0 0 0 0 0 0];
            end
            if f(j) == 1 && fire == 1;       % if node j fires and firing has already occured at this time. 
                in=in+1;interactions(in,1:3) = [0 i j];
            end
           
            if f(j) == 1 && fire == 0;       % if node j fires and firing has not already occured at this time.
                fire = 1; %set flag for fire event at time i
                phi(i,j) = 1;   % j restrict self to max 1
               
                
                for k = setdiff(randperm(n),j,'stable') %k = all nodes but j

                    if phi(i,k) > phiref(k) && rand < listenProb %if the phase of node k is above phiref AND the fly hears the fire event

                        %CHOOSE A PHASE COUPLING FUNCTION
                        pf = 5;

                        switch pf
                            case 1 %standard excitatory: phi + alpha * phi
                               phi(i,k) = phi(i,k)+Alpha*phi(i,k)*adapting(k);                                  
                            case 2 %standard, bidirectional
                               phi(i,k) = phi(i,k)+Alpha*phi(i,k)*(phi(i,k)-0.5)/abs(phi(i,k)-0.5)*adapting(k); 
                            case 3 %sinewave * phase
                               phi(i,k) = phi(i,k)+Alpha*phi(i,k)*-sin(phi(i,k)*2*pi)*adapting(k);              
                            case 4 %sinewave 
                               phi(i,k) = phi(i,k)+Alpha*-sin(phi(i,k)*2*pi)*adapting(k);              
                            case 5 %bidirectional sinewave^2 
                               phi(i,k) = phi(i,k)+Alpha*-sin(phi(i,k)*2*pi)*abs(-sin(phi(i,k)*2*pi))*adapting(k); 
                            case 6 %bidirectional sinewave^2 * phase
                               phi(i,k) = phi(i,k)+Alpha*phi(i,k)*-sin(phi(i,k)*2*pi)*abs(-sin(phi(i,k)*2*pi))*adapting(k); 
                            case 7 %bidirectional coswave^2
                               phi(i,k) = phi(i,k)+Alpha*-cos(phi(i,k)*pi)*abs(-cos(phi(i,k)*pi))*adapting(k);
                            case 8
                               phi(i,k) = phi(i,k)+Alpha*(phi(i,k)-0.5)^3*adapting(k);
                        end
                        
                        if phi(i,k)<0%trim phase to [0 1] range
                            phi(i,k)=0;  
                        elseif phi(i,k)>=1
                            phi(i,k)=1;
                        end
                        
                        if phi(i,k) < phiref(k) || phi(i,k) > 1-phiref(k)
                            E{k}(length(E{k})+1) = 0;
                            %E(i:SimulationLength,k) = 0; %Error is 0 if firing happens within phiref
                        else
                            E{k}(length(E{k})+1) = (1-cos(2*pi*phi(i,k)))/2; %else, error is calculated like this
                            %E(i:SimulationLength,k) = (1-cos(2*pi*phi(i,k)))/2; %else, error is calculated like this
                        end
                        Ei{k}(length(Ei{k})+1) = i;

                        
                        %this line makes oscillator frequencies converge tighter at the risk of instability. if disabled, phase coupling will ensure sync anyway
                        %if phi(i,k)==1 && phi(i-1,k) < 0.95,domega(k,1)=domega(k,1)*1.01;end

                     if length(E{k}) < FilterLength
                       C(i:SimulationLength,k) = filterType(E{k}(1:length(E{k})));
                     else
                        C(i:SimulationLength,k) = filterType(E{k}(length(E{k})-FilterLength+1:length(E{k})));
                     end
                    %C(i:SimulationLength,k) = E{k}(length(E{k}));
                     
                    
                    % FREQUENCY COUPLING FUNCTION
                     freqAdjConfidence = 1; %set to 0 to remove confidence from frequency update function
                     
                     if freqAdjConfidence
                        domega(k,1) = ((omega(k) * 2^(adapting(k)*Beta*sin(2*pi*phi(i,k))*C(i,k)))+prod(domega(k,1:2)))/(domega(k,2)+1);
                     else
                        domega(k,1) = ((omega(k) * 2^(adapting(k)*Beta*sin(2*pi*phi(i,k))*0.2))+prod(domega(k,1:2)))/(domega(k,2)+1); %removing confidence measure
                     end
                        %USE FILTERED SIN(E) INSTEAD OF CONFIDENCE domega(k,1) = ((omega(k) * 2^(adapting(k)*Beta*sin(2*pi*phi(i,k))*(filterType(E(length(E(:,k))-FilterLength+1:length(E(:,k)),k)))))+prod(domega(k,1:2)))/(domega(k,2)+1); % using sin(E) instead of 

                        

                        domega(k,2) = domega(k,2)+1;
                        if domega(k,2) > DoublingThreshold
                            omega(k) = 2*omega(k);
                            domega(k,:) = [0 0];
                            disp(['Someone''s stuck. Doubling frequency of fly ' num2str(k) ' Time ' num2str(i/sr) 's'])
                        end
                        in=in+1;interactions(in,1:9) = [1 i j k 1 C(i,k) phi(i,k)-phi(i-1,k) Beta*sin(2*pi*phi(i,k)) domega(k,1)]; 
                        %in=in+1;interactions(in,1:9) = [1 i j k 1 filterType(E(length(E(:,k))-FilterLength+1:length(E(:,k)),k)) phi(i,k)-phi(i-1,k) Beta*sin(2*pi*phi(i,k)) domega(k,1)]; 
                    else
                        in=in+1;interactions(in,1:6) = [1 i j k 0 C(i,k)]; %col 5 == 0 means fire within tref
                        %in=in+1;interactions(in,1:6) = [1 i j k 0 filterType(E(length(E(:,k))-FilterLength+1:length(E(:,k)),k))]; %col 5 == 0 means fire within tref
                    end
                end
            end
            

            
            if domega(j,2) > 0 %if there has been freq adjustments
                omega(j) = domega(j,1);    %Reachback firefly algorithm
                domega(j,:) = [0 0];
                if omega(j) == 0
                        disp('a node has died... omega = zero')
                end
                in=in+1;interactions(in,1:10) = [-2 i j j 1 C(i,j) -1 0 0 omega(j)]; 
            end
%            [j;E(:,j)]

                if length(E{j}) < FilterLength
                    allzeroes = sum(E{j}(1:length(E{j})));
                    %allzeroes = filterType(E{j}(1:length(E{j})));
                else
                    allzeroes = sum(E{j}(length(E{j})-FilterLength+1:length(E{j})));
                    %allzeroes = filterType(E{j}(length(E{j})-FilterLength+1:length(E{j})));
                end

                %tmp(length(tmp)+1)=allzeroes;
            
            if allzeroes == 0 && length(interactions(interactions(:,3)==j & interactions(:,11)==1,1)) > 8
                %length(interactions(interactions(:,3)==j,1))
                disp(['selfmean adjust ' num2str(j) ' ' num2str(i)])
                sortedOmegas = sort(interactions(interactions(:,3)==j & interactions(:,11)==1,2),'descend')';
                omega(j) = -2000/mean(diff(sortedOmegas(1:FilterLength)));
                in=in+1;interactions(in,1:10) = [-3 i j j 1 C(i,j) -1 0 0 omega(j)]; 
            end
          
            
        end
        
        
        
        if sin(pi*E{j}(end)) > activateThreshold 
        %if sin(pi*C(i,j)) > activateThreshold 
            adapting(j) = 1; %Turning on the sync mechanism when confidence is above some threshold
        else
            adapting(j) = 0; %Turning off the sync mechanism when confidence is below some threshold
        end
    end

    for j=1:n
        if phi(i,j) >= 1
            phi(i,j) = 1;
        end
            omegas(i,j) = omega(j);
    end
end

%figure,plot(tmp)

if plotSpectrogram == 1
    

figure(1),clf(1)
    subplot(5,4,1:2)
        plot((1:PlotLength)/sr,phi(1:PlotLength,:)),title('phase: beginning')

    subplot(5,4,3:4)
        plot((SimulationLength-PlotLength:SimulationLength)/sr,phi(SimulationLength-PlotLength:SimulationLength,:)),title('phase end')

    subplot(5,4,[5 10])
        plot((1:SimulationLength)/sr,phi(1:SimulationLength,:)+1.05*cumsum(ones(n,SimulationLength))'-2),title('phase')
        xlim([0 PlotLength]/sr)
        hold on
        scatter(interactions(interactions(:,11)==1,2)/1000,(interactions(interactions(:,11)==1,3)-1)*1.05,'k.')

    subplot(5,4,[7 12])
        plot((1:SimulationLength)/sr,filter([1 -1],1,phi(1:SimulationLength,:))+1.15*cumsum(ones(n,SimulationLength))'-2),title('delta phase')
        xlim([0 PlotLength]/sr)

    subplot(5,4,[13 18])
        plot((1:SimulationLength)/1000,omegas),title('freq')

    subplot(5,4,[15 20])
        tmp = filter([1 -1],1,omegas);tmp(1:n,:)=[];
        plot((n+1:SimulationLength)/1000,tmp),title('delta freq')


%     figure('Name','Spectrograms of firefly activity')
% 	for i = 1:n
% 	subplot(n,1,i)
% 	spectrogram(sin(2*pi*phi(:,i)),4096,2048,10000,sr,'yaxis');ylim([0 8])
%     end
    
   
   
    
 
%     figure('Name','Confidence')
%     for i = 1:n
%     tmpplotvar = interactions(interactions(:,4)==i,[2,6]);
% 	subplot(n,1,i)
%     plot(tmpplotvar(:,1),tmpplotvar(:,2))
%     clear tmpplotvar;
%     end    
     

%     figure('Name','Instantaneous Error')
% 	for i = 1:n
% 	subplot(n,1,i)
% 	plot(E(:,i))
%     end    
    
    

    figure(2)
	for i = 1:n
	subplot(n,1,i)
	%plot((1:size(C(:,i)))/sr,sin(pi*C(:,i)))
    plot(Ei{k}/sr,sin(pi*E{k}),'.')
    ylim([0 1])
    xlim([0 SimulationLength/sr])
    end    
    set(gcf,'Name','Distance from 1 or 0')
    

    figure(3)
	for i = 1:n
	subplot(n,1,i)
	plot((1:size(C(:,i)))/sr,sin(pi*C(:,i)))
    %plot(Ei{k}/sr,sin(pi*E{k}),'.')
    ylim([0 1])
    xlim([0 SimulationLength/sr])
    end    
    set(gcf,'Name','C')
end

if nargout > 0
    phaseOut = phi;
    if nargout > 1 || animate == 1
        omegaOut = omegas;
    end 
end

if animate == 1
    subplot(6,2,1:4);
    plot((1:length(phaseOut(:,1)))/sr,phaseOut),xlabel('time (s)'),ylabel('phase')
    l1 = line([0 0],[0 1],'LineWidth',2,'color',[0 0 0]);
    k3 = subplot(6,2,[6 8]);
    plot((1:length(phaseOut(:,1)))/sr,phaseOut),xlabel('time (s)'),ylabel('phase'),hold on
    l3 = line([0 0],[0 1],'LineWidth',2,'color',[0 0 0]);
    
    k4 = subplot(6,2,[10 12]);
    plot((1:length(omegaOut(:,1)))/sr,omegaOut),xlabel('time (s)'),ylabel('frequency'),hold on
    l4 = line([0 0],[0 4],'LineWidth',2,'color',[0 0 0]);
    
    k2 = subplot(6,2,5:2:11);
    for i = 1+5*animSpeed:animSpeed:length(phaseOut(:,1))
        tic
        xlim(k3,[i-sr,i+sr]/sr)
        xlim(k4,[i-sr,i+sr]/sr)
        polar(pi, 1);hold on
    p1 = polar(k2,phaseOut(i-5*animSpeed:i,:)*2*pi,1./omegaOut(i-5*animSpeed:i,:),'-');
    line([0 3],[0 0],'LineWidth',3,'color',[0 0 0]);
    p2 = polar(k2,phaseOut(i:i+1,:)*2*pi,1./omegaOut(i:i+1,:),'.');
    set(l1,'xdata',[i/sr i/sr])
%    set(l5,'xdata',[i/sr i/sr])
    set(l3,'xdata',[i/sr i/sr])
    set(l4,'xdata',[i/sr i/sr])
    set(p1,'LineWidth',1)
    for q = 1:n
        r = interactions(interactions(:,1)==2 & interactions(:,2)>i & interactions(:,3)==q,11);
        if ~isempty(r)
            if r(1) == 1 
                set(p2(q),'MarkerSize',5,'Marker','o')
                %p2 = polar(k2,phaseOut(i,q)*2*pi,1./omegaOut(i,q),'o');
            else
                set(p2(q),'MarkerSize',18)
            end
        end
    end
    hold off
    drawnow
    pause(0.2-toc)
    end
end

end