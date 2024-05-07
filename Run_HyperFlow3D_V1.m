%% Run code for torsion and analysis

clear
cla
DataFN = 'Hyperbolic_Flow3D_Test';        % File name for multi run data
Code_Name = 'Hyperbolic_Flow_V3D.m';       % Defining file name for the 
%%
row = 0; % initializing "row"value for data output
for q = 1:1
    for w = 1:14
        for run = 1:20
            for f = 1:1
                row = row+1;
                [q,w,run,f]
                cla
                close all

                eta = .1;                           % viscosity in pascal seconds
                kappa = 1E-15;                      % relaxed bend / unit length
                tau0 = 1E-15;                       % relaxed twist / unit length
                len = 4.0 + (w-1)*0.5;              % total length
                rad = 0.05.*len/5.0;                % Radius of the filament (um)
                
                GridNum = round(len*30+1,0);        % number of nodes
                ds = len./(GridNum-1);              % length of ds segment
                kT = .004114;                       % blotzmann
                A1 = 100.*kT;                       % Bend modulus while twisting about the e1 basis
                A2 = A1;                            % Bend modulus while twisting about the e2 basis
                Cc = A1;                            % Twist modulus while twisting about the e3 basis
                sig0 = .8E4.*((A1+A2)./(2*0.2));    % non-linear spring force along line connecting points
                sig0 = sig0.*ds;                    % redefine sig0 based on ds
                               
                dt = 3.0*10^-6.*(eta./0.1);       % time step
                ttime = .12;      % total time of simulation (used when no motor action.... (else use the final pulse of spiroplasma motor to set final time)
                samplet = 2E-4;   % how often to collect the data and plot
                
                gamx = -0.05E3;      % x-direction rate of hyperbolic fluid flow (sec^-1). negative is compressile
                gamy = -gamx;        % y-direction rate of hyperbolic fluid flow (sec^-1). negative is compressile
                gamz = 0;            % z-direction rate of hyperbolic fluid flow (sec^-1). negative is compressile   

                THF = 1;        % turns thermal force on(1) or off(0)
                RNGmod = run;   % either default or shuffle in order to initialize the random number generator. specifying a number sets 'seed' of Mersenne Twister RNG

                
                %%

                [info,compDataEnd,myRunTime]=Hyperbolic_Flow_V3D(gamx,gamy,gamz,eta,kappa,tau0,A1,A2,Cc,rad,THF,GridNum,dt,sig0,RNGmod,len,ttime,samplet);
                
                IterData1(row,:) = [compDataEnd];
                IterData2(row,:) = [info,run,myRunTime];

                
                
            end 
        end
    end
end
AggTable1 = table(IterData1(:,1),IterData1(:,2),IterData1(:,3),IterData1(:,4),IterData1(:,5),...
                  IterData1(:,6),IterData1(:,7),IterData1(:,8),IterData1(:,9),IterData1(:,10),...
                  'VariableNames',{'Time','LeeDot','ymax', 'zmax', 'Reff','bFact', 'epsDotT','Lee','R/L','muHat'});
writetable(AggTable1,strcat(DataFN,'_data.xlsx'),'Sheet','CompData');

AggTable2 = table(IterData2(:,1),IterData2(:,2),IterData2(:,3),IterData2(:,4),IterData2(:,5),IterData2(:,6),IterData2(:,7),...
                  IterData2(:,8),IterData2(:,9),IterData2(:,10),IterData2(:,11),IterData2(:,12),IterData2(:,13),...
                  IterData2(:,14),IterData2(:,15),IterData2(:,16),IterData2(:,17),IterData2(:,18),...
                  'VariableNames',{'eta','GammaX','GammY','GammaZ', 'OrigLength','radius','Ac1','Ac2','Cc','Version', ...
                  'Sigma','GN','dt','GoldSigma','1-min(t2tL/len)','FinalLength','run','RunTime(s)'});
writetable(AggTable2,strcat(DataFN,'_data.xlsx'),'Sheet','FlowData');

