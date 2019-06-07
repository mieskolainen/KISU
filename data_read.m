%% Read ALICE HEPDATA datasets
%
% Call this from data_main.m
%
% mikael.mieskolainen@cern.ch, 2019

tic;

% Save INEL and NSD
INEL = {};
NSD  = {};

for t = 1:45 % Over all datafiles

    filename = sprintf('./HEPDATA/HEPData-ins1614477-v1-csv/Table%d.csv', t);

    fileid = fopen(filename,'r');
    for i = 1:11
    a   = fgets(fileid);
    a = a(1:end-1); % Remove \n

    % Find out event topology class
    if (i == 3)
        if     (contains(a,'INEL'))
            type = 'INEL';
        elseif (contains(a,'NSD'))
            type = 'NSD';
        else
            type = ''; % Unknown
        end
    end
    
    % Pseudorapidity interval
    if (i == 10)
        eta = sscanf(sscanf(a,'%s'), '#:ETARAP(P=5),,,%f-%f');
    end

    % CMS-energy
    if (i == 11)
        sqrts = sscanf(sscanf(a,'%s'), '#:SQRT(S)[GeV],,,%f');
        sqrts = sqrts/1000; % GeV to TeV
    end
    end
    fclose(fileid);
    
    % ------------------------------------------------------------------------
    % Read the data

    % Skip asymmetric
    %{
    if (eta(1) ~= -eta(2))
        continue;
    end
    %}

    % Read data
    R1 = 12;
    C1 = 0;
    M = csvread(filename, R1, C1);

    % x and y values
    N     = M(1,1); % Is N >= 1 or N >= 0 class
    n     = M(:,1); 
    g     = M(:,4); 
    g_pos = abs(M(:,5)); % + error
    g_neg = abs(M(:,6)); % - error
    
    % Set zero bin if not in data
    if (n(1) == 1)
        n = [0; n];
        g = [0; g];
        
        g_pos = [0; g_pos];
        g_neg = [0; g_neg];
    end
    
    % Skip zero bin
    if (SKIP0BIN && n(1) == 0)
        n = n(2:end);
        g = g(2:end);
        g_pos  = g_pos(2:end);
        g_neg  = g_neg(2:end);
        
        binsum = sum(g);
        g      = g / binsum; % Renormalize
        g_pos  = g_pos / binsum;
        g_neg  = g_neg / binsum;
    end
    
    % Number of triggered events from papers (here modulo efficiency corrections)
    if (sqrts == 0.9)
        EVENTS = 7.4E6;
    end
    if (sqrts == 7.0)
        EVENTS = 61E6;
    end
    if (sqrts == 8.0)
        EVENTS = 26E6;
    end
    
    % Save data
    if     (strcmp(type, 'INEL')) % true if identical
        INEL{end+1} = struct('type','INEL','N',N,'sqrts',sqrts,'eta',eta,'n',n,'g',g,'g_pos',g_pos,'g_neg',g_neg,'EVENTS',EVENTS);
    elseif (strcmp(type, 'NSD'))  % true if identical
        NSD{end+1}  = struct('type','NSD', 'N',N,'sqrts',sqrts,'eta',eta,'n',n,'g',g,'g_pos',g_pos,'g_neg',g_neg,'EVENTS',EVENTS);
    end
end
