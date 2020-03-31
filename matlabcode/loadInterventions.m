function constraints = loadInterventions(p_workopen,type)
diagonal = ones(1,16);
d2 = diagonal;
d2(1:4)=0.5;

% constraints under a DO-NOTHING scenario
switch type
    case "base"
        base=struct;
        base.home = diag(diagonal);
        base.work = diag(diagonal);
        base.school = diag(diagonal);
        base.others = diag(diagonal);
        constraints=base;
        
        
        % Wuhan's lockdown--assume from XX Jan to XX Feb
    case "wuhanlockdown"
        wuhanlockdown=struct;
        wuhanlockdown.home = diag(diagonal);
        wuhanlockdown.work = diag(diagonal*0.1);
        wuhanlockdown.school = diag(diagonal*0);
        wuhanlockdown.others = diag(diagonal*0.1);
        constraints = wuhanlockdown;
        
        
        % constraints under school closure + some social distancing for school-age going children but 100% workplace
    case "schcloseonly"
        schcloseonly=struct;
        schcloseonly.home = diag(diagonal);
        schcloseonly.work = diag(diagonal);
        schcloseonly.school = diag(diagonal*0);
        schcloseonly.others = diag(d2);
        constraints =schcloseonly;
        
        % constraints under work place distancing only (MAYBE UNREALISTIC, should close schools too)
    case "workplacedistonly"
        workplacedistonly=struct;
        workplacedistonly.home = diag(diagonal);
        workplacedistonly.work = diag(diagonal*0.5);
        workplacedistonly.school = diag(diagonal);
        workplacedistonly.others = diag(diagonal*0.1);
         constraints =workplacedistonly;
        % constraints under work place distancing + schoolclosure
    case "schcloseworkplacedist"
        schcloseworkplacedist=struct;
        schcloseworkplacedist.home = diag(diagonal);
        schcloseworkplacedist.work = diag(diagonal*p_workopen);
        schcloseworkplacedist.school = diag(diagonal*0);
        schcloseworkplacedist.others = diag(diagonal*0.1);
        constraints =schcloseworkplacedist;
        % Post Outbeak, people still cautious
    case "postoutbreak"
        postoutbreak=struct;
        postoutbreak.home = diag(diagonal);
        postoutbreak.work = diag(diagonal);
        postoutbreak.school = diag(diagonal);
        postoutbreak.others = diag(diagonal);
        constraints =postoutbreak;
end
        
end