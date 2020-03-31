function output = simulateOutbreakSEIcIscR(R0t,rho=c(rep(0.4,4),rep(0.8,12)), R0tpostoutbreak = 1.5,dateEndIntenseIntervention, #date we begin relaxing intense intervention 
                                    pWorkOpen = c(0.1,0.25,0.5,0.9), # pWorkOpen: proportion of the work force that is working (will be time-varying)
                                    dateStartSchoolClosure = as.Date('2020-01-15') , # cause winter term break 
                                    dateStartIntenseIntervention = as.Date('2020-01-23') , #Intense intervention: starts at Wuhan Lockdown
                                    dateStart = as.Date('2019-11-01'),POP = wuhanpop,numWeekStagger=c(2,4,6),pInfected=0.0002,durInf = 7,contacts_china=contacts)

  % debug dateStartIntenseIntervention = as.Date('2020-01-23')  
  % debug dateEndIntenseIntervention = as.Date('2020-03-01')
  % debug R0est = rep(2,3660) 
  % debug rho = rep(0.8,3660) 
  % debug pWorkOpen =  c(0.1,0.25,0.5,1)
  
  
  % Load population information
  % pop = loadPopInfo(POP = wuhanpop)
  N = sum(POP(:,2));
  p_age = POP(:,3);
  N_age = N*p_age ;
  
  
 % Specify epi info
  durLat = 6.4;   	                                             % Mean latent period (days) from Backer, et al (2020)
  durInf = durInf;                                               % Mean duration of infectiousness (days)
  gamma = 1-exp(-1/durInf);                                      % removal rate
  alpha = 1-exp(-1/durLat);                                      % infection rate
  dt = 1;                                                        % Time step (days)
  tmax = 428;                                                    % Time horizon (days) 366 days in 2020 cause of leap year
  numSteps = tmax/dt;  	                                         % Total number of simulation time steps
  % dateStart = as.Date('2019-12-01')                            % included as a function argument 
  dateEnd = dateStart+(tmax-1);
  dateStartCNY = datetime(2020,01,25); 
  dateEndCNY = datetime(2020,01,31) ;
  
  
  % Declare the state variables and related variables:
  % The values of these variables change over time
  [S,E,Isc,Ic,R] = deal(zeros(numSteps,length(p_age)));
  [lambda, incidence,reported, cumulativeIncidence] = deal(zeros(numSteps,length(p_age)));
  time = zeros(1,numSteps);
  
  % Initialise the time-dependent variables, i.e. setting the values of the variables at time 0
  E(1,:) = 0 ;
  Ic(1,:) =  pInfected*sum(N_age)/16;%rpois(length(N_age),lambda = pInfected*sum(N_age)/16)  % 100 # Assign 100 infected person in each age group (TODO RELAX?)
  Isc(1,:)=0;
  R(1,:) = 0 ;
  S(1,:) = N_age-E(1,:)-Ic(1,:)-Isc(1,:)-R(1,:);
  incidence(1,:) = 0;
  subclinical(1,:) = 0;
  time(1) = 0;
  
 
  % INTERVENTIONS 
  % School closed 2020-02-10, lockdown (intense intervention) started 2020-01-23, end of intense intervention: user-specified 
  % note that intense intervention is time-varying control by pWorkOpen: proportion of the work force that is working
  % debug pWorkOpen = c(0.1,0.25,0.5,1)
  tStartSchoolClosure = hours(dateStartSchoolClosure - dateStart)/24+1;
  tStartIntenseIntervention = hours(dateStartIntenseIntervention - dateStart)/24+1; % for pw = 0.1
  tEndIntenseIntervention = hours(dateEndIntenseIntervention - dateStart)/24+1; % for pw = 0.1
  tRelaxIntervention1 = tEndIntenseIntervention + numWeekStagger(1)*7;           % for pw = 0.25
  tRelaxIntervention2 = tEndIntenseIntervention + numWeekStagger(2)*7;           % for pw = 0.5
  tRelaxIntervention3 = tEndIntenseIntervention + numWeekStagger(3)*7;           % for pw = 1
  % tStartEndClosure = as.vector(dateEndSchoolClosure - dateStart)+1;
  
  
  pwork = ones(1,numSteps);
  p1 = ones(1,tStartIntenseIntervention-0);
  p2 = pWorkOpen(1)*ones(1,tEndIntenseIntervention-tStartIntenseIntervention);
  p3 = pWorkOpen(2)*ones(1,tRelaxIntervention1-tEndIntenseIntervention);
  p4 = pWorkOpen(3)*ones(1,tRelaxIntervention2-tRelaxIntervention1);
  p5 = pWorkOpen(4)*ones(1,tRelaxIntervention3-tRelaxIntervention2);
  pwork(1:tRelaxIntervention3) = [ p1 p2 p3 p4 p5];
  
 R0tpostoutbreak = R0t;
  constraintsIntervention_base = loadInterventions(0.5,"base");
 beta = getbeta(R0t, constraintsIntervention_base,gamma,p_age,1,contacts_china);
  if pWorkOpen(2)<1 
      beta_postfirstwave = getbeta(R0tpostoutbreak, constraintsIntervention_base,gamma,p_age,1,contacts_china);
  else
    if pWorkOpen(2)>=1 
        beta_postfirstwave = beta;
    end
  end
  for stepIndex = 1: (numSteps-1)
    
     % load plausible intervetions 
    %constraintsIntervention = loadInterventions(p_workopen = pwork[stepIndex])
    
    % Age- and location-specific contact rates for the given interventions 
    
    % I0: before school winter break intervention period, use base-case
    if time(stepIndex) < tStartSchoolClosure  
      INTERVENTION = "base"
      CONSTRAINT = loadInterventions(0.5,INTERVENTION);%constraintsIntervention$base
    end
    % I1:  When school winter break but before lockdown period, use 'schcloseonly'
    if time(stepIndex) >= tStartSchoolClosure & time(stepIndex) < tStartIntenseIntervention 
      INTERVENTION = "schcloseonly"   
      CONSTRAINT = loadInterventions(0.5,INTERVENTION);
    end  
     % I2:  Wuhan lockdown period, use 'schcloseonly'
    if time(stepIndex) >= tStartIntenseIntervention & time(stepIndex) < tRelaxIntervention3 
      INTERVENTION = "schcloseworkplacedist" ;
      CONSTRAINT = loadInterventions(pwork(stepIndex),INTERVENTION);
    end
    if time(stepIndex) >= tRelaxIntervention3  
        INTERVENTION = "postoutbreak" ;
        CONSTRAINT = loadInterventions(0.5,INTERVENTION);
    end
     
    
     C = CONSTRAINT.home*contacts_china.home + CONSTRAINT.work*contacts_china.work+ CONSTRAINT.school*contacts_china.school+ CONSTRAINT.others*contacts_china.others;
    
       
   % calculate the force of infection
    
    % beta = getbeta(R0t = R0t[stepIndex],constraints = constraintsIntervention$base,gamma = gamma,p_age = pop$p_age)

    if time(stepIndex)  < tEndIntenseIntervention+0
        lambda(stepIndex,:) = beta*C*I(stepIndex,:)'/N_age;
    else 
        lambda(stepIndex,:) = beta_postfirstwave*C*I(stepIndex,:)'/N_age;
    end
    % lambda[stepIndex,] = as.numeric(beta)*(as.matrix(C)%*%as.matrix(I[stepIndex,]/N_age));
    % calculate the number of infections and recoveries between time t and t+dt
    
    numStoE   = lambda(stepIndex,:).*S(stepIndex,:)*dt;                  # S to E
    numEtoIc  = alpha*rho*E(stepIndex,:)*dt;                          # E to Ic
    numEtoIsc = alpha*(1-rho)*E(stepIndex,:)*dt;                         # E to Isc
    numIctoR  = gamma*Ic(stepIndex,:)*dt;                              # Ic to R
    numIsctoR = gamma*Isc(stepIndex,:)*dt;                          # Isc to R
    
    % Difference equations 
    S(stepIndex+1,:)   = S(stepIndex,:)-numStoE;
    E(stepIndex+1,:)   = E(stepIndex,:)+numStoE-numEtoIc-numEtoIsc;
    Ic(stepIndex+1,:)  = Ic(stepIndex,:)+numEtoIc-numIctoR;
    Isc(stepIndex+1,:) = Ic(stepIndex,:)+numEtoIsc-numIsctoR;
    R(stepIndex+1,:)   = R(stepIndex,:)+numIctoR+numIsctoR;
    
    incidence(stepIndex+1,:) = numEtoIc/dt;
    subclinical(stepIndex+1,:) = numEtoIsc/dt;
    time(stepIndex+1) = time(stepIndex)+dt;

    
  end
  output=struct;
  output.S= S;
  output.E= E;
  output.Ic= Ic;
  output.Isc= Isc;
  output.R= R;
  output.time=time;
  output.lambda=lambda;
  output.incidence = incidence;
  output.N_age= N_age;
  output.R0t = R0t;
  output.subclinical=subclinical;
  output.dateStart = dateStart;
  output.dateEnd = dateEnd;
  output.dateStartIntenseIntervention = dateStartIntenseIntervention;
  output.dateEndIntenseIntervention = dateEndIntenseIntervention;
  output.dateStartSchoolClosure = dateStartSchoolClosure;
  output.dateStartCNY = dateStartCNY;
  output.dateEndCNY = dateEndCNY;
  
  end
