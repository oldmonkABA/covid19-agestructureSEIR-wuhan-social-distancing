clear

wuhanpop = csvread('../data/wuhanpop.csv',1);
load '../data/contacts_china.mat' 
results = getbeta(1,1,1,wuhanpop(:,3),1,contacts);
POP=wuhanpop;
 N = sum(POP(:,2));
  p_age = POP(:,3);
  N_age = N*p_age ;
CHECKMODEL  = true;

if CHECKMODEL

  % Quick checks: Simulate an outbreak for sanity checks
  %set.seed(666)

  
  % test for an R0 value of 2.2
  R0est = 2.2;
  % R0est = sample(x = r0posterior,size = 100)
  
  nsim = 2;
  %epi_doNothing = vector('list',nsim)
  %epi_base = vector('list',nsim)
  %epi_march = vector('list',nsim)
  %epi_april = vector('list',nsim)
  
                           
                            
                            
                            
  R0tpostoutbreak = 1.5;
  dateStartSchoolClosure_nothing = datetime(2019,12,01);
  dateStartSchoolClosure_default = datetime(2020,01,15); % cause winter term break 
  dateStartIntenseIntervention_nothing = datetime(2019,12,01);
  dateStartIntenseIntervention_default = datetime(2020,01,23); %Intense intervention: starts at Wuhan Lockdown
  dateEndIntenseIntervention_nothing = datetime(2019,12,01);
  dateEndIntenseIntervention_base= datetime(2020,1,27);
  dateEndIntenseIntervention_march=datetime(2020,3,1);
  dateEndIntenseIntervention_april=datetime(2020,4,1);
  dateStart = datetime(2019,11,01);
  pWorkOpen_default = [0.1,0.25,0.5,0.9];
  POP = wuhanpop;
  rho = ones(1,3660)*0.5;
  pWorkOpen_nothing = ones(1,4);%(1,1,1,1)
  pWorkOpen_base = [0.5,1,1,1];
  numWeekStagger_base = [0,0,0,0];
  numWeekStagger_nothing = zeros(1,4);%(0,0,0,0)
  numWeekStagger_default=[2,4,6];
  pInfected=0.0002;
  durInf = 7;
  contacts_china=contacts;
  
  
%   R0t,rho, R0tpostoutbreak = 1.5,dateEndIntenseIntervention, %date we begin relaxing intense intervention 
%                             pWorkOpen = c(0.1,0.25,0.5,0.9), % pWorkOpen: proportion of the work force that is working (will be time-varying)
%                             dateStartSchoolClosure = as.Date('2020-01-15') , % cause winter term break 
%                             dateStartIntenseIntervention = as.Date('2020-01-23') , %Intense intervention: starts at Wuhan Lockdown
%                             dateStart = as.Date('2019-11-01'),POP = wuhanpop,numWeekStagger=c(2,4,6),pInfected=0.0002,durInf = 7,contacts_china=contacts
  for sim = 1:nsim
  
    epi_doNothing(sim) = simulateOutbreakSEIR(R0est,rho,R0tpostoutbreak,dateEndIntenseIntervention_nothing,...
                                                    pWorkOpen_nothing,dateStartSchoolClosure_nothing,....
                                                    dateStartIntenseIntervention_nothing,...
                                                    dateStart,POP,numWeekStagger_nothing,...
                                                    pInfected,durInf,contacts_china);
    epi_base(sim) = simulateOutbreakSEIR(R0est,rho,R0tpostoutbreak,dateEndIntenseIntervention_base,...
                                               pWorkOpen_base,dateStartSchoolClosure_default,....
                                                    dateStartIntenseIntervention_default,...
                                                    dateStart,...
                                                    POP,numWeekStagger_base,pInfected,durInf,contacts_china)
    epi_march(sim) = simulateOutbreakSEIR(R0est,rho,R0tpostoutbreak,dateEndIntenseIntervention_march,...
                                                    pWorkOpen_default,dateStartSchoolClosure_default,....
                                                    dateStartIntenseIntervention_default,...
                                                    dateStart,...
                                                    POP,numWeekStagger_default,pInfected,durInf,contacts_china)
    epi_april(sim) = simulateOutbreakSEIR(R0est,rho,R0tpostoutbreak,dateEndIntenseIntervention_april,...
                                                    pWorkOpen_default,dateStartSchoolClosure_default,....
                                                    dateStartIntenseIntervention_default,...
                                                    dateStart,...
                                                    POP,numWeekStagger_default,pInfected,durInf,contacts_china)
  
  end
  
%   par(mfrow=c(2,1))
%   
%   # incidence over time
  agegp =3;
  figure(1);
  subplot(2,1,1);  
  plot(epi_doNothing(1).time, epi_doNothing(1).incidence(:,agegp),'linewidth',1,'Color','k');
  title("Incidence for age 10-15");
  xlabel("Time(days)");
  ylabel("Daily no. of infections");
  hold on;
  plot(epi_base(1).time,epi_base(1).incidence(:,agegp),'linewidth',1,'Color',[0 0 0]+0.05*10);
  plot(epi_march(1).time,epi_march(1).incidence(:,agegp),'linewidth',1,'Color','b');
  plot(epi_april(1).time,epi_april(1).incidence(:,agegp),'linewidth',1,'Color','r');
   legend("Do Nothing", "Base","Lockdown->March","Lockdown->April");

  subplot(2,1,2);
  % cumulative incidence over time
  plot(epi_doNothing(1).time, (epi_doNothing(1).N_age(agegp)-epi_doNothing(1).S(:,agegp))./epi_doNothing(1).N_age(agegp),'linewidth',1,'Color','k'); 
  title("Cum incidence for age 10-15");
  xlabel("Time(days)");
  ylabel("Cum incidence");
  ylim = [0,1];
  hold;
plot(epi_base(1).time, (epi_base(1).N_age(agegp)-epi_base(1).S(:,agegp))./epi_base(1).N_age(agegp),'linewidth',1,'Color',[0 0 0]+0.05*10);
plot(epi_march(1).time, (epi_march(1).N_age(agegp)-epi_march(1).S(:,agegp))./epi_march(1).N_age(agegp),'linewidth',1,'Color','b');
plot(epi_april(1).time, (epi_april(1).N_age(agegp)-epi_april(1).S(:,agegp))./epi_april(1).N_age(agegp),'linewidth',1,'Color','r');
  legend("Do Nothing", "Base","Lockdown->March","Lockdown->April");
         
  agegp =13;
  figure(2);
  subplot(2,1,1);  
  plot(epi_doNothing(1).time, epi_doNothing(1).incidence(:,agegp),'linewidth',1,'Color','k');
  title("Incidence for age 60-65");
  xlabel("Time(days)");
  ylabel("Daily no. of infections");
  hold on;
  plot(epi_base(1).time,epi_base(1).incidence(:,agegp),'linewidth',1,'Color',[0 0 0]+0.05*10);
  plot(epi_march(1).time,epi_march(1).incidence(:,agegp),'linewidth',1,'Color','b');
  plot(epi_april(1).time,epi_april(1).incidence(:,agegp),'linewidth',1,'Color','r');
   legend("Do Nothing", "Base","Lockdown->March","Lockdown->April");

  subplot(2,1,2);
  % cumulative incidence over time
  plot(epi_doNothing(1).time, (epi_doNothing(1).N_age(agegp)-epi_doNothing(1).S(:,agegp))./epi_doNothing(1).N_age(agegp),'linewidth',1,'Color','k'); 
  title("Cum incidence for age 60-65");
  xlabel("Time(days)");
  ylabel("Cum incidence");
  ylim = [0,1];
  hold;
plot(epi_base(1).time, (epi_base(1).N_age(agegp)-epi_base(1).S(:,agegp))./epi_base(1).N_age(agegp),'linewidth',1,'Color',[0 0 0]+0.05*10);
plot(epi_march(1).time, (epi_march(1).N_age(agegp)-epi_march(1).S(:,agegp))./epi_march(1).N_age(agegp),'linewidth',1,'Color','b');
plot(epi_april(1).time, (epi_april(1).N_age(agegp)-epi_april(1).S(:,agegp))./epi_april(1).N_age(agegp),'linewidth',1,'Color','r');
  legend("Do Nothing", "Base","Lockdown->March","Lockdown->April");

end