clear;
close all;
clc;
%%  load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
loadPopData = true;
loadContactMatrices = true;
loadCaseData =true;
loadR0posterior =true;

% 1) population data
if loadPopData 
  wuhanpop = csvread('../data/wuhanpop.csv',1);
end

if(loadContactMatrices)

  load '../data/contacts_china.mat' 
  %normalize_contact_matrices(contacts,popv, makesym=true)
end
popv = wuhanpop(:,2);

% case age distribution
if loadCaseData

  wuhancaseraw = csvread('../data/wuhan_pop_case_dist.csv',1);
  caseage = repelem(wuhancaseraw(:,4),2)/2;
  wuhancase = [caseage(1:15); sum(caseage(16:20))];
  
end

if(loadR0posterior)

  % --- read in R0 posterior
  R0_plot =csvread("../data/out_R0.csv",1);
  fid = fopen('../data/out_date.csv');
  data1 = textscan(fid, '%s',  'headerLines', 1);
  R0_dates =datetime(data1{1,1},'InputFormat','yyyy-MM-dd');
  start_date =R0_dates(1); % first case
  end_date = R0_dates(end); % period to forecast ahead
  date_range =start_date:end_date;
  
  % extract all estimates from 01.01.2020 - 23.01.2020
  R0_posterior = R0_plot(find(date_range == datetime(2020,1,1)):find(date_range == datetime(2020,01,23)),:);
  range(R0_posterior)
  r0posterior = reshape(R0_posterior,[],1);
  %figure(1)
  subplot(2,1,1);
  [R0_densey,R0_densex] = ksdensity(r0posterior);
  plot(R0_densex,R0_densey,'k','LineWidth',1)
  xlabel('R0');
  ylabel('Density');
  
  [R0_densey,R0_densex] = ksdensity(log(r0posterior));
  subplot(2,1,2);
  plot(R0_densex,R0_densey,'k','LineWidth',1)
  xlabel('ln(R0)');
  ylabel('Density');
 
end


%% To simulate n_simSEIR outbreaks

nsim = 200;

%set.seed(123)
r0postCrI = r0posterior;
% hist(r0postCrI)
% summary(r0postCrI)
R0est = randsample(r0postCrI,nsim);
% print(R0est)

% To simulate n_simSEIR outbreaks: duration of infection = 3 days, initial infected  n=~200 infected
%start = Sys.time()
%% durInfSim =3
durInfSim = 3;
initialI = 0.0002;


  R0tpostoutbreak = 1.5;
  dateStartSchoolClosure_nothing = datetime(2019,11,01);
  dateStartIntenseIntervention_nothing = datetime(2019,11,01);
  dateEndIntenseIntervention_nothing = datetime(2019,11,01);
  
  dateStartSchoolClosure_default = datetime(2020,01,15); % cause winter term break 
  dateStartIntenseIntervention_default = datetime(2020,01,23); %Intense intervention: starts at Wuhan Lockdown
  dateEndIntenseIntervention_base= datetime(2020,1,31);
  dateEndIntenseIntervention_march=datetime(2020,3,1);
  dateEndIntenseIntervention_april=datetime(2020,4,1);
  dateStart = datetime(2019,11,01);
  pWorkOpen_default = [0.1,0.25,0.5,0.9];
  POP = wuhanpop;
  rho = [ones(1,4)*0.4 ones(1,12)*0.8];%(0.4,4),rep(0.8,12))
  pWorkOpen_nothing = ones(1,4);%(1,1,1,1)
  pWorkOpen_base = [0.1,0.75,1,1];
  numWeekStagger_base = [10/7,10/7,10/7];
  numWeekStagger_nothing = zeros(1,4);%(0,0,0,0)
  numWeekStagger_default=[2,4,6];
  pInfected=initialI;
  durInf = durInfSim;
  contacts_china=contacts;
tic
for sim = 1:nsim
    
  epi_doNothingDurInf3(sim) = simulateOutbreakSEIcIscR(R0est(sim),rho,R0tpostoutbreak,dateEndIntenseIntervention_nothing,...
                                                    pWorkOpen_nothing,dateStartSchoolClosure_nothing,....
                                                    dateStartIntenseIntervention_nothing,...
                                                    dateStart,POP,numWeekStagger_nothing,...
                                                    pInfected,durInf,contacts_china);
   epi_baseDurInf3(sim) = simulateOutbreakSEIcIscR(R0est(sim),rho,R0tpostoutbreak,dateEndIntenseIntervention_base,...
                                               pWorkOpen_base,dateStartSchoolClosure_default,....
                                                    dateStartIntenseIntervention_default,...
                                                    dateStart,...
                                                    POP,numWeekStagger_base,pInfected,durInf,contacts_china);
    epi_marchDurInf3(sim) = simulateOutbreakSEIcIscR(R0est(sim),rho,R0tpostoutbreak,dateEndIntenseIntervention_march,...
                                                    pWorkOpen_default,dateStartSchoolClosure_default,....
                                                    dateStartIntenseIntervention_default,...
                                                    dateStart,...
                                                    POP,numWeekStagger_default,pInfected,durInf,contacts_china);
    epi_aprilDurInf3(sim) = simulateOutbreakSEIcIscR(R0est(sim),rho,R0tpostoutbreak,dateEndIntenseIntervention_april,...
                                                    pWorkOpen_default,dateStartSchoolClosure_default,....
                                                    dateStartIntenseIntervention_default,...
                                                    dateStart,...
                                                    POP,numWeekStagger_default,pInfected,durInf,contacts_china);


  if mod(sim,10)==0 
      display('Done with simulation ',num2str(sim));
  end
  end
toc
%array=CreateAStruct(1);
%covid_SDurInf3 = struct(1); 
covid_SDurInf3(1) = summariseSimulations('S',50,epi_doNothingDurInf3);
covid_SDurInf3(2) = summariseSimulations('S',50,epi_baseDurInf3);
covid_SDurInf3(3) = summariseSimulations('S',50,epi_marchDurInf3);
covid_SDurInf3(4) = summariseSimulations('S',50,epi_aprilDurInf3);

 %covid_IDurInf3 = list() 
covid_IDurInf3(1) = summariseSimulations('incidence',50,epi_doNothingDurInf3)
covid_IDurInf3(2) = summariseSimulations('incidence',50,epi_baseDurInf3)
covid_IDurInf3(3) = summariseSimulations('incidence',50,epi_marchDurInf3)
covid_IDurInf3(4) = summariseSimulations('incidence',50,epi_aprilDurInf3)

%peaktime_DurInf3 = list()
peaktime_DurInf3(1) = summarisePeakTimePeakSize(epi_doNothingDurInf3)
peaktime_DurInf3(2) = summarisePeakTimePeakSize(epi_baseDurInf3)
peaktime_DurInf3(3) = summarisePeakTimePeakSize(epi_marchDurInf3)
peaktime_DurInf3(4) = summarisePeakTimePeakSize(epi_aprilDurInf3)

%covid_DurInf3 = list() 
covid_DurInf3(1) = summariseSimulations_mid(50,epi_doNothingDurInf3)
covid_DurInf3(2) = summariseSimulations_mid(50,epi_baseDurInf3)
covid_DurInf3(3) = summariseSimulations_mid(50,epi_marchDurInf3)
covid_DurInf3(4) = summariseSimulations_mid(50,epi_aprilDurInf3)

%AGEcovid_IDurInf3 = list()
AGEcovid_IDurInf3(1) = summariseSimulationsAGE('incidence',50,epi_doNothingDurInf3)
AGEcovid_IDurInf3(2) = summariseSimulationsAGE('incidence',50,epi_baseDurInf3)
AGEcovid_IDurInf3(3) = summariseSimulationsAGE('incidence',50,epi_marchDurInf3)
AGEcovid_IDurInf3(4) = summariseSimulationsAGE('incidence',50,epi_aprilDurInf3)
% 
epiFirstSimDurInf3=struct;
epiFirstSimDurInf3.epi_doNothingDurInf3 = epi_doNothingDurInf3(1);
epiFirstSimDurInf3.epi_baseDurInf3= epi_baseDurInf3(1);
epiFirstSimDurInf3.epi_marchDurInf3 = epi_marchDurInf3(1);
epiFirstSimDurInf3.epi_aprilDurInf3 = epi_aprilDurInf3(1);

%% durInf =7
  durInf = 7;
  contacts_china=contacts;
tic
for sim = 1:nsim
    
  epi_doNothingDurInf7(sim) = simulateOutbreakSEIcIscR(R0est(sim),rho,R0tpostoutbreak,dateEndIntenseIntervention_nothing,...
                                                    pWorkOpen_nothing,dateStartSchoolClosure_nothing,....
                                                    dateStartIntenseIntervention_nothing,...
                                                    dateStart,POP,numWeekStagger_nothing,...
                                                    pInfected,durInf,contacts_china);
   epi_baseDurInf7(sim) = simulateOutbreakSEIcIscR(R0est(sim),rho,R0tpostoutbreak,dateEndIntenseIntervention_base,...
                                               pWorkOpen_base,dateStartSchoolClosure_default,....
                                                    dateStartIntenseIntervention_default,...
                                                    dateStart,...
                                                    POP,numWeekStagger_base,pInfected,durInf,contacts_china);
    epi_marchDurInf7(sim) = simulateOutbreakSEIcIscR(R0est(sim),rho,R0tpostoutbreak,dateEndIntenseIntervention_march,...
                                                    pWorkOpen_default,dateStartSchoolClosure_default,....
                                                    dateStartIntenseIntervention_default,...
                                                    dateStart,...
                                                    POP,numWeekStagger_default,pInfected,durInf,contacts_china);
    epi_aprilDurInf7(sim) = simulateOutbreakSEIcIscR(R0est(sim),rho,R0tpostoutbreak,dateEndIntenseIntervention_april,...
                                                    pWorkOpen_default,dateStartSchoolClosure_default,....
                                                    dateStartIntenseIntervention_default,...
                                                    dateStart,...
                                                    POP,numWeekStagger_default,pInfected,durInf,contacts_china);


  if mod(sim,10)==0 
      display('Done with simulation ',num2str(sim));
  end
  end
toc

% covid_SDurInf3 = struct(1); 
covid_SDurInf7(1) = summariseSimulations('S',50,epi_doNothingDurInf7);
covid_SDurInf7(2) = summariseSimulations('S',50,epi_baseDurInf7);
covid_SDurInf7(3) = summariseSimulations('S',50,epi_marchDurInf7);
covid_SDurInf7(4) = summariseSimulations('S',50,epi_aprilDurInf7);

% covid_IDurInf7 = list() 
covid_IDurInf7(1) = summariseSimulations('incidence',50,epi_doNothingDurInf7);
covid_IDurInf7(2) = summariseSimulations('incidence',50,epi_baseDurInf7);
covid_IDurInf7(3) = summariseSimulations('incidence',50,epi_marchDurInf7);
covid_IDurInf7(4) = summariseSimulations('incidence',50,epi_aprilDurInf7);

%peaktime_DurInf7 = list()
peaktime_DurInf7(1) = summarisePeakTimePeakSize(epi_doNothingDurInf7);
peaktime_DurInf7(2) = summarisePeakTimePeakSize(epi_baseDurInf7);
peaktime_DurInf7(3) = summarisePeakTimePeakSize(epi_marchDurInf7)
peaktime_DurInf7(4) = summarisePeakTimePeakSize(epi_aprilDurInf7);;

%covid_DurInf7 = list() 
covid_DurInf7(1) = summariseSimulations_mid(50,epi_doNothingDurInf7);
covid_DurInf7(2) = summariseSimulations_mid(50,epi_baseDurInf7);
covid_DurInf7(3) = summariseSimulations_mid(50,epi_marchDurInf7);
covid_DurInf7(4) = summariseSimulations_mid(50,epi_aprilDurInf7);

%AGEcovid_IDurInf7 = list()
AGEcovid_IDurInf7(1) = summariseSimulationsAGE('incidence',50,epi_doNothingDurInf7);
AGEcovid_IDurInf7(2) = summariseSimulationsAGE('incidence',50,epi_baseDurInf7);
AGEcovid_IDurInf7(3) = summariseSimulationsAGE('incidence',50,epi_marchDurInf7);
AGEcovid_IDurInf7(4) = summariseSimulationsAGE('incidence',50,epi_aprilDurInf7);

epiFirstSimDurInf7=struct;
epiFirstSimDurInf7.epi_doNothingDurInf7 = epi_doNothingDurInf7(1);
epiFirstSimDurInf7.epi_baseDurInf7= epi_baseDurInf7(1);
epiFirstSimDurInf7.epi_marchDurInf7 = epi_marchDurInf7(1);
epiFirstSimDurInf7.epi_aprilDurInf7 = epi_aprilDurInf7(1);

% save('outputs/SEIR/covid_IDurInf3.mat','-struct','covid_IDurInf3');
% save('outputs/SEIR/covid_SDurInf3.mat','-struct','covid_SDurInf3');
% save('outputs/SEIR/covid_IDurInf7.mat','-struct','covid_IDurInf7');
% save('outputs/SEIR/covid_SDurInf7.mat','-struct','covid_SDurInf7');
% 
% save('outputs/SEIR/peaktime_DurInf3.mat','-struct','peaktime_DurInf3');
% save('outputs/SEIR/peaktime_DurInf7.mat','-struct','peaktime_DurInf7');
% 
% save('outputs/SEIR/covid_DurInf3.mat','-struct','covid_DurInf3');
% save('outputs/SEIR/covid_DurInf7.mat','-struct','covid_DurInf7');
% 
% save('outputs/SEIR/AGEcovid_IDurInf3.mat','-struct','AGEcovid_IDurInf3');
% save('outputs/SEIR/AGEcovid_IDurInf7.mat','-struct','AGEcovid_IDurInf7');
% 
% save('outputs/SEIR/epiFirstSimDurInf3.mat','-struct','epiFirstSimDurInf3');
% save('outputs/SEIR/epiFirstSimDurInf7.mat','-struct','epiFirstSimDurInf7');

%% durInfSim = 3 , initialI = 0.002

durInfSim = 3;
initialI = 0.002;

R0tpostoutbreak = 1.5;
  dateStartSchoolClosure_nothing = datetime(2019,11,01);
  dateStartIntenseIntervention_nothing = datetime(2019,11,01);
  dateEndIntenseIntervention_nothing = datetime(2019,11,01);
  
  dateStartSchoolClosure_default = datetime(2020,01,15); % cause winter term break 
  dateStartIntenseIntervention_default = datetime(2020,01,23); %Intense intervention: starts at Wuhan Lockdown
  dateEndIntenseIntervention_base= datetime(2020,1,31);
  dateEndIntenseIntervention_march=datetime(2020,3,1);
  dateEndIntenseIntervention_april=datetime(2020,4,1);
  dateStart = datetime(2019,11,01);
  pWorkOpen_default = [0.1,0.25,0.5,0.9];
  POP = wuhanpop;
  %rho = ones(1,3660)*0.5;
  pWorkOpen_nothing = ones(1,4);%(1,1,1,1)
  pWorkOpen_base = [0.1,0.75,1,1];
  numWeekStagger_base = [10/7,10/7,10/7];
  numWeekStagger_nothing = zeros(1,4);%(0,0,0,0)
  numWeekStagger_default=[2,4,6];
  pInfected=initialI;
  durInf = durInfSim;
  contacts_china=contacts;
tic
for sim = 1:nsim
    
  epi_doNothingDurInf3I2000(sim) = simulateOutbreakSEIcIscR(R0est(sim),rho,R0tpostoutbreak,dateEndIntenseIntervention_nothing,...
                                                    pWorkOpen_nothing,dateStartSchoolClosure_nothing,....
                                                    dateStartIntenseIntervention_nothing,...
                                                    dateStart,POP,numWeekStagger_nothing,...
                                                    pInfected,durInf,contacts_china);
   epi_baseDurInf3I2000(sim) = simulateOutbreakSEIcIscR(R0est(sim),rho,R0tpostoutbreak,dateEndIntenseIntervention_base,...
                                               pWorkOpen_base,dateStartSchoolClosure_default,....
                                                    dateStartIntenseIntervention_default,...
                                                    dateStart,...
                                                    POP,numWeekStagger_base,pInfected,durInf,contacts_china);
    epi_marchDurInf3I2000(sim) = simulateOutbreakSEIcIscR(R0est(sim),rho,R0tpostoutbreak,dateEndIntenseIntervention_march,...
                                                    pWorkOpen_default,dateStartSchoolClosure_default,....
                                                    dateStartIntenseIntervention_default,...
                                                    dateStart,...
                                                    POP,numWeekStagger_default,pInfected,durInf,contacts_china);
    epi_aprilDurInf3I2000(sim) = simulateOutbreakSEIcIscR(R0est(sim),rho,R0tpostoutbreak,dateEndIntenseIntervention_april,...
                                                    pWorkOpen_default,dateStartSchoolClosure_default,....
                                                    dateStartIntenseIntervention_default,...
                                                    dateStart,...
                                                    POP,numWeekStagger_default,pInfected,durInf,contacts_china);


  if mod(sim,10)==0 
      display('Done with simulation ',num2str(sim));
  end
  end
toc
%array=CreateAStruct(1);
%covid_SDurInf3I2000 = struct(1); 
covid_SDurInf3I2000(1) = summariseSimulations('S',50,epi_doNothingDurInf3I2000);
covid_SDurInf3I2000(2) = summariseSimulations('S',50,epi_baseDurInf3I2000);
covid_SDurInf3I2000(3) = summariseSimulations('S',50,epi_marchDurInf3I2000);
covid_SDurInf3I2000(4) = summariseSimulations('S',50,epi_aprilDurInf3I2000);

 %covid_IDurInf3I2000 = list() 
covid_IDurInf3I2000(1) = summariseSimulations('incidence',50,epi_doNothingDurInf3I2000)
covid_IDurInf3I2000(2) = summariseSimulations('incidence',50,epi_baseDurInf3I2000)
covid_IDurInf3I2000(3) = summariseSimulations('incidence',50,epi_marchDurInf3I2000)
covid_IDurInf3I2000(4) = summariseSimulations('incidence',50,epi_aprilDurInf3I2000)

%peaktime_DurInf3I2000 = list()
peaktime_DurInf3I2000(1) = summarisePeakTimePeakSize(epi_doNothingDurInf3I2000)
peaktime_DurInf3I2000(2) = summarisePeakTimePeakSize(epi_baseDurInf3I2000)
peaktime_DurInf3I2000(3) = summarisePeakTimePeakSize(epi_marchDurInf3I2000)
peaktime_DurInf3I2000(4) = summarisePeakTimePeakSize(epi_aprilDurInf3I2000)

%covid_DurInf3I2000 = list() 
covid_DurInf3I2000(1) = summariseSimulations_mid(50,epi_doNothingDurInf3I2000)
covid_DurInf3I2000(2) = summariseSimulations_mid(50,epi_baseDurInf3I2000)
covid_DurInf3I2000(3) = summariseSimulations_mid(50,epi_marchDurInf3I2000)
covid_DurInf3I2000(4) = summariseSimulations_mid(50,epi_aprilDurInf3I2000)

%AGEcovid_IDurInf3I2000 = list()
AGEcovid_IDurInf3I2000(1) = summariseSimulationsAGE('incidence',50,epi_doNothingDurInf3I2000)
AGEcovid_IDurInf3I2000(2) = summariseSimulationsAGE('incidence',50,epi_baseDurInf3I2000)
AGEcovid_IDurInf3I2000(3) = summariseSimulationsAGE('incidence',50,epi_marchDurInf3I2000)
AGEcovid_IDurInf3I2000(4) = summariseSimulationsAGE('incidence',50,epi_aprilDurInf3I2000)
% 
epiFirstSimDurInf3I2000=struct;
epiFirstSimDurInf3I2000.epi_doNothingDurInf3I2000 = epi_doNothingDurInf3I2000(1);
epiFirstSimDurInf3I2000.epi_baseDurInf3I2000= epi_baseDurInf3I2000(1);
epiFirstSimDurInf3I2000.epi_marchDurInf3I2000 = epi_marchDurInf3I2000(1);
epiFirstSimDurInf3I2000.epi_aprilDurInf3I2000 = epi_aprilDurInf3I2000(1);

%% durInfSim = 7 initialI = 0.002

durInfSim = 7;
initialI = 0.002;

R0tpostoutbreak = 1.5;
  dateStartSchoolClosure_nothing = datetime(2019,11,01);
  dateStartIntenseIntervention_nothing = datetime(2019,11,01);
  dateEndIntenseIntervention_nothing = datetime(2019,11,01);
  
  dateStartSchoolClosure_default = datetime(2020,01,15); % cause winter term break 
  dateStartIntenseIntervention_default = datetime(2020,01,23); %Intense intervention: starts at Wuhan Lockdown
  dateEndIntenseIntervention_base= datetime(2020,1,31);
  dateEndIntenseIntervention_march=datetime(2020,3,1);
  dateEndIntenseIntervention_april=datetime(2020,4,1);
  dateStart = datetime(2019,11,01);
  pWorkOpen_default = [0.1,0.25,0.5,0.9];
  POP = wuhanpop;
  %rho = ones(1,3660)*0.5;
  pWorkOpen_nothing = ones(1,4);%(1,1,1,1)
  pWorkOpen_base = [0.1,0.75,1,1];
  numWeekStagger_base = [10/7,10/7,10/7];
  numWeekStagger_nothing = zeros(1,4);%(0,0,0,0)
  numWeekStagger_default=[2,4,6];
  pInfected=initialI;
  durInf = durInfSim;
  contacts_china=contacts;
tic
for sim = 1:nsim
    
  epi_doNothingDurInf7I2000(sim) = simulateOutbreakSEIcIscR(R0est(sim),rho,R0tpostoutbreak,dateEndIntenseIntervention_nothing,...
                                                    pWorkOpen_nothing,dateStartSchoolClosure_nothing,....
                                                    dateStartIntenseIntervention_nothing,...
                                                    dateStart,POP,numWeekStagger_nothing,...
                                                    pInfected,durInf,contacts_china);
   epi_baseDurInf7I2000(sim) = simulateOutbreakSEIcIscR(R0est(sim),rho,R0tpostoutbreak,dateEndIntenseIntervention_base,...
                                               pWorkOpen_base,dateStartSchoolClosure_default,....
                                                    dateStartIntenseIntervention_default,...
                                                    dateStart,...
                                                    POP,numWeekStagger_base,pInfected,durInf,contacts_china);
    epi_marchDurInf7I2000(sim) = simulateOutbreakSEIcIscR(R0est(sim),rho,R0tpostoutbreak,dateEndIntenseIntervention_march,...
                                                    pWorkOpen_default,dateStartSchoolClosure_default,....
                                                    dateStartIntenseIntervention_default,...
                                                    dateStart,...
                                                    POP,numWeekStagger_default,pInfected,durInf,contacts_china);
    epi_aprilDurInf7I2000(sim) = simulateOutbreakSEIcIscR(R0est(sim),rho,R0tpostoutbreak,dateEndIntenseIntervention_april,...
                                                    pWorkOpen_default,dateStartSchoolClosure_default,....
                                                    dateStartIntenseIntervention_default,...
                                                    dateStart,...
                                                    POP,numWeekStagger_default,pInfected,durInf,contacts_china);


  if mod(sim,10)==0 
      display('Done with simulation ',num2str(sim));
  end
  end
toc
%array=CreateAStruct(1);
%covid_SDurInf7I2000 = struct(1); 
covid_SDurInf7I2000(1) = summariseSimulations('S',50,epi_doNothingDurInf7I2000);
covid_SDurInf7I2000(2) = summariseSimulations('S',50,epi_baseDurInf7I2000);
covid_SDurInf7I2000(3) = summariseSimulations('S',50,epi_marchDurInf7I2000);
covid_SDurInf7I2000(4) = summariseSimulations('S',50,epi_aprilDurInf7I2000);

 %covid_IDurInf7I2000 = list() 
covid_IDurInf7I2000(1) = summariseSimulations('incidence',50,epi_doNothingDurInf7I2000)
covid_IDurInf7I2000(2) = summariseSimulations('incidence',50,epi_baseDurInf7I2000)
covid_IDurInf7I2000(3) = summariseSimulations('incidence',50,epi_marchDurInf7I2000)
covid_IDurInf7I2000(4) = summariseSimulations('incidence',50,epi_aprilDurInf7I2000)

%peaktime_DurInf7I2000 = list()
peaktime_DurInf7I2000(1) = summarisePeakTimePeakSize(epi_doNothingDurInf7I2000)
peaktime_DurInf7I2000(2) = summarisePeakTimePeakSize(epi_baseDurInf7I2000)
peaktime_DurInf7I2000(3) = summarisePeakTimePeakSize(epi_marchDurInf7I2000)
peaktime_DurInf7I2000(4) = summarisePeakTimePeakSize(epi_aprilDurInf7I2000)

%covid_DurInf7I2000 = list() 
covid_DurInf7I2000(1) = summariseSimulations_mid(50,epi_doNothingDurInf7I2000)
covid_DurInf7I2000(2) = summariseSimulations_mid(50,epi_baseDurInf7I2000)
covid_DurInf7I2000(3) = summariseSimulations_mid(50,epi_marchDurInf7I2000)
covid_DurInf7I2000(4) = summariseSimulations_mid(50,epi_aprilDurInf7I2000)

%AGEcovid_IDurInf7I2000 = list()
AGEcovid_IDurInf7I2000(1) = summariseSimulationsAGE('incidence',50,epi_doNothingDurInf7I2000)
AGEcovid_IDurInf7I2000(2) = summariseSimulationsAGE('incidence',50,epi_baseDurInf7I2000)
AGEcovid_IDurInf7I2000(3) = summariseSimulationsAGE('incidence',50,epi_marchDurInf7I2000)
AGEcovid_IDurInf7I2000(4) = summariseSimulationsAGE('incidence',50,epi_aprilDurInf7I2000)
% 
epiFirstSimDurInf7I2000=struct;
epiFirstSimDurInf7I2000.epi_doNothingDurInf7I2000 = epi_doNothingDurInf7I2000(1);
epiFirstSimDurInf7I2000.epi_baseDurInf7I2000= epi_baseDurInf7I2000(1);
epiFirstSimDurInf7I2000.epi_marchDurInf7I2000 = epi_marchDurInf7I2000(1);
epiFirstSimDurInf7I2000.epi_aprilDurInf7I2000 = epi_aprilDurInf7I2000(1);
