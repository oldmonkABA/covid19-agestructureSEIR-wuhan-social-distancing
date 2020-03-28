clear
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


function Cnorm = normalize_contact_matrices(C,popv,makesym)
% FUN normalize them so that
% the matrix with all the contacts is normalized so that its dominant eigenvalue is 1 
% and other matrices keep their contributions 
    matmul = popv*(1./popv)';
    Csym=struct;
    if makesym % make sure contacts are reciprocal
        Csym.home = (C.home + C.home'.*matmul)/2;
        Csym.work = (C.work + C.work'.*matmul)/2;
        Csym.school = (C.school + C.school'.*matmul)/2;
        Csym.others = (C.others + C.others'.*matmul)/2;
        Csym.all = (C.all + C.all'.*matmul)/2;
    else  % if makesym = F leave it as is
        Csym = C;
    end
    eigen_values = eig(Csym.all);
    eig1 = real(eigen_values(1)); % save dominant eigenvalue of the matrix with all contacts
    Cnorm=struct;
    % divide all matrices by the real part of the dominant matrix of all of the contacts
    % that way all the rescaled matrices still sum to C_all = C_work + C_home + C_school + C_other
    Cnorm.home = Csym.home/eig1;
    Cnorm.work = Csym.work/eig1;
    Cnorm.school = Csym.school/eig1;
    Cnorm.others = Csym.others/eig1;
    Cnorm.all = Csym.all/eig1;
end

