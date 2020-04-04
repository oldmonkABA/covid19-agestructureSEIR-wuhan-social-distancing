function results =summariseSimulations_mid(CI,SIMS)

  temp1 =[]; 
  for i = 1: length(SIMS)
      epi = sum(SIMS(i).S,2);
      temp1=[temp1 epi];
  end
  lastrow = temp1(end,:);
  [m,i] = min(abs(lastrow-quantile(lastrow,0.5)));
  [m,j] = min(abs(lastrow-quantile(lastrow,0.25)));
  [m,k] = min(abs(lastrow-quantile(lastrow,0.75)));
% size(temp1(:,i))
%size(SIMS(1).time)
  S = [temp1(:,i),temp1(:,j),temp1(:,k),SIMS(1).time'];
  S_age=struct;
  S_age.med=SIMS(i).S;
  S_age.lci=SIMS(k).S;
  S_age.uci=SIMS(j).S;
  inc = struct;
  inc.med=SIMS(i).incidence;
  inc.lci=SIMS(k).incidence;
  inc.uci=SIMS(j).incidence;
  time = SIMS(1).time;
  N_age=SIMS(1).N_age;
  results = struct;
  results.S=S;
  results.inc=inc;
  results.time=time;
  results.N_age=N_age;
  results.S_age=S_age;
end
