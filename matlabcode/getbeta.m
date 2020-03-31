function results = getbeta(R0t,constraints,gamma,p_age,calculate_transmission_probability,CONTACTMATRIX)
% 1) R0
% 2) gamma = removal rate  
% 3) f = population age proportion 
% 4) constraints = a scale matrix contstraint age- and location-specific contact matrices (a linear combination over all locations; TODO to specify carefully based on interventions)
% 5) calculate_transmission_probability if this is 1, then calculate the transmission probability from R0 otherwise, assume it is beta=0.05 
% 6) npop = population size 
  
% constraints for age-specific contacts at home, work, school, others
  if nargin < 5
      calculate_transmission_probability=1;
      CONTACTMATRIX = contacts;
  end
  n = 16; %length(p_age)
  diagonal = ones(1,16);
  constraints_base=struct;
  constraints_base.home = diag(diagonal);
  constraints_base.work = diag(diagonal);
  constraints_base.school = diag(diagonal);
  constraints_base.others = diag(diagonal);
  
  matmul = p_age*(1./p_age)';
  Csym=struct;
  Csym.home = (CONTACTMATRIX.home + CONTACTMATRIX.home'.*matmul)/2
  Csym.work = (CONTACTMATRIX.work + CONTACTMATRIX.work'.*matmul)/2
  Csym.school = (CONTACTMATRIX.school + CONTACTMATRIX.school'.*matmul)/2
  Csym.others = (CONTACTMATRIX.others + CONTACTMATRIX.others'.*matmul)/2
  Csym.all = (CONTACTMATRIX.all + CONTACTMATRIX.all'.*matmul)/2
  
  CONTACTMATRIX=Csym;
  
  C = constraints_base.home*CONTACTMATRIX.home + constraints_base.work*CONTACTMATRIX.work+ constraints_base.school*CONTACTMATRIX.school+ constraints_base.others*CONTACTMATRIX.others;
  
  if calculate_transmission_probability==1
      M = C;
      for i= 1:n
          for j = 1:n
              M(i,j) = C(i,j)*p_age(i)/p_age(j)
          end
      end
      eigen = eig(M);
      beta = R0t*gamma/max(real(eigen));  % reverse engineer beta from the R0 and gamma
      
  else
      beta = 0.025;%0.05
  end
  results=beta;
end