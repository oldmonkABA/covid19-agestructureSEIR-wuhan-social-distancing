function results=summariseSimulations(VAR,CI,SIMS)
  temp1 =[]; 
  switch VAR
      case "S"
          for i = 1: length(SIMS)
              epi = sum(SIMS(i).S,2);
              temp1=[temp1 epi];
          end
      case 'incidence'
          for i = 1: length(SIMS)
              epi = sum(SIMS(i).incidence,2);
              temp1=[temp1 epi];
          end
  end

  var_p_median = quantile(temp1,0.5,2);
  var_p_lci = quantile(temp1,(1-CI/100)/2,2);
  var_p_uci = quantile(temp1,1-(1-CI/100)/2,2);
  SUMMARY=struct;
  SUMMARY.median = var_p_median;
  SUMMARY.lci = var_p_lci;
  SUMMARY.uci=var_p_uci;
  results = struct; 
  results.summary=SUMMARY;
  results.Sim1 = SIMS(1);
  
  
  end