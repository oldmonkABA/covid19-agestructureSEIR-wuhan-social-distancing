function results = summariseSimulationsAGE(VAR,CI,SIMS)

  [var_p_median , var_p_lci , var_p_uci] = deal(zeros(length(SIMS(1).time),16));
  for age = 1:16
    temp=[];
    temp1=[];
    for i = 1: length(SIMS)
      temp=[temp SIMS(i).incidence];%row sums
    end  
    temp1 = temp;
    %temp = lapply(SIMS,FUN = function(x) (x[[VAR]][,age]))
    %temp1 = do.call(cbind.data.frame, temp)
    var_p_median(:,age) = quantile(temp1,0.5,2);
    var_p_lci(:,age) = quantile(temp1,(1-CI/100)/2,2);
    var_p_uci(:,age) = quantile(temp1,1-(1-CI/100)/2,2);
             
  SUMMARY=struct;
  SUMMARY.median=var_p_median;
  SUMMARY.lci = var_p_lci;
  SUMMARY.uci = var_p_uci;
  results = struct;
  results.summary = SUMMARY;
  results.Sim1 = SIMS(1);
  end