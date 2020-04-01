function results=summarisePeakTimePeakSize(SIMS)
temp1=[];
  time = SIMS(1).time;
  for i = 1: length(SIMS)
      epi = sum(SIMS(i).incidence,2);%row sums
      temp1=[temp1 epi];
  end
  [peaksize, ind] = max(temp1,[],1);
  peaksize = time(ind);
  results = struct;
  results.time=time;
  results.peaktime=peaktime;
  results.peaksize=peaksize;
  end