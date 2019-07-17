function eqindex = extract_eqindex(t,V,Veq)
% function to extract earthquake indices
% INPUT - 
% t - time array (seconds)
% V - space-time array of velocities for all patches (m/s)
% Veq - velocity threshold for detecting an earthquake (generally set as 1e-3 m/s
% OUTPUT - 
% eqindex - indices of detected earthquakes
% Rishav Mallick, EOS, 2018

% set criteria for detecting an earthquake
% Veq = 1e-3;
j = 1;
eqindex =[];

for i = 1:length(t)-1
    if i > 1 && max(V(i,:)) >= Veq && max(V(i-1,:)) < Veq 
        eqindex(j)=i;
        j = j+1;
    end
end

writetable(table([1:length(eqindex)]',eqindex'),'results/eqindex.dat','Delimiter','\t','WriteVariableNames',0)
end