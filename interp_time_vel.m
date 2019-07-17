Veq = 1e-3;
eqnum = 3;
eqindex = extract_eqindex(t,V,Veq);
eqindex(end+1) = t(end);


s2d = 60*60*24;
td = t(eqindex(eqnum)):s2d:t(eqindex(eqnum+1));
for i = 1:ss.M
    Vint(:,i) = interp1(t,V