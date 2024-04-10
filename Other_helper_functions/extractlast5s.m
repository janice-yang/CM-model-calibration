function [extract] = extractlast5s(k19cm)
tlast = k19cm.state.t(end);
tstart = tlast - 5500;
tstartdex = find(k19cm.state.t >= tstart);
tstartdex = tstartdex(1);

tend = tlast - 499;
tenddex = find(k19cm.state.t >= tend);
tenddex = tenddex(1);

t = k19cm.state.t(tstartdex:tenddex)-k19cm.state.t(tstartdex);
V = k19cm.state.Y(tstartdex:tenddex,1);
Cai = k19cm.state.Y(tstartdex:tenddex,3);
init = k19cm.state.Y(tstartdex,:)';
extract = {init, t, V, Cai};
end