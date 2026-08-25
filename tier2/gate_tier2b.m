%% gate_tier2b.m -- THE go/no-go test for tier 2.
% Take grid points that tier 1 ACCEPTED while sitting on the outer wall -- the
% actual escaping ridge, not a guessed direction -- and ask whether tier 2
% rejects them. If it does, tier 2 closes the region. If not, it does not.
clear; clc;
G=load('grid4c.mat'); theta0=G.theta0; task=G.task; N=200; tS=[16 23 30 37 44 51];
sz=[numel(G.su_range) numel(G.ss_range) numel(G.d_range) numel(G.r_range)];
Sig1=statsCovariance(theta0,N,400,60000,task); i3=[1 3 4]; S1=Sig1(i3,i3);
thr1=hotelling(3,400);
s_obs1=ofcSummary(theta0,N,70001,task);
P3=reshape(G.Spred(i3,:,:,:,:),3,[]);
e=P3-s_obs1(i3); D1=reshape(sum(e.*(S1\e),1),sz);
[A,B,C,E]=ndgrid(1:sz(1),1:sz(2),1:sz(3),1:sz(4));
wall = A==1|B==1|E==1;                       % the LOW walls that escape
cand = find(D1<=thr1 & wall);
[~,o]=sort(D1(cand)); cand=cand(o);
pick = cand(round(linspace(1,numel(cand),min(5,numel(cand)))));
fprintf('subject 1: %d accepted points, %d of them on a low wall\n\n', ...
        sum(D1(:)<=thr1), numel(cand));
Sig2=statsCovariance2(theta0,N,400,60000,task,tS); thr2=hotelling(18,400);
s_obs2=ofcSummary2(theta0,N,70001,task,tS);
fprintf('%-30s %10s %10s %10s %10s\n','theta on the escaping ridge', ...
        'D_tier1','/thr1','D_tier2','/thr2');
for q=1:numel(pick)
    [a,b,c,ee]=ind2sub(sz,pick(q));
    th=[G.su_range(a); G.ss_range(b); G.d_range(c); G.r_range(ee)];
    s2=ofcSummary2(th,20000,12345,task,tS);
    D2=discrepancy(s_obs2,s2,Sig2);
    fprintf('su%.3f ss%.3f d%d r%.2e %10.2f %10.2f %10.1f %10.2f\n', ...
            th(1),th(2),th(3),th(4), D1(pick(q)), D1(pick(q))/thr1, D2, D2/thr2);
end
fprintf('\nthr1 %.2f  thr2 %.2f\n', thr1, thr2);
fprintf(['PASS if D_tier2/thr2 > 1 where D_tier1/thr1 < 1: tier 2 rejects what\n' ...
         'tier 1 accepted, so the ridge is cut and the region closes.\n']);
