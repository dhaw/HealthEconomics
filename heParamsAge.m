function [piligs,phgili,pdgh]=heParamsAge%(tab)
%Then hard code outputs into prep file!
%%
%Work in row vectors
lx=64;
%England and Wales, ONS 2019 midpoint estimate:
nn=[3463,3726,3538,3260,3693,4022,4011,3926,3586,3919,4129,3890,3308,2982,2960,2069,1531+933+414+117+12];
nntot=[nn(1),sum(nn(2:4)),sum(nn(5:13)),sum(nn(13:end))];
ranges=[1,3,9,4];%Last group to 80+
nntot=repelem(nntot,ranges);
nnprop=nn./nntot;%5-year group in larger group
subs=1:4;
subs=repelem(subs,ranges);
%%
%ILI given infection:
piligs=[0.467,0.467,0.485,0.485,0.550,0.550,0.601,0.601,0.639,0.639,0.757,0.757,0.860,0.860,0.958,0.958,1.000];
piligs=accumarray(subs',piligs'.*nnprop');
%Hosp given ILI:
phgili=[0.009,0.007,0.013,0.027,0.053,0.101,0.147,0.211,0.246,0.328,0.414,0.519,0.641,0.738,0.827,1.000,0.976];
phgili=accumarray(subs',phgili'.*nnprop');
%Death given hosp:
pdgh=[0.022,0.022,0.022,0.022,0.022,0.022,0.022,0.023,0.026,0.033,0.047,0.072,0.118,0.181,0.257,0.351,1.000];
pdgh=accumarray(subs',pdgh'.*nnprop');
%%
tinfh=5;
thdeath=10.81;%Median 8
threc=12.73;%Median 9
pdeath=.39;

%p2=pili;%/tinfh;
%h=pdeath/thdeath;%mu
%gX=(1-pdeath)/threc
%g3=(1-pdeath)/threc;%gamma3
end
