%% Frequency of Firing
% Generate U first from T and Conta

count=sum(U{1}.R_ntk(:)); %summation of all spikes for cell
time=length(U{1}.R_ntk(:)); %all time in ms 
freq=(count/time)*1000 %spks/s

%% L3 
0.9263 - 87C
0.7634 - 86B
4.7944 - 61B
0.4808 - 61A
1.5738 - 82B
2.8084 - 76B
10.2248 - 75A
0.2258 - 62A
53.6573 - 60E
0.3472 - 60D
0.7025 - 60C
0.3250 - 52B
0.4665 - 35A
0.5405 - 36B
0.8969 - 37E
0.1964 - 37C
1.1298 - 37A
AVERAGE = 4.98
%% L5B
4.1559 - 85B
12.4439 - 77C
19.7050 - 81F
23.7888 - 81D
2.1168 - 72E
4.8012 - 52E
38.4590 - 52D
5.4130 - 51D
4.8830 - 51C
16.2775 - 50D
5.1850 - 38B
AVERAGE = 12.48