%
% Script that contains Shaw's 1985 data and plots a horizontal diffuse-field
% response. Really should take sqrt of average squared magnitude...
%
% Data is from:
%
% E. A. G. Shaw and M. Vaillancourt, "Transformation of sound-pressure level
% from the free field to the eardrum presented in numerical from,"
% J. Acoust. Soc. Am., vol 78, pg 1120-1123, 1985.
%

T1 = [ 
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 ; 
0.4	0.4	0.5	0.5	0.8	1.1	1.4	1.4	1.4	1.3	1.4	1.4	1.5	1.5	1.4 ; 
0.7	0.8	1.0	1.1	1.7	2.2	2.6	2.6	2.6	2.5	2.6	2.7	2.9	2.8	2.7 ; 
1.0	1.2	1.4	1.5	2.4	3.2	3.6	3.6	3.5	3.5	3.7	3.8	4.0	3.9	3.7 ; 
1.3	1.5	1.8	1.9	2.8	3.8	4.3	4.3	4.3	4.4	4.5	4.7	4.9	4.9	4.7 ; 
1.5	1.7	2.1	2.2	3.2	4.3	4.8	4.8	4.8	4.9	5.1	5.3	5.7	5.7	5.5 ; 
1.6	1.8	2.2	2.3	3.4	4.5	5.0	5.0	5.0	5.1	5.3	5.6	6.1	6.1	5.7 ; 
1.6	1.8	2.1	2.2	3.3	4.3	4.8	4.9	4.9	5.0	5.2	5.5	6.1	6.0	5.5 ; 
1.4	1.6	1.9	2.0	3.0	3.9	4.4	4.5	4.6	4.7	5.0	5.4	6.0	5.9	5.4 ; 
1.1	1.3	1.5	1.6	2.4	3.1	3.5	3.6	3.8	4.1	4.5	5.0	5.9	5.9	5.4 ; 
0.7	0.8	0.9	0.9	1.4	2.0	2.4	2.5	2.7	3.2	3.7	4.3	5.2	5.2	5.0 ; 
0.2	0.1	0	0	0.3	0.7	1.1	1.2	1.4	1.9	2.6	3.3	4.1	4.1	3.9 ; 
-.2	-.5	-.7	.7	-.7	-.6	-.4	-.3	-.1	0.4	1.1	1.8	2.3	2.3	2.1 ; 
-.5	-.9	-1.2	-1.3	-1.6	-1.7	-1.6	-1.6	-1.4	-1.2	-.8	-.3	0.2	0.2	0.2 ; 
-.8	-1.3	-1.7	-1.8	-2.2	-2.5	-2.6	-2.6	-2.7	-2.8	-2.9	-2.8	-2.3	-2.2	-2.0 ; 
-1.1	-1.5	-1.8	-1.9	-2.2	-2.5	-3.0	-3.1	-3.5	-3.9	-4.3	-4.5	-4.5	-4.5	-4.5 ; 
-1.4	-1.6	-1.7	-1.7	-1.9	-2.0	-2.5	-2.6	-2.9	-3.2	-3.4	-3.6	-3.9	-4.1	-5.2 ; 
-1.6	-1.6	-1.5	-1.5	-1.4	-1.4	-1.7	-1.8	-1.9	-2.0	-2.0	-2.0	-2.1	-2.4	-3.8 ; 
-1.6	-1.5	-1.4	-1.4	-1.2	-1.0	-1.2	-1.3	-1.5	-1.5	-1.4	-1.2	-1.0	-1.1	-1.7 ; 
-1.5	-1.4	-1.4	-1.4	-1.4	-1.4	-1.6	-1.7	-2.0	-2.3	-2.4	-2.5	-2.9	-3.1	-4.0 ; 
-1.3	-1.3	-1.4	-1.4	-1.8	-2.2	-2.5	-2.6	-2.9	-3.4	-3.9	-4.4	-5.9	-6.3	-7.1 ; 
-1.0	-1.1	-1.3	-1.4	-1.9	-2.4	-2.9	-3.0	-3.4	-4.1	-4.8	-5.4	-6.3	-6.3	-6.1 ; 
-0.7	-0.8	-1.0	-1.1	-1.5	-2.0	-2.4	-2.5	-2.8	-3.2	-3.6	-3.9	-4.0	-3.9	-3.6 ; 
-0.4	-0.4	-0.5	-0.5	-0.8	-1.1	-1.3	-1.3	-1.4	-1.6	-1.8	-1.9	-1.8	-1.8	-1.6 ; 
];

T2 = [
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 ; 
1.3	1.2	1.2	1.3	1.5	1.8	2.0	2.1	2.2	2.2	2.1	2.1	2.4	2.6	2.6 ;
2.5	2.3	2.3	2.4	2.8	3.1	3.5	3.7	3.9	3.7	3.5	3.8	4.5	5.0	5.1 ;
3.4	3.1	3.0	3.2	3.6	4.1	4.6	4.8	4.8	4.3	4.0	4.6	5.8	6.8	7.4 ;
4.2	3.7	3.3	3.5	4.0	4.6	4.9	4.8	4.6	3.8	3.3	4.4	6.2	7.7	8.8 ;
4.6	3.6	2.9	3.0	3.5	4.0	4.1	4.0	3.5	2.5	1.8	3.0	5.6	7.9	9.4 ;
4.5	3.1	2.0	1.7	2.1	2.4	2.4	2.3	1.9	0.7	-0.3	0.8	4.0	7.0	9.3 ;
4.1	2.2	0.8	-0.2	-0.2	0	0.1	0.1	-0.1	-1.2	-2.4	-1.7	1.7	4.8	7.7 ;
3.9	2.0	0.2	-1.3	-1.5	-1.4	-1.2	-1.2	-1.4	-2.3	-4.3	-4.2	-1.0	2.0	4.8 ;
4.0	2.2	0.4	-1.6	-1.9	-1.8	-1.8	-1.8	-1.9	-2.7	-5.2	-6.0	-3.8	-0.8	1.7 ;
3.8	2.2	0.6	-0.9	-1.5	-1.4	-1.3	-1.3	-1.4	-2.5	-4.9	-5.9	-4.9	-3.0	-1.0 ;
2.8	1.2	-0.2	-1.2	-1.5	-1.6	-1.5	-1.4	-1.5	-2.4	-4.3	-5.3	-4.9	-4.2	-3.3 ;
1.3	-0.2	-1.6	-2.5	-2.6	-2.5	-2.5	-2.5	-2.6	-3.4	-5.0	-5.4	-4.5	-5.0	-5.1 ;
-0.4	-2.0	-3.4	-4.3	-4.4	-4.3	-4.2	-4.2	-4.4	-5.3	-6.7	-6.6	-5.8	-6.4	-6.7 ;
-2.4	-3.8	-5.4	-6.5	-6.6	-6.6	-6.7	-6.7	-6.9	-7.8	-9.1	-9.2	-8.7	-9.0	-9.4 ;
-4.9	-6.2	-7.8	-9.0	-9.1	-9.1	-9.1	-9.2	-9.5	-10.6	-12.3	-12.2	-11.4	-11.6	-11.9 ;
-7.1	-8.9	-10.0	-10.5	-10.6	-10.6	-10.8	-10.9	-11.3	-12.5	-13.8	-13.6	-12.8	-12.3	-12.5 ;
-5.8	-7.7	-9.0	-9.7	-10.0	-10.2	-10.6	-10.8	-11.3	-12.5	-13.8	-13.6	-12.6	-12.0	-12.2 ;
-3.1	-5.1	-6.9	-8.0	-8.3	-8.4	-8.8	-9.1	-9.9	-11.2	-13.2	-13.0	-11.9	-11.5	-11.4 ;
-5.9	-6.8	-7.0	-6.8	-7.0	-7.2	-7.7	-8.1	-8.9	-10.7	-12.9	-12.7	-11.7	-11.3	-11.1 ;
-7.3	-6.6	-6.0	-6.2	-6.6	-6.9	-7.0	-7.2	-8.0	-10.0	-11.9	-11.9	-11.1	-11.0	-10.9 ;
-5.3	-4.5	-4.3	-5.1	-5.8	-6.0	-6.0	-6.1	-6.7	-8.1	-9.3	-9.0	-8.6	-8.8	-8.9 ;
-3.1	-2.8	-2.7	-3.4	-3.9	-4.2	-4.4	-4.5	-4.8	-5.2	-5.8	-5.7	-5.6	-5.6	-5.5 ;
-1.4	-1.3	-1.3	-1.6	-1.9	-2.1	-2.3	-2.4	-2.4	-2.5	-2.6	-2.6	-2.6	-2.6	-2.6 ;
];

T3 = [
0	0	0	0	0	0	0	0	0	0	0	0	0 ;
2.6	2.6	2.4	2.3	2.1	2.1	2.1	2.2	2.3	2.4	2.3	2.0	1.2 ;
5.1	5.0	4.8	4.4	4.0	3.8	3.7	3.7	3.7	3.6	3.3	2.8	1.6 ;
7.4	7.3	7.0	6.3	5.6	5.1	4.7	4.4	4.4	4.3	3.7	3.0	1.6 ;
8.9	8.9	8.5	7.9	7.0	6.3	5.6	5.0	4.7	4.4	4.0	3.3	2.3 ;
9.9	10.0	9.7	9.1	8.1	7.2	6.4	5.9	5.4	5.2	5.1	5.0	5.0 ;
10.1	1.03	10.0	9.5	8.7	7.7	7.0	6.4	6.2	6.4	6.8	7.4	8.2 ;
8.8	9.2	9.0	8.6	7.9	7.2	6.7	6.4	6.4	6.6	7.3	8.0	9.0 ;
6.3	7.0	7.2	6.9	6.3	5.9	5.5	5.5	5.7	6.0	6.7	7.3	8.0 ;
3.3	4.0	4.6	4.3	3.8	3.4	3.3	3.6	4.0	4.5	5.1	5.6	6.2 ;
0.1	0.7	1.8	1.6	1.0	0.7	0.8	1.1	1.6	2.1	2.6	3.1	3.7 ;
-2.7	-2.2	-0.8	-1.0	-1.4	-1.8	-1.8	-1.5	-1.0	-0.7	-0.2	-0.2	0.5 ;
-4.8	-4.5	-3.2	-3.2	-3.8	-4.2	-4.4	-4.4	-4.0	-3.7	-3.3	-2.9	-2.5 ;
-6.6	-6.4	-5.7	-5.6	-5.9	-6.4	-6.8	-7.0	-6.9	-6.7	-6.3	-5.9	-5.2 ;
-9.3	-9.0	-8.1	-8.0	-8.2	-8.8	-9.3	-9.6	-9.7	-9.6	-9.4	-9.0	-8.5 ;
-11.7	-11.4	-10.6	-10.2	-10.5	-11.0	-11.9	-12.5	-12.9	-13.0	-12.9	-12.4	-11.6 ;
-12.4	-12.1	-11.5	-11.2	-11.5	-12.3	-13.3	-14.1	-14.6	-14.8	-14.8	-14.3	-13.9 ;
-12.0	-11.7	-10.7	-10.4	-10.9	-11.9	-13.0	-13.9	-14.5	-14.8	-15.0	-14.8	-14.5 ;
-11.1	-10.6	-9.3	-8.9	-9.0	-9.6	-10.7	-12.0	-13.1	-13.7	-14.0	-14.2	-14.2 ;
-10.8	-10.4	-8.8	-8.0	-8.0	-8.9	-10.0	-11.0	-11.8	-12.1	-12.3	-12.2	-12.0 ;
-10.6	-10.2	-8.8	-8.2	-8.3	-9.0	-10.0	-11.0	-11.3	-11.3	-11.0	-10.0	-9.0 ;
-8.7	-8.4	-7.8	-7.6	-7.7	-8.7	-9.8	-10.7	-10.8	-10.4	-9.6	-8.1	-6.5 ;
-5.5	-5.5	-5.4	-5.2	-5.0	-5.4	-6.3	-7.0	-7.3	-7.0	-6.4	-5.1	-4.0 ;
-2.6	-2.6	-2.6	-2.5	-2.4	-2.4	-2.7	-3.0	-3.1	-3.0	-2.8	-2.3	-1.6 ;
];

FF1 = [
0.5	1.0	1.3	1.4	1.5	1.8	2.3	2.4	2.8	3.1	3.0	2.6	2.7	3.0	4.1
];

FF2 = [
6.1	9.0	12.0	15.9	16.8	16.8	15.8	15.4	14.9	14.7	14.3	12.8	10.7	8.9	7.3
];

FF3 = [
6.4	5.8	4.3	3.1	1.8	0.5	-0.6	-1.7	-1.7	2.5	6.8	8.4	8.5
];

F1 = [
0.2	0.25	0.3	0.32	0.4	0.5	0.6	0.63	0.7	0.8	0.9	1.0	1.2	1.25	1.4
];

F2 = [
1.6	1.8	2.0	2.3	2.5	2.7	2.9	3.0	3.2	3.5	4.0	4.5	5.0	5.5	6.0
];

F3 = [
6.3	6.5	7.0	7.5	8.0	8.5	9.0	9.5	10.0	10.5	11.0	11.5	12.0
];

a = [0:15:345];			% Azimuth

D = [T1 T2 T3];			% Data matrix from Shaw
FF = [FF1 FF2 FF3];		% Free field to eardrum at 0 azimuth
F = [F1 F2 F3];			% Frequency

DF = mean(D)+FF;		% "Diffuse field response" 
semilogx(F,DF,'-',F,FF,'x');