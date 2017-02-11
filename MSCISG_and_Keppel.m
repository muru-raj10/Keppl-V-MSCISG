                  
clear;
clc;
format compact;

%Missing data on row 1043, 9th March 1992 has been input manually as 2.4936
%for keppel corp
%MSDLSGF data is more accurate in the
%sense that there isn't any repeat values on days that the market is
%closed. Hence, we matched the two data sets on the days of MSDLSGF to
%eliminate all the excess days where Keppel corp wasn't trading. 
%note that some of the days there was trading and they have similar closing prices.

%%% Read in dates into d1 and prices into p1
filename1 = 'MSDLSGF--Index.csv';
[d1 , p1] = textread(filename1, '%s %f', 'delimiter', ',', 'headerlines', 1);
filename2 = 'KEP-SP-Equity.csv';
[d1e , p1e] = textread(filename2, '%s %f', 'delimiter', ',', 'headerlines', 1);

%let e or eq at the end of the variable to denote keppel corp equity
d1n=datenum(d1,'dd/mm/yyyy');
d1en=datenum(d1e,'dd/mm/yyyy');  

Ind = table(d1n,d1,p1,'VariableNames',{'DateNum','Date','MSDLSGF'});
Eq = table(d1en,d1e,p1e,'VariableNames',{'DateNum','Date','KeppelC'});

DT=innerjoin(Ind,Eq);  %to sort out the data and allign prices for msdlsg and keppel
[WeekNbr]=weeknum(DT.DateNum);
[MthNbr,MthString]=month(DT.DateNum);

WkMthTbl=table(WeekNbr,MthNbr,MthString,'VariableNames', {'WeekNbr','MonthNbr','MonthString'});

DT=[DT WkMthTbl];

ts = timeseries;
ts.Data = DT.MSDLSGF;
ts = setabstime(ts, DT.Date);
ts.Name = ' Daily MSCI SG Index';
ts.TimeInfo.Format = 'mmmyy';

tse = timeseries;
tse.Data = DT.KeppelC;
tse = setabstime(tse, DT.Date);
tse.Name = ' Daily Kep Corp equity price';
tse.TimeInfo.Format = 'mmmyy';

%figure (1);
FigHandle1 = figure('Position', [40, 20, 1200, 600]);
plot(ts, 'blue');
grid on;

%figure (2);
FigHandle2 = figure('Position', [40, 20, 1200, 600]);
plot(tse, 'blue');
grid on;

%%% Plot the daily log return  for mscisg%%%
logprice = log(DT.MSDLSGF);
logret = diff(logprice);
ts_r = timeseries;
ts_r.Data = logret;
ts_r = setabstime(ts_r, DT.Date(2:end));
ts_r.Name = ' Daily Log Return on MSCI SG Index';
ts_r.TimeInfo.Format = 'mmmyy';

%%% Plot the daily log return  for keppel corp%%%
logpriceEq = log(DT.KeppelC);
logretEq = diff(logpriceEq);
ts_rEq = timeseries;
ts_rEq.Data = logretEq;
ts_rEq = setabstime(ts_rEq, DT.Date(2:end));
ts_rEq.Name = ' Daily Log Return on Kep Corp Equity';
ts_rEq.TimeInfo.Format = 'mmmyy';

FigHandle3 = figure('Position', [100, 30, 1200, 600]);
plot(ts_r, 'red');
grid on;

FigHandle4 = figure('Position', [100, 30, 1200, 600]);
plot(ts_rEq, 'red');
grid on;

L_ind=length(logret);         %nbr or elements
L_eq=length(logretEq);         %nbr or elements should be the same

%moments
momentsDaily=table;
momentsDaily.DailyLogReturnsInd=[mean(logret); std(logret,1); skewness(logret); kurtosis(logret)];
momentsDaily.DailyLogReturnsEq=[mean(logretEq); std(logretEq,1); skewness(logretEq); kurtosis(logretEq)];
momentsDaily.Properties.RowNames={'Mean', 'SD', 'Skewness', 'Kurtosis'};

%Jarque Bera test statistic
JB_ind_Daily=(L_ind/6)*(momentsDaily.DailyLogReturnsInd(3))^2+(L_ind/24)*(momentsDaily.DailyLogReturnsInd(4)-3)^2;
JB_eq_Daily=(L_eq/6)*(momentsDaily.DailyLogReturnsEq(3))^2+(L_eq/24)*(momentsDaily.DailyLogReturnsEq(4)-3)^2;

%H_0=Data is normally distributed, H_1 = Data is not normal
%reject H_0 if JB > chi_sq(2,0.95). reject H_0 if test = 1
[Nor_hyp_indDaily,pval_indDaily]=jbtest(logret);
[Nor_hyp_eqDaily,pval_eqDaily]=jbtest(logretEq);

%for correlation coefficient
StDevDaily=table;
StDevDaily.DailyLogReturnsInd=sqrt((sum((logret-momentsDaily.DailyLogReturnsInd(1)).^2))/(L_ind-1));
StDevDaily.DailyLogReturnsEq=sqrt((sum((logretEq-momentsDaily.DailyLogReturnsEq(1)).^2))/(L_eq-1));
StDevDaily.DailyCovariance=(sum((logret-momentsDaily.DailyLogReturnsInd(1)).*(logretEq-momentsDaily.DailyLogReturnsEq(1))))/(L_ind-1);
StDevDaily.DailyCorrelation= StDevDaily.DailyCovariance(1)/(StDevDaily.DailyLogReturnsInd(1)*StDevDaily.DailyLogReturnsEq(1));

r_capDaily=StDevDaily.DailyCorrelation(1);
%for hypothesis testing
t_statDaily=sqrt(L_ind-2)*r_capDaily/(sqrt(1-r_capDaily^2));
%H_0 : rho =0, H_a = rho != 0, 
%At alpha=5%, reject H_0 when t_stat>t(0.975,L_ind-2), reject H_0 if pval
%is small
[rho_capDaily, pvalDaily]=corrcoef(logret,logretEq);

%for regression coefficient
mdl=fitlm(logret,logretEq);

%% Linear Regression Plots
% Daily Returns
figure
plot(logret,logretEq,'k.')
hold on
xhatd=linspace(min(logret),max(logret));
yhatd= mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2)*xhatd;
plot(xhatd,yhatd,'r-');
title('Plot of Daily Log Returns of Keppel Corp Ltd against MSCI Singapore');
xlabel('Log Returns (MSCI Singapore)');
ylabel('Log Returns (Keppel Corp Ltd.');
%saveas(gcf,'Fig 7 - Linear Regression Daily')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MSDWeek=zeros(L_ind,1);
for i = 1 : L_ind
    if DT.WeekNbr(i+1)==DT.WeekNbr(i)
        MSDWeek(i)=0;
    else
        MSDWeek(i)=DT.MSDLSGF(i);
    end
end

WkIndex=find(MSDWeek);
MSDWeeklyPrice=MSDWeek(WkIndex);
KepEqWeeklyPrice=DT.KeppelC(WkIndex);

Wklogprice = log(MSDWeeklyPrice);
Wklogret = diff(Wklogprice);

WklogpriceEq = log(KepEqWeeklyPrice);
WklogretEq = diff(WklogpriceEq);

tsW = timeseries;
tsW.Data = MSDWeeklyPrice;
tsW = setabstime(tsW, DT.Date(WkIndex));
tsW.Name = ' Weekly MSCI SG Index';
tsW.TimeInfo.Format = 'mmmyy';

tsWe = timeseries;
tsWe.Data = KepEqWeeklyPrice;
tsWe = setabstime(tsWe, DT.Date(WkIndex));
tsWe.Name = ' Weekly Kep Corp equity price';
tsWe.TimeInfo.Format = 'mmmyy';

%figure (5);
FigHandle5 = figure('Position', [40, 20, 1200, 600]);
plot(tsW, 'blue');
grid on;

%figure (6);
FigHandle6 = figure('Position', [40, 20, 1200, 600]);
plot(tsWe, 'blue');
grid on;

%%% Plot the weekly log return  for mscisg%%%
ts_Wr = timeseries(Wklogret,datestr(DT.DateNum(WkIndex(2:end))));
% ts_Wr.Data = Wklogret;
% ts_Wr.Time=WkDate;
%ts_Wr = setabstime(ts_Wr,WkDate);
ts_Wr.Name = ' Weekly Log Return on MSCI SG Index';
ts_Wr.TimeInfo.Format = 'mmmyy';

%%% Plot the weekly log return  for keppel corp%%%
ts_WrEq = timeseries(WklogretEq,datestr(DT.DateNum(WkIndex(2:end))));
% ts_WrEq.Data = WklogretEq;
% ts_WrEq = setabstime(ts_WrEq, WkDate);
ts_WrEq.Name = ' Weekly Log Return on Kep Corp Equity';
ts_WrEq.TimeInfo.Format = 'mmmyy';

FigHandle7 = figure('Position', [100, 30, 1200, 600]);
plot(ts_Wr, 'red');
grid on;

FigHandle8 = figure('Position', [100, 30, 1200, 600]);
plot(ts_WrEq, 'red');
grid on;

L_indWk=length(Wklogret);         %nbr or elements
L_eqWk=length(WklogretEq);         %nbr or elements should be the same

momentsWeekly=table;
momentsWeekly.WeeklyLogReturnsInd=[mean(Wklogret); std(Wklogret,1); skewness(Wklogret); kurtosis(Wklogret)];
momentsWeekly.WeeklyLogReturnsEq=[mean(WklogretEq); std(WklogretEq,1); skewness(WklogretEq); kurtosis(WklogretEq)];
momentsWeekly.Properties.RowNames={'Mean', 'SD', 'Skewness', 'Kurtosis'};

%Jarque Bera test statistic
JB_ind_Weekly=(L_indWk/6)*(momentsWeekly.WeeklyLogReturnsInd(3))^2+(L_indWk/24)*(momentsWeekly.WeeklyLogReturnsInd(4)-3)^2;
JB_eq_Weekly=(L_eqWk/6)*(momentsWeekly.WeeklyLogReturnsEq(3))^2+(L_eqWk/24)*(momentsWeekly.WeeklyLogReturnsEq(4)-3)^2;

%H_0=Data is normally distributed, H_1 = Data is not normal
%reject H_0 if JB > chi_sq(2,0.95). reject H_0 if test = 1
[Nor_hyp_indWeekly,pval_indWeekly]=jbtest(Wklogret);
[Nor_hyp_eqWeekly,pval_eqWeekly]=jbtest(WklogretEq);

%for correlation coefficient
StDevWeekly=table;
StDevWeekly.WeeklyLogReturnsInd=sqrt((sum((logret-momentsWeekly.WeeklyLogReturnsInd(1)).^2))/(L_indWk-1));
StDevWeekly.WeeklyLogReturnsEq=sqrt((sum((logretEq-momentsWeekly.WeeklyLogReturnsEq(1)).^2))/(L_eqWk-1));
StDevWeekly.WeeklyCovariance=(sum((logret-momentsWeekly.WeeklyLogReturnsInd(1)).*(logretEq-momentsWeekly.WeeklyLogReturnsEq(1))))/(L_indWk-1);
StDevWeekly.WeeklyCorrelation= StDevWeekly.WeeklyCovariance(1)/(StDevWeekly.WeeklyLogReturnsInd(1)*StDevWeekly.WeeklyLogReturnsEq(1));

r_capWeekly=StDevWeekly.WeeklyCorrelation(1);
%for hypothesis testing
t_statWeekly=sqrt(L_indWk-2)*r_capWeekly/(sqrt(1-r_capWeekly^2));
%H_0 : rho =0, H_a = rho != 0, 
%At alpha=5%, reject H_0 when t_stat>t(0.975,L_ind-2), reject H_0 if pval is small
[rho_capWeekly, pvalWeekly]=corrcoef(Wklogret,WklogretEq);

Wkmdl=fitlm(Wklogret,WklogretEq);

% Weekly Returns
figure
plot(Wklogret,WklogretEq,'k.')
hold on
xhatw=linspace(min(Wklogret),max(WklogretEq));
yhatw= Wkmdl.Coefficients.Estimate(1) + Wkmdl.Coefficients.Estimate(2)*xhatw;
plot(xhatw,yhatw,'r-');
title('Plot of Weekly Log Returns of Keppel Corp Ltd against MSCI Singapore');
xlabel('Log Returns (MSCI Singapore)');
ylabel('Log Returns (Keppel Corp Ltd.');
%saveas(gcf,'Fig 8 - Linear Regression Weekly')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MSDMonth=zeros(L_ind,1);
for i = 1 : L_ind
    if DT.MonthNbr(i+1)==DT.MonthNbr(i)
        MSDMonth(i)=0;
    else
        MSDMonth(i)=DT.MSDLSGF(i);
    end
end

MonthIndex=find(MSDMonth);
MSDMonthlyPrice=MSDMonth(MonthIndex);
KepEqMonthlyPrice=DT.KeppelC(MonthIndex);

Mthlogprice = log(MSDMonthlyPrice);
Mthlogret = diff(Mthlogprice);

MthlogpriceEq = log(KepEqMonthlyPrice);
MthlogretEq = diff(MthlogpriceEq);


tsM = timeseries;
tsM.Data = MSDMonthlyPrice;
tsM = setabstime(tsM, datestr(DT.DateNum(MonthIndex)));
tsM.Name = ' Monthly MSCI SG Index';
tsM.TimeInfo.Format = 'mmmyy';

tsMe = timeseries;
tsMe.Data = KepEqMonthlyPrice;
tsMe = setabstime(tsMe, datestr(DT.DateNum(MonthIndex)));
tsMe.Name = ' Monthly Kep Corp equity price';
tsMe.TimeInfo.Format = 'mmmyy';

%figure (9);
FigHandle9 = figure('Position', [40, 20, 1200, 600]);
plot(tsM, 'blue');
grid on;

%figure (10);
FigHandle10 = figure('Position', [40, 20, 1200, 600]);
plot(tsMe, 'blue');
grid on;

%%% Plot the monthly log return  for mscisg%%%
ts_Mr = timeseries(Mthlogret,datestr(DT.DateNum(MonthIndex(2:end))));
% ts_Mr.Data = Mthlogret;
% ts_Mr = setabstime(ts_Wr, DT.Date(MonthIndex(2:end)));
ts_Mr.Name = ' Monthly Log Return on MSCI SG Index';
ts_Mr.TimeInfo.Format = 'mmmyy';

%%% Plot the monthly log return  for keppel corp%%%
ts_MrEq = timeseries(MthlogretEq,datestr(DT.DateNum(MonthIndex(2:end))));
% ts_MrEq.Data = MthlogretEq;
% ts_MrEq = setabstime(ts_MrEq,DT.Date(MonthIndex(2:end)));
ts_MrEq.Name = ' Monthly Log Return on Kep Corp Equity';
ts_MrEq.TimeInfo.Format = 'mmmyy';

FigHandl1 = figure('Position', [100, 30, 1200, 600]);
plot(ts_Mr, 'red');
grid on;

FigHandle12 = figure('Position', [100, 30, 1200, 600]);
plot(ts_MrEq, 'red');
grid on;



L_indMth=length(Mthlogret);         %nbr or elements
L_eqMth=length(MthlogretEq);         %nbr or elements should be the same

momentsMonthly=table;
momentsMonthly.MonthlyLogReturnsInd=[mean(Mthlogret); std(Mthlogret,1); skewness(Mthlogret); kurtosis(Mthlogret)];
momentsMonthly.MonthlyLogReturnsEq=[mean(MthlogretEq); std(MthlogretEq,1); skewness(MthlogretEq); kurtosis(MthlogretEq)];
momentsMonthly.Properties.RowNames={'Mean', 'SD', 'Skewness', 'Kurtosis'};


%Jarque Bera test statistic
JB_ind_Monthly=(L_indMth/6)*(momentsMonthly.MonthlyLogReturnsInd(3))^2+(L_indMth/24)*(momentsMonthly.MonthlyLogReturnsInd(4)-3)^2;
JB_eq_Monthly=(L_eqMth/6)*(momentsMonthly.MonthlyLogReturnsEq(3))^2+(L_eqMth/24)*(momentsMonthly.MonthlyLogReturnsEq(4)-3)^2;

%H_0=Data is normally distributed, H_1 = Data is not normal
%reject H_0 if JB > chi_sq(2,0.95). reject H_0 if test = 1
[Nor_hyp_indMonthly,pval_indMonthly]=jbtest(Mthlogret);
[Nor_hyp_eqMonthly,pval_eqMonthly]=jbtest(MthlogretEq);

%for correlation coefficient
StDevMonthly=table;
StDevMonthly.MonthlyLogReturnsInd=sqrt((sum((logret-momentsMonthly.MonthlyLogReturnsInd(1)).^2))/(L_indMth-1));
StDevMonthly.MonthlyLogReturnsEq=sqrt((sum((logretEq-momentsMonthly.MonthlyLogReturnsEq(1)).^2))/(L_eqMth-1));
StDevMonthly.MonthlyCovariance=(sum((logret-momentsMonthly.MonthlyLogReturnsInd(1)).*(logretEq-momentsMonthly.MonthlyLogReturnsEq(1))))/(L_indMth-1);
StDevMonthly.MonthlyCorrelation= StDevMonthly.MonthlyCovariance(1)/(StDevMonthly.MonthlyLogReturnsInd(1)*StDevMonthly.MonthlyLogReturnsEq(1));

r_capMonthly=StDevMonthly.MonthlyCorrelation(1);
%for hypothesis testing
t_statMonthly=sqrt(L_indMth-2)*r_capMonthly/(sqrt(1-r_capMonthly^2));
%H_0 : rho =0, H_a = rho != 0, 
%At alpha=5%, reject H_0 when t_stat>t(0.975,L_ind-2), reject H_0 if pval is small
[rho_capMonthly, pvalMonthly]=corrcoef(Mthlogret,MthlogretEq);

Mthmdl=fitlm(Mthlogret,MthlogretEq);

% Monthly Returns
figure
plot(Mthlogret,MthlogretEq,'k.')
hold on
xhatm=linspace(min(Mthlogret),max(Mthlogret));
yhatm= Mthmdl.Coefficients.Estimate(1) + Mthmdl.Coefficients.Estimate(2)*xhatm;
plot(xhatm,yhatm,'r-');
title('Plot of Monthly Log Returns of Keppel Corp Ltd against MSCI Singapore');
xlabel('Log Returns (MSCI Singapore)');
ylabel('Log Returns (Keppel Corp Ltd.');
%saveas(gcf,'Fig 9 - Linear Regression Monthly');