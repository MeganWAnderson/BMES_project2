%BMES Anderson Project 2 Code 
close all; clear all; clc; 
%% Load and Sort Data 
raw_data = readtable('framinghamData.xls'); 

%Sort out P1 data 
rows = height(raw_data);
count1 = 1; 
count2 = 1;
count3 = 1; 
for i = 1:rows
    if raw_data.PERIOD(i) == 1 
        p1_data(count1,:) = raw_data(i,:);
        count1 = count1 + 1;
    end
end 

%Sort Men and Women 
rows = height(p1_data); 
count1 = 1; 
count2 = 1; 
for i = 1:rows
    if p1_data.SEX(i) == 1 %1 = men 
        men_data1(count1,:) = p1_data(i,:);
        count1 = count1 + 1;  
    else %2 = women 
        women_data1(count2,:) =p1_data(i,:); 
        count2 = count2+1; 
end 
end

%% Part 1 - Linear Regression 
%#1 raw or log data?

%Fit Values - Normal
model_p1 = fitlm(p1_data.BMI, p1_data.SYSBP) %creates model
fit_vals_p1 = predict(model_p1, p1_data.BMI); %creates predicted vals fit to model
residuals = p1_data.SYSBP - fit_vals_p1; %calculatred residuals based on model

%Fit Values - Log
model_p1_log = fitlm(log(p1_data.BMI), log(p1_data.SYSBP))
fit_vals_p1_log = predict(model_p1_log, log(p1_data.BMI));
residuals_log = p1_data.SYSBP - fit_vals_p1_log; 

%Linearity and Homoscedasticity
%EQ for y = 0 line
%for normal  data
x=[100:0.1:200];
y = zeros(length(x),1);

%for log data
x2=[4:0.1:6];
y2 = zeros(length(x2),1);

% Raw Data  
figure(1)
scatter(fit_vals_p1,residuals) %residuals vs. fitted data 
hold on
plot(x,y) %line y = 0 for comparison 
hold off
xlabel('Fitted Vals')
ylabel('Residuals')
title('Residuals v. Fitted Vals for Normal Data')

%Log Data
figure(2)
scatter(fit_vals_p1_log,residuals_log)
hold on
plot(x2,y2)
hold off
xlabel('Fitted Vals')
ylabel('Residuals')
title('Residuals v. Fitted Vals for Log Data')

%Normality of Residuals 
figure()
subplot(1,2,1)
probplot('normal', residuals)
xlabel('Raw Data')
subplot(1,2,2)
probplot('normal', residuals_log)
xlabel('Natural log of Raw Data')

%#2 Linear Models - men and women
%take out NaN values from data sets 
men_data_index = ~isnan(men_data1.BMI) & ~isnan(men_data1.SYSBP);
women_data_index = ~isnan(women_data1.BMI) & ~isnan(women_data1.SYSBP);

data_men_p1 = men_data1(men_data_index,:);
data_women_p1 = women_data1(women_data_index,:);

%men 
p1_men_model = fitlm(data_men_p1.BMI, data_men_p1 .SYSBP) %creates model 
fit_vals_men = predict(p1_men_model, data_men_p1.BMI); %creates predicted values fit to model 

%women 
p1_women_model = fitlm(data_women_p1.BMI, data_women_p1 .SYSBP)
fit_vals_women = predict(p1_women_model, data_women_p1.BMI); 


%#3 Plot Regressions 
figure()
subplot (1,2,1)
hold on
scatter(data_men_p1.BMI, data_men_p1.SYSBP) %scatter plot of data 
plot(data_men_p1.BMI, fit_vals_men) %line representing linear regression 
xlabel('BMI (kg/m^2)')
ylabel('SBP (mmHg)')
legend('Men Data', 'Men Fit')
title('SBP (mmHg) vs. BMI (kg/m^2) Regression (Men)')

subplot (1,2,2)
hold on
scatter(data_women_p1.BMI, data_women_p1.SYSBP)
plot(data_women_p1.BMI, fit_vals_women)
xlabel('BMI (kg/m^2)')
ylabel('SBP (mmHg)')
legend('Women Data', 'Women Fit')
title('SBP (mmHg) vs BMI (kg/m^2) Regression (Women)')

%#4 Analyze Assumptions 
%linearity and homoscedasticity
%men
figure()
subplot(2, 1, 1);
plotResiduals(p1_men_model, 'probability'); %normality of residuals 
title('Normality of Residuals (Men)');

subplot(2, 1, 2);
plotResiduals(p1_men_model, 'fitted'); %residuals vs. fitted vals 
title('Residuals vs. Fitted Values (Men)'); 

%women
figure()
subplot(2, 1, 1);
plotResiduals(p1_men_model, 'probability');
title('Normality of Residuals (Women)');

subplot(2, 1, 2);
plotResiduals(p1_men_model, 'fitted');
title('Residuals vs. Fitted Values (Women)'); 


%#5 Predictions of man and woman of BMI = 33
SBP_predict_men = predict(p1_men_model, 33)
SBP_predict_women = predict(p1_women_model, 33)

%% Part 2 - Multiple Regression 
predictors = {'AGE', 'BMI', 'GLUCOSE', 'TOTCHOL'}; 
%predictors = age, BMI, glucose levels, totoal cholesterol levels 
xvals = men_data1{:, predictors};
yvals= men_data1.SYSBP;

xvals2 = women_data1{:, predictors};
yvals2= women_data1.SYSBP;

%create regression models 
mult_regress_men = fitlm(xvals, yvals, 'VarNames', [predictors, {'SBP'}])
mult_regress_women = fitlm(xvals2, yvals2, 'VarNames', [predictors, {'SBP'}])

%Find Fit Vals 
fit_vals_mult_men = predict(mult_regress_men, xvals); 
fit_vals_mult_women = predict(mult_regress_women, xvals2); 

%Vizualization - plot matrix 
%men 
figure()
plotmatrix([xvals yvals]) 
title('Plot Matrix - Men Predictors and SBP')

%women
figure()
plotmatrix([xvals2 yvals2])
title('Plot Matrix - Women Predictors and SBP')

% normality, homoscedasciity, and linearity plots - Men 
figure()
subplot(2, 1, 1);
plotResiduals(mult_regress_men, 'probability'); %normality of residuals 
title('Normality of Residuals Men');

subplot(2, 1, 2);
plotResiduals(mult_regress_men, 'fitted'); %residuals vs. fitted vals 
title('Homoscedasticity of Residuals Men'); 

% Normality and Homoscedasciity plots - Women 
figure()
subplot(2, 1, 1);
plotResiduals(mult_regress_women, 'probability');
title('Normality of Residuals Women');

subplot(2, 1, 2);
plotResiduals(mult_regress_women, 'fitted');
title('Homoscedasticity of Residuals Women');

