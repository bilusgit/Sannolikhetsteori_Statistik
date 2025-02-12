%% Problem 1: Simulering av konfidensintervall
% Parametrar:
n = 100; %Antal matningar
mu = 2; %Vantevardet
sigma = 1; %Standardavvikelsen
alpha = 0.05;
%Simulerar n observationer for varje intervall
x = normrnd(mu, sigma,n,100); %n x 100 matris med varden
%Skattar mu med medelvardet
xbar = mean(x); %vektor med 100 medelvarden.
%Beraknar de undre och ovre granserna
undre = xbar - norminv(1-alpha/2)*sigma/sqrt(n);
ovre = xbar + norminv(1-alpha/2)*sigma/sqrt(n);

%% Problem 1: Simulering av konfidensintervall (forts.)
%Ritar upp alla intervall
figure(1)
hold on
for k=1:100
if ovre(k) < mu % Rodmarkerar intervall som missar mu
plot([undre(k) ovre(k)],[k k],'r')
elseif undre(k) > mu
plot([undre(k) ovre(k)],[k k],'r')
else
plot([undre(k) ovre(k)],[k k],'b')
end
end
%b1 och b2 ar bara till for att figuren ska se snygg ut.
b1 = min(xbar - norminv(1 - alpha/2)*sigma/sqrt(n));
b2 = max(xbar + norminv(1 - alpha/2)*sigma/sqrt(n));
axis([b1 b2 0 101]) %Tar bort outnyttjat utrymme i figuren
%Ritar ut det sanna vardet
plot([mu mu],[0 101],'g')
hold off

%% Problem 2: Maximum likelihood/Minsta kvadrat
M = 1e4;
b = 4;
x = raylrnd(b, M, 1);
hist_density(x, 40)
hold on
est_ml = sqrt(sum(x.*x)/(2*M)); 
est_mk = (sum(x)/M)*sqrt(2/pi);
plot(est_ml, 0, 'r*')
plot(est_mk, 0, 'g*')
plot(b, 0, 'ro')

%% Problem 2: Maximum likelihood/Minsta kvadrat (forts.)
plot(0:0.5:20, raylpdf(0:0.5:20, est_ml), 'r')
hold off

%% Problem 3: Konfidensintervall for Rayleighfordelning
y = load('wave_data.dat');
subplot(2,1,1), plot(y(1:end))
subplot(2,1,2), hist_density(y)

y_mean = mean(y)
n_y = length(y)
est_ml = sqrt(sum(y.*y)/(2*n_y));
est = (sum(y)/n_y)*sqrt(2/pi);

alpha = 0.05
s = sqrt(sum((y - y_mean).^2) / (n_y - 1)); % samma sak som std(y)?
t = tinv(1 - alpha/2, n_y - 1);

mu_lower = mean(y) - t * (s / sqrt(n_y));
mu_upper = mean(y) + t * (s / sqrt(n_y));


lower_bound = mu_lower / sqrt(pi/2);
upper_bound = mu_upper / sqrt(pi/2);

%% Problem 3: Konfidensintervall (forts.)
hold on
plot(lower_bound, 0, 'g*')
plot(upper_bound, 0, 'g*')
plot(est, 0, "ro")

%% Problem 3: Konfidensintervall (forts.)
plot(0:0.1:6, raylpdf(0:0.1:6, est), 'r')
hold off

%% Problem 4:
data = load("birth.dat")

barn_vikt = data(:, 3)
moder_alder = data(:, 4)
moder_langd = data(:, 16)
moder_vikt = data(:, 15)

figure; 
subplot(2, 2, 1); 
hist_density(barn_vikt);
title('Barnets vikt');

subplot(2, 2, 2);  
hist_density(moder_alder)
title('Moderns ålder');

subplot(2, 2, 3);
hist_density(moder_langd)
title('Moderns längd');

subplot(2, 2, 4);
hist_density(moder_vikt)
title('Moderns vikt');

%% Problem 4: Fordelningar av givna data
x = data(data(:, 20) < 3, 3); % nonsmoking
y = data(data(:, 20) == 3, 3); % smoking 

size_x = size(x) % antal mammor som inte röker
size_y = size(y) % antal som röker

%% Problem 4: Fordelningar av givna data (forts.)
subplot(2,2,1), boxplot(x); % till vänster: icke-rökande
axis([0 2 500 5000]);
subplot(2,2,2), boxplot(y); % till höger: rökande
axis([0 2 500 5000]);
%% Problem 4: Fordelningar av givna data (forts.)
subplot(2,2,3:4), ksdensity(x); % blå = icke-rökande
hold on
[fy, ty] = ksdensity(y);
plot(ty, fy, 'r'); % röd = rökande
hold off

%% Problem 4: fortsättning - förhållande mellan underviktiga barn och underviktiga moder?
height = data(:, 16)*10^(-2) % längd i m
weight = data(:, 15) % vikt i kg

BMI = weight./(height.^2)
% weight/(height*10^-2)^2
% overweight if 25 <= BMI < 30
overweight = BMI >=25

barn_vikt_not_overweight = data(BMI < 25, 3)
barn_vikt_overweight = data(BMI >= 25, 3)

subplot(2,2,1), boxplot(barn_vikt_not_overweight); % till vänster: inte övervikt
axis([0 2 500 5000]);
subplot(2,2,2), boxplot(barn_vikt_overweight); % till höger: övervikt
axis([0 2 500 5000]);
size(overweight)

%% Problem 5: undersökning om moder vikt, längd, ålder och barn vikt normalfördelade?
mother_height = data(:, 16)*10^(-2) % längd i m
mother_weight = data(:, 15) % vikt i kg
mother_age = data(:, 4) % moderns ålder i år
baby_weight = data(:, 3) % barnets vikt i gram

figure; 
subplot(2, 2, 1); 
normplot(baby_weight);
title('Barnets vikt');

subplot(2, 2, 2);  
normplot(mother_age)
title('Moderns ålder');

subplot(2, 2, 3);
normplot(mother_height)
title('Moderns längd');

subplot(2, 2, 4);
normplot(mother_weight)
title('Moderns vikt');

figure; 
subplot(2, 2, 1); 
qqplot(baby_weight);
title('Barnets vikt');

subplot(2, 2, 2);  
qqplot(mother_age)
title('Moderns ålder');

subplot(2, 2, 3);
qqplot(mother_height)
title('Moderns längd');

subplot(2, 2, 4);
qqplot(mother_weight)
title('Moderns vikt');

%% Problem 5: forts
jbtest(baby_weight, 0.05) % 1
jbtest(mother_age, 0.05) % 1
jbtest(mother_height, 0.05) % 0
jbtest(mother_weight, 0.05) % 1
%% Problem 6: linjär regression
moore = load("moore.dat");

plot(moore(:,1), moore(:,2));
xlabel('Årtal');
ylabel('Transistorer per ytenhet');

%% plotta exponentiella regressionen
figure
scatter(moore(:,1), moore(:,2)); % lägg punkter i x-värden som alla rader i första kolumnen mot y-värden som alla rader i andra kolumnen

w = log(moore(:,2)); % Logatimera antalet transistorer per ytenhet
ones_column = ones(length(moore(:,1)), 1) % skapar en nx1 matris av ettor
X = [ones_column, moore(:,1)]; % skapar en matris vars första kolumn är ettor and andra kolumn är årtal
[b, bint, r] = regress(w, X); 

hold on;
y_fit = X * b; % Beräkna de skattade y-värdena
plot(moore(:,1), exp(y_fit), 'r'); % Exponera y_fit för att ta bort logatimeringen
hold off;
%% plotta linjärt
figure;
hold on;
scatter(moore(:,1), log(moore(:,2))); % lägg punkter i x-värden som alla rader i första kolumnen mot y-värden som alla rader i andra kolumnen

w = log(moore(:,2)); % Logatimera antalet transistorer per ytenhet
ones_column = ones(length(moore(:,1)), 1) % skapar en nx1 matris av ettor
X = [ones_column, moore(:,1)]; % skapar en matris vars första kolumn är ettor and andra kolumn är årtal
[b, bint, r] = regress(w, X); 


y_fit = X * b; % Beräkna de skattade y-värdena
plot(moore(:,1), y_fit, 'r'); % Exponera y_fit för att ta bort logatimeringen
hold off;
%% Problem 6: Regression
res = X * b - w;
subplot(2,1,1), normplot(res);
subplot(2,1,2), hist(res);
jbtest(res, 0.05)
%% Problem 6: forts. prediktion år 2025
x_2025 = [1, 2025] % samma format som skapade x:en
y_prediction = x_2025*b

exp(y_prediction) % exponera för att få bort logatimeringen => y_pred = 1.3599e+08
%% Problem 7: multipel linjär regression av barnvikt
% Vi använde den kategoriska variabeln moderns vikt och omvandlade den till
% en BMI skala: om BMI >= 25 sätt 1 annars 0 (skala för övervikt)
data = load("birth.dat")

mother_height = data(:, 16) % längd i cm
mother_height_m = data(:, 16)*10^(-2) % längd i m
mother_weight = data(:, 15) % vikt i kg
mother_age = data(:, 4) % moderns ålder i år
baby_weight = data(:, 3) ./ 1000 % barnets vikt i kg
mothers_smoking = data(:, 20) == 3

BMI = mother_weight./(mother_height_m.^2)

overweight = BMI >=25
%% Enkel modell hur barnets vikt beror på moderns längd
X = [ones(length(mother_height), 1), mother_height]; % Precis som förut, en massa ettor på första kolumnen och moder längd i andra kol

[b,bint,r,rint,stats] = regress(baby_weight, X); % källa: https://se.mathworks.com/help/stats/regress.html

figure;
scatter(mother_height, baby_weight); % rita ut moderlängd på x och barnvikt på y som punkter
hold on;
y_fit = X * b;
plot(mother_height, y_fit)
hold off;
%% multipel moderns vikt, moderns rökvanor och BMI
Z = [ones(length(mothers_smoking), 1), mother_weight, mothers_smoking, overweight]; % alla vektorer har samma "längd" (747)
[b,bint,r,rint,stats] = regress(baby_weight,Z);
normplot(r);
jbtest(r)

% 'regress treats NaN values as missing values. regress omits observations with missing values from the regression fit.'
% bint - konfidensintervall för de skattade koefficienterna
% 2.1516    2.7791 - intervallet för vart b0, dvs interceptet ligger
% 0.0114    0.0219 - intervall för moderns vikt - positiv koefficient
% -0.2477   -0.0712 - intervall för rökandet, minskning
% -0.4032   -0.0945 - intervall för övervikt, minskning?

bint
stats

% stats = 
%   R2 värde  f-statistic  p-värde standrad error
%    0.0666   17.1715    0.0000    0.3040

%% Utför ytterliggare test
data_new = table(baby_weight, mother_weight, mothers_smoking, overweight, ...
    'VariableNames', {'BabyWeight', 'MotherWeight', 'SmokingStatus', 'Overweight'});
 
% Fit linear model
modelspec = 'BabyWeight ~ MotherWeight + SmokingStatus + Overweight';
mdl = fitlm(data_new, modelspec);
 
% Plotting
plot(mdl);
mdl 
