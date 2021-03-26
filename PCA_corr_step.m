clear all; close all; clc;

% load model posteriors
C = readtable('HadCRUT5_Cheng.csv');
C = table2array(C(:,2:end-1)); % convert to array, removing title and end-bracket rows
N = readtable('HadCRUT5_NODC.csv');
N = table2array(N(:,2:end-1));
I = readtable('HadCRUT5_noinfilling_Cheng.csv');
I = table2array(I(:,2:end-1));

% define weights
w = [];
w(1:length(C)) = 1/length(C); 
w(end+1:end+length(N)) = 1/length(N);
w(end+1:end+length(I)) = 1/length(I);

% standardize
P = [C I N]';
P = P-mean(P);
P = P./std(P);

% do PCA
[eivecs,score,~,~,pcts,~] = pca(P,'Weights',w);

%%

% plot PCA results

%
figure
subplot(2,3,1)
scatter(1:25,pcts./100,50,[82 105 79]/256,'filled')
hold on
plot(1:25,pcts./100,'color',[82 105 79]/256,'linewidth',2)
set(gca,'fontsize',15,'ticklabelinterpreter','latex','xtick',[1 5:5:25])
ylabel('Variance Explained','interpreter','latex')
xlabel('PC','interpreter','latex')
box on
axis([0 25 -.005 .31])

% only show highest loadings
subplot(2,3,2)
bar(eivecs(1,[2 14 16]),'FaceColor',[82 105 79]/256)
set(gca,'fontsize',15,'ticklabelinterpreter','latex','xtick',[1:25],'xticklabel',{'$\lambda_{fast}^{equil}$','$a_{CO2}$','$\gamma_{aero-SOx}$'})
ylabel('Loading','interpreter','latex')
title('PC 1','interpreter','latex')
box on
xtickangle(90)
axis([.5 3.5 -.8 .8])


subplot(2,3,3)
bar(eivecs(2,[2 3 16 20]),'FaceColor',[82 105 79]/256)
set(gca,'fontsize',15,'ticklabelinterpreter','latex','xtick',[1:25],'xticklabel',{'$\lambda_{fast}^{equil}$','$\lambda_{multidecadal}^{equil}$','$\gamma_{aero-SOx}$','$\gamma_{aero-NH3}$'})
title('PC 2','interpreter','latex')
box on
xtickangle(90)
axis([.5 4.5 -.8 .8])

subplot(2,3,4)
bar(eivecs(3,[2 3 16 20]),'FaceColor',[82 105 79]/256)
set(gca,'fontsize',15,'ticklabelinterpreter','latex','xtick',[1:25],'xticklabel',{'$\lambda_{fast}^{equil}$','$\lambda_{multidecadal}^{equil}$','$\gamma_{aero-SOx}$','$\gamma_{aero-NH3}$'})
ylabel('Loading','interpreter','latex')
title('PC 3','interpreter','latex')
box on
xtickangle(90)
axis([.5 4.5 -.8 .8])

subplot(2,3,5)
bar(eivecs(4,[5 8 12 13]),'FaceColor',[82 105 79]/256)
set(gca,'fontsize',15,'ticklabelinterpreter','latex','xtick',[1:25],'xticklabel',{'$\tau_{multidecadal}$','$\tau_{mixed}$','$\tau_{bottom}$','$I_b$',})
title('PC 4','interpreter','latex')
box on
xtickangle(90)
axis([.5 4.5 -.8 .8])

subplot(2,3,6)
bar(eivecs(5,[6 7 10]),'FaceColor',[82 105 79]/256)
set(gca,'fontsize',15,'ticklabelinterpreter','latex','xtick',[1:25],'xticklabel',{'ratio 1','ratio 2','$\tau_{inter}$'})
title('PC 5','interpreter','latex')
box on
xtickangle(90)
axis([.5 3.5 -.8 .8])

% same figures with full PCs
figure
bar(eivecs(1,:),'FaceColor',[82 105 79]/256)
set(gca,'fontsize',15,'ticklabelinterpreter','latex','xtick',[1:25],'xticklabel',{'$\lambda_{Planck}$','$\lambda_{fast}^{equil}$','$\lambda_{multidecadal}^{equil}$','$\tau_{fast}$','$\tau_{multidecadal}$','ratio','ratio 2','$\tau_{mixed}$','$\tau_{upper}$','$\tau_{inter}$','$\tau_{deep}$','$\tau_{bottom}$','$I_b$','$a_{CO2}$','$R_{volcanic}$ coeff','$\gamma_{aero-SOx}$','$\gamma_{aero-NOx}$','$\gamma_{aero-OC}$','$\gamma_{aero-BC}$','$\gamma_{aero-NH3}$','AeroRF-ind-2011','$\gamma_{aero-VOC}$','Uncert-CH4','Uncert-N2O','Uncert-Halocarbons'})
ylabel('Loading','interpreter','latex')
title('Principal Component 1','interpreter','latex')
box on
xtickangle(90)

figure
bar(eivecs(2,:),'FaceColor',[82 105 79]/256)
set(gca,'fontsize',15,'ticklabelinterpreter','latex','xtick',[1:25],'xticklabel',{'$\lambda_{Planck}$','$\lambda_{fast}^{equil}$','$\lambda_{multidecadal}^{equil}$','$\tau_{fast}$','$\tau_{multidecadal}$','ratio','ratio 2','$\tau_{mixed}$','$\tau_{upper}$','$\tau_{inter}$','$\tau_{deep}$','$\tau_{bottom}$','$I_b$','$a_{CO2}$','$R_{volcanic}$ coeff','$\gamma_{aero-SOx}$','$\gamma_{aero-NOx}$','$\gamma_{aero-OC}$','$\gamma_{aero-BC}$','$\gamma_{aero-NH3}$','AeroRF-ind-2011','$\gamma_{aero-VOC}$','Uncert-CH4','Uncert-N2O','Uncert-Halocarbons'})
ylabel('Loading','interpreter','latex')
title('Principal Component 2','interpreter','latex')
box on
xtickangle(90)

figure
bar(eivecs(3,:),'FaceColor',[82 105 79]/256)
set(gca,'fontsize',15,'ticklabelinterpreter','latex','xtick',[1:25],'xticklabel',{'$\lambda_{Planck}$','$\lambda_{fast}^{equil}$','$\lambda_{multidecadal}^{equil}$','$\tau_{fast}$','$\tau_{multidecadal}$','ratio','ratio 2','$\tau_{mixed}$','$\tau_{upper}$','$\tau_{inter}$','$\tau_{deep}$','$\tau_{bottom}$','$I_b$','$a_{CO2}$','$R_{volcanic}$ coeff','$\gamma_{aero-SOx}$','$\gamma_{aero-NOx}$','$\gamma_{aero-OC}$','$\gamma_{aero-BC}$','$\gamma_{aero-NH3}$','AeroRF-ind-2011','$\gamma_{aero-VOC}$','Uncert-CH4','Uncert-N2O','Uncert-Halocarbons'})
ylabel('Loading','interpreter','latex')
title('Principal Component 3','interpreter','latex')
box on
xtickangle(90)

figure
bar(eivecs(4,:),'FaceColor',[82 105 79]/256)
set(gca,'fontsize',15,'ticklabelinterpreter','latex','xtick',[1:25],'xticklabel',{'$\lambda_{Planck}$','$\lambda_{fast}^{equil}$','$\lambda_{multidecadal}^{equil}$','$\tau_{fast}$','$\tau_{multidecadal}$','ratio','ratio 2','$\tau_{mixed}$','$\tau_{upper}$','$\tau_{inter}$','$\tau_{deep}$','$\tau_{bottom}$','$I_b$','$a_{CO2}$','$R_{volcanic}$ coeff','$\gamma_{aero-SOx}$','$\gamma_{aero-NOx}$','$\gamma_{aero-OC}$','$\gamma_{aero-BC}$','$\gamma_{aero-NH3}$','AeroRF-ind-2011','$\gamma_{aero-VOC}$','Uncert-CH4','Uncert-N2O','Uncert-Halocarbons'})
ylabel('Loading','interpreter','latex')
title('Principal Component 4','interpreter','latex')
box on
xtickangle(90)

figure
bar(eivecs(5,:),'FaceColor',[82 105 79]/256)
set(gca,'fontsize',15,'ticklabelinterpreter','latex','xtick',[1:25],'xticklabel',{'$\lambda_{Planck}$','$\lambda_{fast}^{equil}$','$\lambda_{multidecadal}^{equil}$','$\tau_{fast}$','$\tau_{multidecadal}$','ratio','ratio 2','$\tau_{mixed}$','$\tau_{upper}$','$\tau_{inter}$','$\tau_{deep}$','$\tau_{bottom}$','$I_b$','$a_{CO2}$','$R_{volcanic}$ coeff','$\gamma_{aero-SOx}$','$\gamma_{aero-NOx}$','$\gamma_{aero-OC}$','$\gamma_{aero-BC}$','$\gamma_{aero-NH3}$','AeroRF-ind-2011','$\gamma_{aero-VOC}$','Uncert-CH4','Uncert-N2O','Uncert-Halocarbons'})
ylabel('Loading','interpreter','latex')
title('Principal Component 5','interpreter','latex')
box on
xtickangle(90)
%}

%%

% load model outputs
Co = readtable('HadCRUT5_Cheng_outputs.csv');
Co = table2array(Co(:,2:end-1));
No = readtable('HadCRUT5_NODC_outputs.csv');
No = table2array(No(:,2:end-1));
Io = readtable('HadCRUT5_noinfill_Cheng_outputs.csv');
Io = table2array(Io(:,3:end-2));

%%

% plot  correlation matrix
wcio = weightedcorrs([C N I; Co No Io]',w);
wcio(wcio==1) = 0;
imagesc(wcio);
colormap(flipud(bluewhitered))

%%

% perform stepwise regressions against PCs
Po = [Co Io No]';
slm_20 = stepwiselm(score(:,1:5),Po(:,1),'Criterion','bic','Weights',w);
slm_100 = stepwiselm(score(:,1:5),Po(:,3),'Criterion','bic','Weights',w);
slm_TCR = stepwiselm(score(:,1:5),Po(:,5),'Criterion','bic','Weights',w);

% perform stepwise regressions against variables
i_slm_20 = stepwiselm(P,Po(:,1),'Criterion','bic');
i_slm_100 = stepwiselm(P,Po(:,3),'Criterion','bic');
i_slm_TCR = stepwiselm(P,Po(:,5),'Criterion','bic');
