clear
close 'all'
%the files produced by the main code, names are set to match default
%outputs of that code
load('means.mat')%the file for the means
load('vars.mat')%variances
load('skews.mat')%3rd moments
load('kurts.mat')%4th moments

%adjustment to the usual normalised form of skewness and kurtosis
skews=skews./vars.^(3/2);
kurts=kurts./vars.^2;

figure;
hold on
h=histogram(means(1,:),'BinWidth',0.02,'Normalization','probability');
histogram(means(2,:),'BinWidth',0.02,'Normalization','probability');
hold off
title('Means')
legend('\phi_1','\phi_2')
figure;
hold on
histogram(vars(1,:),'BinWidth',0.03,'Normalization','probability');
histogram(vars(2,:),'BinWidth',0.03,'Normalization','probability');
hold off
title('Variances')
legend('\phi_1','\phi_2')
figure;
hold on
histogram(skews(1,:),'BinWidth',0.3,'Normalization','probability');
histogram(skews(2,:),'BinWidth',0.3,'Normalization','probability');
hold off
title('Skews')
legend('\phi_1','\phi_2')
figure;
hold on
h=histogram(kurts(1,:),'BinWidth',3,'Normalization','probability');
histogram(kurts(2,:),'BinWidth',3,'Normalization','probability');
hold off
title('Kurtosises')
legend('\phi_1','\phi_2')

disp("Median of phi means:")
disp(["phi_1:", median(means(1,:)),"phi_2:", median(means(2,:))])
disp("Median of phi variances:")
disp(["phi_1:", median(vars(1,:)),"phi_2:", median(vars(2,:))])
disp("Median of phi means:")
disp(["phi_1:", median(skews(1,:)),"phi_2:", median(skews(2,:))])
disp("Median of phi means:")
disp(["phi_1:", median(kurts(1,:)),"phi_2:", median(kurts(2,:))])
