start = [1 1];
count = [100 100];
angleArray = h5read("D:\JingALing\Wormpose\Long term tracking data\N2_on_food\N2 on food L_2010_04_20__09_52_38___4___1\processed_results.h5",'/raw/theta',start, count);
numEigWorms = 4;
verbose = 1;

%function [eigenWorms, eigenVals] = findEigenWorms(angleArray, numEigWorms, verbose)
% check for NaN input data
if isempty(~isnan(angleArray))
    error('angleArray contains only NaN values.')
end

% only use non-NaN frames for calculating the covariance matrix
angleArray = angleArray(:,~all(isnan(angleArray)));
covarianceMat = cov(angleArray);

%get the eigenvectors and eigenvalues of the covaraince matrix
[M, eVals] = eig(covarianceMat);

%sort the eigenvalues
eVals = sort(diag(eVals),'descend');

%keep the numEigWorms dimensions that capture most of the variance 
eigenWorms = M(:,end:-1:end - numEigWorms + 1)';
eigenVals = eVals(1:numEigWorms);
projectedAmps = angleArray*(eigenWorms(1:numEigWorms,:)');

if verbose 
    % plot eigenvalues to show fraction of variance captured
    figure
    plot(cumsum(eVals/sum(eVals)),'o','markeredgecolor',[1 0.5 0.1],...
        'markerfacecolor', [1 0.5 0.1],'markersize',8)
    %adjust font and font size
    set(gca, 'FontName', 'Helvetica', 'FontSize', 16);

    % plot the eigenworms
    figure
    for i = 1:numEigWorms
       subplot(ceil(numEigWorms/2),2,i)
        plot(eigenWorms(i,:), 'Color', [1 0.5 0.1], 'LineWidth', 2)
        xlim([0 size(eigenWorms, 2) + 1])
        title(num2str(i))
        %adjust font and font size
        set(gca, 'FontName', 'Helvetica', 'FontSize', 13);
    end
    
    %plot the covariance matrix
    figure
    imagesc(covarianceMat)
    set(gca,'YDir','normal')
    %adjust font and font size
    set(gca, 'FontName', 'Helvetica', 'FontSize', 16);
end
