addpath('C:/Users/feczko/Box Sync/Heterogeneity_overview_paper/figures_in_slides/MVoutlierCode/')
number_dims=3;
nreps = 1000;
sub_range = [10 50 100 200 300 400 500 600 700 800 900 1000];
accuracy = zeros(3,length(sub_range),nreps,number_dims);
for curr_dims = 1:number_dims
    nsub_count = 0;
    if curr_dims > 1
        temp_data = mvnrnd(zeros(curr_dims,1),ones(curr_dims,1)',10000);
        [mu,s,RD,chi_crt] = DetectMultVarOutliers(temp_data(1:end-1,:),[],[],false);
        temp_outliers = find(RD > chi_crt(2))';
        temp_notoutliers = find(RD < chi_crt(2))';
        temp_data(end+1,:) = zeros(size(temp_data,2),1);
        euc_mat = squareform(pdist(temp_data));
        outlier_thresh = min(euc_mat(temp_outliers,end));
    else
        outlier_thresh = 2;
    end    
    for nsubs = sub_range
        nsub_count = nsub_count + 1;
        for iter = 1:nreps
            temp_data = mvnrnd(zeros(curr_dims,1),ones(curr_dims,1)',nsubs);
            temp_data(end+1,:) = zeros(size(temp_data,2),1);
            euc_mat = squareform(pdist(temp_data));
            true_outliers = find(euc_mat(end,1:end-1) > outlier_thresh);
            true_nonoutliers = find(euc_mat(end,1:end-1) < outlier_thresh);
            if curr_dims > 1
                [mu,s,RD,chi_crt] = DetectMultVarOutliers(temp_data(1:end-1,:),[],[],false);
                observed_outliers = find(RD > chi_crt(2))';
                observed_notoutliers = find(RD < chi_crt(2))';
            else
                observed_outliers = find(abs(temp_data(1:end-1)) > 2*std(temp_data(1:end-1)));
                observed_notoutliers = find(abs(temp_data(1:end-1)) < 2*std(temp_data(1:end-1)));
            end
            accuracy(2,nsub_count,iter,curr_dims) = length(intersect(true_outliers,observed_outliers))/length(true_outliers);
            accuracy(3,nsub_count,iter,curr_dims) = length(intersect(true_nonoutliers,observed_notoutliers))/length(true_nonoutliers);
            accuracy(1,nsub_count,iter,curr_dims) = (length(true_nonoutliers)/nsubs)*accuracy(3,nsub_count,iter) + (length(true_outliers)/nsubs)*accuracy(2,nsub_count,iter);
        end
    end
end
accuracy1D = accuracy(:,:,:,1);
accuracy2D = accuracy(:,:,:,2);
accuracy3D = accuracy(:,:,:,3);
%reshape outliers for good boxplotting
outlier1D(:,:,1) = accuracy1D(2,:,:);
outlier2D(:,:,2) = accuracy2D(2,:,:);
outlier3D(:,:,3) = accuracy3D(2,:,:);
%%create a simple figure -- canceled for now
%plot(mean(accuracy1D(2,:,:),3,'omitnan'),'blue','linewidth',5)
%hold
%plot(mean(accuracy2D(2,:,:),3,'omitnan'),'red','linewidth',5)
%plot(mean(accuracy3D(2,:,:),3,'omitnan'),'black','linewidth',5)
%title('mean outlier detection accuracy by sample size')
%ylabel('percent correctly identified')
%xlabel('sample size')
%set(gca,'XTick',1:length(sub_range))
%set(gca,'XTickLabel',sub_range)
%legend({'1 dimension','2 dimensions','3 dimensions'})

%%oscar's boxplot code modified for this plot
n_fit=sub_range;
n_fit=sort(n_fit,'ascend');
n_n_fit=length(n_fit);
for i=1:3
    
    subplot (3,1,i)
    boxplot(squeeze(accuracy(:,:,i)'),n_fit)
    if i==1
        title({'outlier identification performance:',[num2str(nreps) ' repetitions']})
    end
    if i==3
        xlabel('Sample size')
    end
    ylim([0 1])
    ylabel({'outlier accuracy (%) for',[num2str(i) ' dimensions']})
end