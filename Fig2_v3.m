%% Some variables

rng=32;
gn=0;% gain noise
n=10;% core n (sample size)
nn=[n n*5 n*30 10*n];%(sample size per column
% nn=[n n^2 n^3];%(sample size per column
N=1e4;% population size
n_ticks=101;% ticks for visualazing population distribution
n_cases=3;

lw=3;
my_color=[37,52,148]/255;
my_color=[0 0 1];
my_color1=[186,228,179
    116,196,118
    35,139,69]/255;
yl1=zeros(n_cases+1,2);

my_color1=[65,182,196;65,182,196;65,182,196]/255;
% my_color1=my_color1*0;
% my_color1(:,end)=1;
ms=20;
mt='.';


mu(:,1)=[100 50 1];
covar=ones(3,3,1)*.1;
covar(1:2,1:2,1)=[10 -.7 ;-.7 30];
% covar(1:2,1:2,1)=[10 3 ;3 30];
covar(2,3)=21/3;
covar(3,2)=21/3;
covar(3,3)=10;
data=mvnrnd(mu, covar,N);

pph=.1;%percent plot histogram
%% Fontsizes
set(0,'DefaultTextFontSize',16);
fs_label=14;
fs_title=16;
%% titles
tit1{1}='Pop. dist';
for i=1:n_cases
    tit1{i+1}=[num2str(nn(i)) ' samples'];
end

%% generate rando distributions
X=data;
Xn=data+randn(size(X))*gn;

%% fit 1d data
ftemp=figure;
hist(X(:,1),n_ticks);
axis tight
xl=xlim;
hist(X(:,2),n_ticks);
axis tight
yl2=xlim;
close(ftemp)
pd = fitdist(X(:,1),'Normal');
x1d=linspace(xl(1),xl(2),n_ticks);
y1d = pdf(pd,x1d);

X1d_sample=cell(n_cases,1);
Y1d_sample=cell(n_cases,1);
y1d_fit=zeros(n_cases,n_ticks);

for i=1:n_cases
    ix=randperm(N);
    X1d_sample{i}=Xn(ix(1:nn(i)),:);
    
    pd = fitdist(X1d_sample{i}(:,1),'Normal');
    Y1d_sample{i} = pdf(pd,X1d_sample{i});
    y1d_fit(i,:) = pdf(pd,x1d);
end

%% Mess with 2d data
mu2d = mu(1:2)';
Sigma = covar(1:2,1:2);


x1 = linspace(xl(1),xl(2),n_ticks);
x2 = linspace(yl2(1),yl2(2),n_ticks);
% x1 = -3:.2:3; x2 = -3:.2:3;

[X1,X2] = meshgrid(x1,x2);
F = mvnpdf([X1(:) X2(:)],mu2d,Sigma);
F = reshape(F,length(x2),length(x1));
mvncdf([0 0],[1 1],mu2d,Sigma);
%% Mess with 3d data
FA=0.1;
EC='w';
%[x, y, z] = ellipsoid(mu(1),mu(2),mu(3),covar(1,1),covar(2,2),covar(3,3),17);
factor=2.5;
[x, y, z] = ellipsoid(mu(1),mu(2),mu(3),factor*sqrt(covar(1,1)),factor*sqrt(covar(2,2)),factor*sqrt(covar(3,3)),11);


%% Define plot settings
fig_size=[1 1 24 16];
fig_name='Dim';
afs=9;% axis font size
f = figure('Units','centimeters',...
    'PaperUnits','centimeters',...
    'name',fig_name,...
    'PaperPosition',fig_size,...
    'Position',fig_size,...
    'DefaultAxesFontSize',afs,...
    'color',[1 1 1]);

marg_plot_perc=10/100;

s1{1}=subplot (3,6,1);
s1{2}=subplot (3,6,2);
s1{3}=subplot (3,6,3);
s1{4}=subplot (3,6,4);
s1{5}=subplot (3,6,5);
s1{6}=subplot (3,6,6);
s2{1}=subplot (3,6,7);
s2{2}=subplot (3,6,8);
s2{3}=subplot (3,6,9);
s2{4}=subplot (3,6,10);
s2{5}=subplot (3,6,11);
s2{6}=subplot (3,6,12);
for i=1:6
    s3{i}=subplot (3,6,12+i);
end
%%
% clf
marg_plot_perc=10/100;

% s1{1}=subplot (3,5,1);
% s1{2}=subplot (3,5,2);
% s1{3}=subplot (3,5,3);
% s1{4}=subplot (3,5,4);
% s1{5}=subplot (3,5,5);
% s2{1}=subplot (3,5,6);
% s2{2}=subplot (3,5,7);
% s2{3}=subplot (3,5,8);
% s2{4}=subplot (3,5,9);
% s2{5}=subplot (3,5,10);
% for i=1:4
%     s1{i}=subplot (3,4,i);
%     s2{i}=subplot (3,4,4+i);
%     s3{i}=subplot (3,4,8+i);
% end


% for i=1:4
%     subplot(s2{i});
%     pos=get(gca,'position');
%     s2x{i}=subplot('position',[pos(1)* pos(2) pos(3) pos(4)*marg_plot_perc]);
%     s2y{i}=subplot('position',[pos(1) pos(2) pos(3)*marg_plot_perc pos(4)*marg_plot_perc])
%
%
% end

%%
subplot (s1{1})
set(gca,'FontSize',fs_label);
plot(x1d,y1d,...
    'color',my_color,...
    'linewidth',lw)
yl1(1,:)=ylim;

for i=1:n_cases
    subplot (s1{i+1})
    set(gca,'FontSize',fs_label);    
    plot(x1d,y1d,...
        'color',my_color,...
        'linewidth',lw)
    hold all
    
    %plot(x1d,y1d_fit(i,:),...
     %   'color',my_color1(i,:),...
      %  'linewidth',lw)
    %plot(x1d,y1d_fit(i,:),...
     %   'color','w',...
     %   'linewidth',lw-2)
    
    plot(X1d_sample{i},Y1d_sample{i},'.',...
        'marker',mt,...
        'markersize',ms,...
        'color',my_color1(i,:))
    plot(X1d_sample{i},Y1d_sample{i},'.',...
        'marker',mt,...
        'markersize',ms/3,...
        'color','w')
    
    hold off
    yl1(i+1,:)=ylim;
end
yl1=max(yl1);
for i=1:n_cases+1
    subplot (s1{i})
    set(gca,'FontSize',fs_label);    
    axis square
    ylim(yl1)
    xlim(xl);
    set(gca,'ytick',[])
    if i==1
        ylabel('Relative abundance',...
            'fontsize',fs_label);
    end
    xlabel('Trait 1',...
        'fontsize',fs_label);
    title(tit1{i},...
        'fontsize',fs_title)
    
end

%%
subplot (s2{1})
set(gca,'FontSize',fs_label);
contour(x1,x2,F,n_ticks);
yl1(1,:)=ylim;
xl1(1,:)=xlim;
hold all
%plot(x1d,pph*diff(ylim)*y1d/max(y1d)+min(ylim),...
%        'color',my_color,...
%        'linewidth',lw)
    
    hold off


for i=1:n_cases
    subplot (s2{i+1})
    set(gca,'FontSize',fs_label);
    contour(x1,x2,F,n_ticks);
    hold all
    x_temp = X1d_sample{i}(:,1);
    y_temp = X1d_sample{i}(:,2);
    xpts = linspace(90,110,20);
    ypts = linspace(30,70,40);
    empty_mat = nan(70,110);
    [xG,yG] = meshgrid(-5:5);
    sigma = 2.5;
    g = exp(-xG.^2./(2.*sigma.^2)-yG.^2./(2.*sigma.^2));
    g = g./sum(g(:));
    conv_mat = conv2(histcounts2(x_temp,y_temp,xpts,ypts)',g,'same');
    empty_mat(31:69,91:109) = conv2(histcounts2(x_temp,y_temp,xpts,ypts)',g,'same');
    h_hist = imagesc(empty_mat);
    cmap_data = colormap(jet);
    cmap_data(1,:) = [1 1 1];
    colormap(cmap_data)
%    set(h_hist,'AlphaData',0.7);
    %plot(X1d_sample{i}(:,1),X1d_sample{i}(:,2),'.',...
    %    'marker',mt,...
    %    'markersize',ms,...
    %    'color',my_color1(i,:))
    %plot(X1d_sample{i}(:,1),X1d_sample{i}(:,2),'.',...
    %    'marker',mt,...
    %    'markersize',ms/3,...
    %    'color','w')
    
       
    
    
    yl1(i+1,:)=ylim;
    xl1(i+1,:)=ylim;
%    plot(x1d,pph*diff(ylim)*y1d/max(y1d)+min(yl1(i+1,:)),...
%        'color',my_color,...
%        'linewidth',lw)
%    plot(x1d,pph*diff(ylim)*y1d_fit(i,:)/max(y1d)+min(yl1(i+1,:)),...
%        'color',my_color1(i,:),...
%        'linewidth',lw)
%    plot(x1d,pph*diff(ylim)*y1d_fit(i,:)/max(y1d)+min(yl1(i+1,:)),...
%        'color','w',...
%        'linewidth',lw-2)
    
    hold off
    
end
yl1=max(yl1);
xl1=max(xl1);

for i=1:n_cases+1
    subplot (s2{i})
    set(gca,'FontSize',fs_label);
    axis square
    ylim(yl1)
    xlim(xl);
%     set(gca,'ytick',[])
    if i==1
        ylabel('Trait 2',...
            'fontsize',fs_label);
    end
    xlabel('Trait 1',...
        'fontsize',fs_label);
    xt=get(gca,'xtick');
    yt=get(gca,'ytick');
    %     title(tit1{i},...
    %         'fontsize',fs_title)
    
end
%%
subplot (s3{1})
set(gca,'FontSize',fs_label);
%h=surf(x, y, z);
%set(h,'FaceColor',my_color)
%set(h,'EdgeColor',EC)
%set(h,'FaceAlpha',FA)
%grid off
 yl1(1,:)=ylim;
xl1(1,:)=xlim;
zl1(1,:)=zlim;

for i=1:n_cases
    subplot (s3{i+1})
    set(gca,'FontSize',fs_label);    
%    h=surf(x, y, z);
%    set(h,'FaceColor',my_color)
%    set(h,'EdgeColor',EC)
%    set(h,'FaceAlpha',FA)
%    grid off
    
%    hold all
    
%    plot3(X1d_sample{i}(:,1),X1d_sample{i}(:,2),X1d_sample{i}(:,3),'.',...
%        'marker',mt,...
%        'markersize',ms,...
%        'color',my_color1(i,:))
%    plot3(X1d_sample{i}(:,1),X1d_sample{i}(:,2),X1d_sample{i}(:,3),'.',...
%        'marker',mt,...
%        'markersize',ms/3,...
%        'color','w')
    
%         plot3(x1d,-10*y1d/max(y1d)+yl1(1),zeros(size(x1d)),'color',my_color,'linewidth',lw)

%    hold off
%     yl1(i+1,:)=ylim;
%    xl1(i+1,:)=ylim;
%    zl1(i+1,:)=zlim;
end

%yl1=max(yl1);
%xl1=max(xl1);
%zl1=max(zl1);

for i=1:n_cases+1
    subplot (s3{i})
    set(gca,'FontSize',fs_label);    
%    axis square
%    ylim(yl1)
%    xlim(xl1);
%    zlim(zl1);
%    axis  tight
%    xlim(xl1);
    
%    if i==1
        
%        zlabel('Trait 3',...
%            'fontsize',fs_label);
%        ylabel('Trait 2',...
%        'fontsize',fs_label);
%    end
    
%    xlabel('Trait 1',...
%        'fontsize',fs_label);
%    set(gca,'xtick',xt);
%    set(gca,'ytick',yt);
%    view(-15,45)
    %     title(tit1{i},...
    %         'fontsize',fs_title)
    
end
%colormap jet

%%
s5{1}=s1{5};
s5{2}=s2{5};
s5{3}=s3{5};
s6{1}=s3{1};
s6{2}=s3{2};
s6{3}=s3{3};
%%
load('boxplot_data.mat');
to_plot_flag=2;
to_plot{1}='MSE';
to_plot{2}='MAPE';
to_plot{3}='Accuracy';

switch to_plot_flag
    case 1
        Y=mse_params;
    case 2
        Y=mape_params;
    case 3
        Y=accuracy;
end
%
delta=5;
clear yl1
for i=1:3
    subplot(s5{i})
    set(gca,'FontSize',fs_label);
    
    my=nanmean(Y(:,:,i),2);
    my_end(i)=my(end);
    plot(n_fit,my,...
        'color',my_color1(i,:),...
        'linewidth',lw)
    hold all
    
    
    
    ptiles=prctile(Y(:,:,i)',[delta 100-delta]);
    ix_patch=sum(ptiles);
%     ix_patch=
%     ix_patch=find(~isnan(ix_patch))
    xpatch=[n_fit n_fit(end:-1:1)];
    ypatch=[my' ptiles(1,end:-1:1)];
    cpatch=ones(size(ypatch));
    cpatch(end/2+1:end)=0;
    patch(xpatch,ypatch,cpatch)
    
    ypatch=[my' ptiles(2,end:-1:1)];
    cpatch=ones(size(ypatch));
    cpatch(end/2+1:end)=0;
    patch(xpatch,ypatch,cpatch)
    plot(n_fit,my,...
        'color',my_color1(i,:),...
        'linewidth',lw)
     plot(n_fit,my,...
        'color','w',...
        'linewidth',lw-2)
    
    
    yl1(i,:)=ylim;
    hold off
end

yl(1)=min(yl1(:,1));
yl(2)=max(yl1(:,2));

for i=1:3
    subplot(s5{i})
    set(gca,'FontSize',fs_label);
    axis square
    colormap(s5{i},winter)
    ylim([0 1])
    currpos=get(gca,'Position');
    set(gca,'Position',[currpos(1)-0.002,currpos(2),currpos(3),currpos(4)])
    set(gca,'YTick',[0, 0.5 1])
    hold
    plot(1:1000,repmat(0.5,1000,1),'k-','LineWidth',3);
    hold    
    set(gca,'yticklabel',num2str(get(gca,'ytick')','%4.1f'))
    if i==1
        title(to_plot{to_plot_flag},...
        'fontsize',fs_title)
    end
    xlabel('Sample size',...
        'fontsize',fs_label);
    text(n_fit(end),my_end(i),num2str(my_end(i),'%4.2f'),...
            'fontsize',fs_label-2)
    ylabel('MAE parameter estimates',...
        'fontsize',fs_label)
end
to_plot_flag = 3;
switch to_plot_flag
    case 1
        Y=mse_params;
    case 2
        Y=mape_params;
    case 3
        Y=accuracy;
end
%
delta=5;
clear yl1
for i=1:3
    subplot(s3{i+1})
    set(gca,'FontSize',fs_label);
    my=nanmean(Y(:,:,i),2);
    my_end(i)=my(end);
    plot(n_fit,my,...
        'color',my_color1(i,:),...
        'linewidth',lw)
    hold all
    
    
    
    ptiles=prctile(Y(:,:,i)',[delta 100-delta]);
    ix_patch=sum(ptiles);
%     ix_patch=
%     ix_patch=find(~isnan(ix_patch))
    xpatch=[n_fit n_fit(end:-1:1)];
    ypatch=[my' ptiles(1,end:-1:1)];
    cpatch=ones(size(ypatch));
    cpatch(end/2+1:end)=0;
    patch(xpatch,ypatch,cpatch)
    
    ypatch=[my' ptiles(2,end:-1:1)];
    cpatch=ones(size(ypatch));
    cpatch(end/2+1:end)=0;
    patch(xpatch,ypatch,cpatch)
    plot(n_fit,my,...
        'color',my_color1(i,:),...
        'linewidth',lw)
     plot(n_fit,my,...
        'color','w',...
        'linewidth',lw-2)
    
    
    yl1(i,:)=ylim;
    hold off
end

yl(1)=min(yl1(:,1));
yl(2)=max(yl1(:,2));
posshift=[0,0,0];
for i=1:3
    subplot(s3{i+1})
    set(gca,'FontSize',fs_label);
    axis square
    colormap(s3{i+1},winter)
    ylim(yl)
%    currpos=get(gca,'Position');
%    set(gca,'Position',[currpos(1)+posshift(i),currpos(2),currpos(3),currpos(4)])    
    hold
    plot(1:1000,repmat(0.5,1000,1),'k-','LineWidth',3);
    hold
    set(gca,'YTick',[])
    if i==1
        ylabel('out of sample acc.',...
        'FontSize',fs_label);
        set(gca,'YTick',[0,0.5,1]) 
        set(gca,'yticklabel',num2str(get(gca,'ytick')','%4.1f'))
        title({'1 dimension'},...
        'fontsize',fs_title)
    elseif i==2
        xlabel('Sample size',...
        'fontsize',fs_label);
    end
    text(n_fit(end),my_end(i),num2str(my_end(i),'%4.2f'),...
            'fontsize',fs_label-2)
    title({[num2str(i) ' ' 'dimensions']},...
        'fontsize',fs_title)
end
delete(subplot(3,6,[6 12 18]))
delete(subplot(3,6,[5 11 17]))
delete(subplot(3,6,13))