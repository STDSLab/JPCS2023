%% Load data

data = importdata('../Data/Matlab/Oregon_ADHD_data_exg.csv')
data_exg = zscore(data.data);

data = importdata('../Data/Matlab/Oregon_ADHD_data_dti.csv')
data_dti = zscore(data.data);

data = importdata('../Data/Matlab/Oregon_ADHD_data_area.csv')
data_area = zscore(data.data);

data = importdata('../Data/Matlab/Oregon_ADHD_data_thick.csv')
data_thick = zscore(data.data);

data = importdata('../Data/Matlab/Oregon_ADHD_data_vol.csv')
data_vol = zscore(data.data);

data = importdata('../Data/Matlab/Oregon_ADHD_data_sulc.csv')
data_sulc = zscore(data.data);

%% compute networks

corr_exg = (1-pdist(data_exg,'correlation'));
corr_dti = (1-pdist(data_dti,'correlation'));
corr_area = (1-pdist(data_area,'correlation'));
corr_vol = (1-pdist(data_vol,'correlation'));
corr_thick = (1-pdist(data_thick,'correlation'));
corr_sulc = (1-pdist(data_sulc,'correlation'));


net_exg = binary_conn_net(corr_exg);
net_dti = binary_conn_net(corr_dti);
net_area = binary_conn_net(corr_area);
net_vol = binary_conn_net(corr_vol);
net_thick = binary_conn_net(corr_thick);
net_sulc = binary_conn_net(corr_sulc);


%% find communities in real data - uses BCT for community_louvain.m (https://sites.google.com/site/bctnet/)

[M_exg, Q_exg, Q_exg_track] = community_louvain(net_exg);
[M_dti, Q_dti, Q_dti_track] = community_louvain(net_dti);
[M_area, Q_area, Q_area_track] = community_louvain(net_area);
[M_vol, Q_vol, Q_vol_track] = community_louvain(net_vol);
[M_thick, Q_thick, Q_thick_track] = community_louvain(net_thick);
[M_sulc, Q_sulc, Q_sulc_track] = community_louvain(net_sulc);


histc(M_exg, unique(M_exg))
histc(M_dti, unique(M_dti))
histc(M_area, unique(M_area))
histc(M_vol, unique(M_vol))
histc(M_thick, unique(M_thick))


%% Randomize data and find communities

data_exg_rand = data_exg;
data_dti_rand = data_dti;
data_area_rand = data_area;
data_vol_rand = data_vol;
data_thick_rand = data_thick;
data_sulc_rand = data_sulc;

Q_exg_vec = zeros(1000,1);
Q_dti_vec = zeros(1000,1);
Q_area_vec = zeros(1000,1);
Q_vol_vec = zeros(1000,1);
Q_thick_vec = zeros(1000,1);
Q_sulc_vec = zeros(1000,1);


for i=1:1000
    for j=1:size(data_exg,2)
        data_exg_rand(:,j) = data_exg_rand(randperm(size(data_exg_rand,1)),j);
    end

    for j=1:size(data_dti,2)
        data_dti_rand(:,j) = data_dti_rand(randperm(size(data_dti_rand,1)),j);
    end

    for j=1:size(data_area,2)
        data_area_rand(:,j) = data_area_rand(randperm(size(data_area_rand,1)),j);
    end    

    for j=1:size(data_vol,2)
        data_vol_rand(:,j) = data_vol_rand(randperm(size(data_vol_rand,1)),j);
    end       

    for j=1:size(data_thick,2)
        data_thick_rand(:,j) = data_thick_rand(randperm(size(data_thick_rand,1)),j);
    end  

    for j=1:size(data_sulc,2)
        data_sulc_rand(:,j) = data_sulc_rand(randperm(size(data_sulc_rand,1)),j);
    end  


    net_exg_rand = binary_conn_net(1-pdist(data_exg_rand,'correlation'));
    [M_exg_rand, Q_exg_rand] = community_louvain(net_exg_rand);
    Q_exg_vec(i) = Q_exg_rand;

    net_dti_rand = binary_conn_net(1-pdist(data_dti_rand,'correlation'));
    [M_dti_rand, Q_dti_rand] = community_louvain(net_dti_rand);
    Q_dti_vec(i) = Q_dti_rand;

    net_area_rand = binary_conn_net(1-pdist(data_area_rand,'correlation'));
    [M_area_rand, Q_area_rand] = community_louvain(net_area_rand);
    Q_area_vec(i) = Q_area_rand;

    net_vol_rand = binary_conn_net(1-pdist(data_vol_rand,'correlation'));
    [M_vol_rand, Q_vol_rand] = community_louvain(net_vol_rand);
    Q_vol_vec(i) = Q_vol_rand;

    net_thick_rand = binary_conn_net(1-pdist(data_thick_rand,'correlation'));
    [M_thick_rand, Q_thick_rand] = community_louvain(net_thick_rand);
    Q_thick_vec(i) = Q_thick_rand;

    net_sulc_rand = binary_conn_net(1-pdist(data_sulc_rand,'correlation'));
    [M_sulc_rand, Q_sulc_rand] = community_louvain(net_sulc_rand);
    Q_sulc_vec(i) = Q_sulc_rand;

end



%% Plot Q for real and null models

figure;
subplot(3,2,1)
histogram(Q_exg_vec,10)
hold on;
xline(Q_exg,'r','LineWidth',2); title(['RTV (p-val =', num2str(round(sum(Q_exg_vec>=Q_exg)/1000,2)),')']);
xlabel('Modularity');

subplot(3,2,2)
histogram(Q_dti_vec,10)
hold on;
xline(Q_dti,'r','LineWidth',2); title(['DTI (p-val =', num2str(round(sum(Q_dti_vec>=Q_dti)/1000,2)),')']);
xlabel('Modularity');


subplot(3,2,3)
histogram(Q_area_vec,10)
hold on;
xline(Q_area,'r','LineWidth',2); title(['Area (p-val =', num2str(round(sum(Q_area_vec>=Q_area)/1000,2)),')']);
xlabel('Modularity');

subplot(3,2,4)
histogram(Q_vol_vec,10)
hold on;
xline(Q_vol,'r','LineWidth',2); title(['Volume (p-val =', num2str(round(sum(Q_vol_vec>=Q_vol)/1000,2)),')']);
xlabel('Modularity');


subplot(3,2,5)
histogram(Q_thick_vec,10)
hold on;
xline(Q_thick,'r','LineWidth',2); title(['Thickness (p-val =', num2str(round(sum(Q_thick_vec>=Q_thick)/1000,2)),')']);
xlabel('Modularity');


subplot(3,2,6)
histogram(Q_sulc_vec,10)
hold on;
xline(Q_sulc,'r','LineWidth',2); title(['SulcDepth (p-val =', num2str(round(sum(Q_sulc_vec>=Q_sulc)/1000,2)),')']);
xlabel('Modularity');

saveas(gcf,'Oregon_ADHD_data_community.png')

%% function binary_conn_net uses another function connectedComponentMeasures from matlab-networks-toolbox (https://github.com/ivanbrugere/matlab-networks-toolbox/tree/master)

function [net, q] = binary_conn_net(A)
q = 0.99;
nc = connectedComponentMeasures(squareform(A>=q));
while(nc>1)
    q = q - 0.01;
    nc = connectedComponentMeasures(squareform(A>=q));
end

net = squareform(A>=q);

