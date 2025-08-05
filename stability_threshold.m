% Ke Feng, 2025
% Stability threshold of PHWG w. i.i.d. Exponential radii.
% https://arxiv.org/abs/2501.10712

clear all
W = 2; % window [-W,W]^2
num_block = 1e5;
block_random_size = zeros(1,num_block);
sigma_hat = zeros(1,num_block); % service time of each block
RACS_R_mean = 1; % Mean RACS Radius
f = 1; % file length

display('probability of max exclusion')
p_W = exp(-sqrt(2)*W/RACS_R_mean)
dist_link = 0; % pairwise distance - need to review if holds for general value
assert(dist_link<=W);

signal_strength = path_loss_func(dist_link);
noise_strength = 0.05;
display('max rate')
log2(1+signal_strength/noise_strength)

max_block_size = 1e5; % max num of arrivals in a block to cap the complexity 
for kk = 1:num_block
    ii = 0;
    rr = exprnd(RACS_R_mean,1);
    RACS_R = zeros(max_block_size,1);

    while rr < sqrt(2)*W
        ii = ii+1;
        RACS_R(ii) = rr; 
        rr = exprnd(RACS_R_mean,1);
    end
    ii = ii+1;
    RACS_R(ii) = rr; 
    assert(ii<max_block_size);

    num_customer = ii; % number of customers in the current block
    block_random_size(kk) = num_customer;
    RACS_R = RACS_R(1:num_customer);
    sz = [1 num_customer];
    
    
    X_1 = unifrnd(0,2*W,sz)'-W;
    X_2 = unifrnd(0,2*W,sz)'-W;
    X = [X_1, X_2];
    % wrap around receivers
    Theta  = unifrnd(0,2*pi,sz)';
    Y = wrap_points_torus([X_1 + dist_link*cos(Theta),X_2 + dist_link*sin(Theta)],W);
    Y_1 = Y(:,1);
    Y_2 = Y(:,2);
    F = exprnd(f,sz); % file size

    ToD = 0; 
    departure_order = zeros(1,num_customer);
    interference = zeros(1,num_customer);
    rate = zeros(1,num_customer);
    
    F_residual  = F;
    status_precendence = ones(1,num_customer); 
    % 1 means waiting, 0 means ground/active, 100 means left.
    count_departure = 0; 
    interference_graph = zeros(num_customer,num_customer); % initialization of the inteference graph
    
    for ii = 2:num_customer
        dist = dist_cal(X(ii,:), X(1:ii-1,:),W); % only w. prior arrivals
        interference_graph(ii,1:ii-1) = (dist <= RACS_R(ii)+RACS_R(1:ii-1));
    end
    DAG = interference_graph;
    while count_departure < num_customer
        for ii = 1:num_customer
            status_precendence(ii) = sum(DAG(ii,1:ii));
        end
        I_ground_nodes = find(status_precendence==0); % index set of the points on the ground
        for vv = I_ground_nodes 
            dist_ground_nodes = dist_cal(Y(vv,:),X(I_ground_nodes,:),W);
            interference(vv) = sum(path_loss_func(dist_ground_nodes)) -  path_loss_func(0);
            rate(vv) = log2(1+signal_strength./(interference(vv)+noise_strength));
        end
        [delta_t,I_depart_ground] = min(F_residual(I_ground_nodes)./rate(I_ground_nodes));
        I_depart = I_ground_nodes(I_depart_ground);
        DAG(I_depart,:) = 0;
        DAG(:,I_depart) = 0;
        rate(I_depart) = 0;
        DAG(I_depart,I_depart) = 100;
        ToD = ToD + delta_t;
        F_residual = F_residual - rate*delta_t;
        count_departure = count_departure + 1;
        departure_order(I_depart) = count_departure;
    end
sigma_hat(kk) = ToD;
end


display('lambda_c')
1/(p_W*mean(sigma_hat)*4*W^2)


