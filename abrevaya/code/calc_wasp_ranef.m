% $$$ keys=c(5, 6, 7, 9, 10, 13) 
% $$$ [1] "covRanef[1,1]" "covRanef[2,1]" "covRanef[3,1]" "covRanef[2,2]"
% $$$ [5] "covRanef[3,2]" "covRanef[3,3]"
function calc_wasp_ranef(dd, nsub, ndim, keys, dirp)

addpath('/opt/gurobi/6.0.4/linux64/matlab/')

covs = {};

for jj = 1:nsub
    tmp = csvread(strcat(dirp, 'data/wasp_samp/', num2str(dd), '/samp_', num2str(jj), '.csv'));        
    for ii = 1:ndim
        covs{jj, ii} = tmp(:, keys(ii)); 
    end
end

rtime = zeros(ndim, 1);

% calculate the pair-wise sq. euclidean distance between the atoms of subset
% posteriors and WASP atoms
for dims = 1:ndim
    subsetPost = {};    
    for jj = 1:nsub
        subsetPost{jj} = covs{jj, dims}(randi([1 1000], 200, 1));
    end         

    lbd = min(cellfun(@(x) min(x), subsetPost));
    ubd = max(cellfun(@(x) max(x), subsetPost));    
    overallPost = linspace(lbd, ubd, 1000)';

    distMatCell = {};
    
    m00 = diag(overallPost * overallPost');
    for ii = 1:nsub
        mm = diag(subsetPost{ii} * subsetPost{ii}');    
        mm1 = overallPost * subsetPost{ii}'; 
        distMatCell{ii} = bsxfun(@plus, bsxfun(@plus, -2 * mm1, mm'), m00);    
    end
    
    % constants
    K  = nsub;
    Ni = cell2mat(cellfun(@(x) size(x, 2), distMatCell, 'UniformOutput', false));
    N  = size(overallPost, 1);
    nx = N * (N+1);
    mx = K * N + N + 1;
    In = eye(N);
    En = ones(1, N);

    % Generate matrix A0.
    A0  = sparse([]);
    for p = 1:K
        cc = (1:N)';                  % terribly fast version of 
        idx = cc(:, ones(Ni(p), 1));  % repmat(In, 1, Ni(p)) / Ni(p)
        Rp  = In(:, idx(:)) / Ni(p);  % in 3 steps
        A0  = blkdiag(A0, Rp); 
    end
    cc = (1:N)';                  % terribly fast version of 
    idx = cc(:, ones(K, 1));      % repmat(-In, K, 1) 
    A00  = -In(idx(:), :);        % in 3 steps
    
    A0 = sparse([A00, A0]);
    b0 = zeros(size(A0, 1), 1);
    disp('done generating A ...');        
    
    % Generate matrix B from simplex constraints.
    B = sparse([]);
    for p = 0:(sum(Ni))
        B = blkdiag(B, En);
    end
    disp('done generating B ...');        
    
    % The hold matrix C.
    A = sparse([A0; B]);

    % Generate the right hand size vector b.
    b = sparse([zeros(K * N, 1); ones(sum(Ni) + 1, 1)]);
    
    % Generate the cost vector
    costCell = cellfun(@(x) x(:) / size(x, 2), distMatCell, 'UniformOutput', false);
    costVec = [zeros(size(overallPost, 1), 1); cell2mat(costCell(:))];
    
    c = sparse(costVec);
    
    tic;
    lpsol = callLpSolver('gurobi', A, b, c, 10000, 1e-10);
    rtime(dims, 1) = toc;
    
    [tmats, avec] = recoverSolution(lpsol, K, N, Ni);

    save(strcat(dirp, 'result/wasp/ranef_fit_', num2str(dims), '_rep_', num2str(dd), '.mat'), 'tmats','avec', 'subsetPost', 'overallPost');    
    summ = [overallPost avec];
    csvwrite(strcat(dirp, 'result/wasp/ranef_', num2str(dims), '_rep_', num2str(dd),  '.csv'), summ);

    disp(['evaluating sim ' num2str(dd) '...' ' done with dim ' num2str(dims) '... ']);    
end
csvwrite(strcat(dirp, 'result/wasp/ranef_times_rep_', num2str(dd),  '.csv'), rtime);    

quit
