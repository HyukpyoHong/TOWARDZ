function [] = TOWARDZ(ndim, max_order, edgeoccnum)
max_order_gen = max_order;
tmp_mat =  nchoosek(1:(ndim+max_order_gen), ndim)';

total_complexes = nan(ndim, nchoosek(ndim+max_order_gen, max_order_gen));

num_spec = size(total_complexes, 1);
cmplx_num = size(total_complexes, 2);
edge_num = cmplx_num * (cmplx_num - 1);
% Line 12-17: generating all posible complexes under a given maximum
% order, 'max_order_gen.'
for ci = 1:cmplx_num
    total_complexes(1,ci) =  tmp_mat(1,ci)-1;
    for ni = 2:ndim
        total_complexes(ni,ci) =  tmp_mat(ni,ci) - tmp_mat(ni-1,ci) - 1;
    end
end

tmp_mat = nan(ndim, cmplx_num*(cmplx_num-1)/2);
cid = 0;
for ii = 1:cmplx_num
    for jj = (ii+1):cmplx_num
        cid = cid+1;
        tmp_mat(:, cid) = total_complexes(:, jj) - total_complexes(:, ii);
    end
end

% Line 31: Here the variable 'base_edge' contains all the 
% stoichiometric vectors for all possible complexes.
% However, it contains only either -v or v, not both.
base_edge = unique(tmp_mat', 'rows')';

tic    
rng(2) 

% search_num = 0;
cum_num_cand = 0; % the number of total translated network generated.
cum_num_filtered_by_pos_sol = 0; % the number of total translated network generated.

cum_num_cand_filtered_by_def = 0; % Cumulative number of filtered translated CRNs by the condition s+1 <= n <= 2s (Theorem 3.2)
cum_num_filtered_by_num_comp = 0;  % Cumulative number of CRNs filtered by the conditions in Theorem 3.3 and Theorem 3.5.


cum_num_WR_ZD = 0; % Cumulative number of translated CRNs having WR and ZD
cum_num_notrans_WR = 0; % Cumulative number of CRNs having WR without translation.
cum_num_notrans_ZD = 0; % Cumulative number of CRNs having ZD without translation.
cum_num_notrans_WR_ZD = 0; % Cumulative number of CRNs having WR and ZD without translation.
cum_num_trans_WR = 0; % Cumulative number of CRNs having WR with translation.
cum_num_trans_ZD = 0; % Cumulative number of CRNs having ZD with translation.
cum_num_trans_WR_ZD = 0; % Cumulative number of CRNs having WR and ZD with translation.

% mat_size = 10^7

search_num_max = min(10^7, nchoosek(edge_num, edgeoccnum));

for search_num = 1:search_num_max % except "config_idx = 1" to avoid empty graph.
    
    empty_TF = 0;
    wr_TF = 0;
    dz_TF = 0;
    wrdz_TF = 0;
    
    wr_trans_TF = 0;
    dz_trans_TF = 0;
    wrdz_trans_TF = 0;
    
    
    edge_tmp_idx = 1;
    adj_matrix = zeros(cmplx_num, cmplx_num);
    sources = [];
    products = [];
    
    edgeTF = zeros(2, edge_num);
    % edgeTF = edgeTF_tmp < edge_prob;
    % edgeTF(2,:) = 0;
    
    % randomly generate a matrix with a given number of reactions,
    % 'edgeoccnum'
    edgeTF(1,randperm(edge_num, edgeoccnum)) = 1;
    
    % This config_idx identifies the graph configuration. It can be used to
    % check whether two results are from the same CRNs or not.
    config_idx = 0;
    for ee = 1:edge_num
        config_idx = config_idx + edgeTF(1,ee) * 2^(ee-1);
    end
    
    % For a given configuration, generate the source complexes matrix and
    % the product complex matrix.
    for ii = 1:cmplx_num
        for jj = 1:cmplx_num
            if ii ~= jj
                adj_matrix(ii,jj) =  edgeTF(1, edge_tmp_idx);
                edge_tmp_idx = edge_tmp_idx + 1;
            end
            if adj_matrix(ii,jj)
                sources = [sources, total_complexes(:,ii)];
                products = [products, total_complexes(:,jj)];
            end
        end
    end
    
    
    occ_edge_num = sum(edgeTF(1,:));
    occ_node_num = sum(sum(adj_matrix,1)+sum(adj_matrix,2)' > 0);
    
    % the below condition indicate that the given network is the empty
    % network.
    if sum(size(sources)) == 0
        cum_num_empty= cum_num_empty+1;
        % continue
        empty_TF = 1;
    end
    
    
    if empty_TF == 0
        stoi_vec = products - sources;
        stoi_dim = rank(stoi_vec);
        
        [S1,S2] = CRN_countlinkage(sources, products);
        % S1: the number of strongly connected components
        % S2: the number of linkage classes
        num_complexes = size(unique([sources, products]','rows')',2);
        deficiency = num_complexes - S2 - stoi_dim;
        
        % Line 128-141 counts the number of networks having WR, ZD, or both
        % without translation
        if deficiency == 0 && S1 == S2
            cum_num_notrans_WR = cum_num_notrans_WR + 1;
            cum_num_notrans_ZD = cum_num_notrans_ZD + 1;
            cum_num_notrans_WR_ZD = cum_num_notrans_WR_ZD + 1;
            wr_TF = 1;
            dz_TF = 1;
            wrdz_TF = 1;
        elseif deficiency == 0
            cum_num_notrans_ZD = cum_num_notrans_ZD + 1;
            dz_TF = 1;
        elseif S1 == S2
            cum_num_notrans_WR = cum_num_notrans_WR + 1;
            wr_TF = 1;
        end
        
        
        
        % having ZD after translation (Theorem ..) in the paper.
        
        stoi_origin = products - sources;
        
        num_base_edge = size(base_edge, 2);
        
        stoi_origin_extend = [stoi_origin, -stoi_origin];
        
        % The below code (line 154 - 186) checks the necessary condition
        % for having ZD after translation (Theorem 3.3 and 3.5).
        base_edge_id = ismember(base_edge', stoi_origin_extend', 'rows');
        % unique(stoi_origin', 'rows')  
        comp_reject  = 0;
        if sum(base_edge_id) > num_spec * (2*num_spec - 1) %(Theorem 3.5)
            comp_reject = 1; 
        elseif sum(base_edge_id) > 2 * num_spec %(Theorem 3.3)
            comp_reject = 1;
            for ss = 3:(2*num_spec)
                test_list = nchoosek(find(base_edge_id), ss);
                sign_matrix = ones(2^(ss-1), ss-1);
                for rr = 0:(2^(ss-1) - 1)
                    sign_matrix(rr+1,:) = 1 -  2 * de2bi_hh(rr, ss-1);
                end
                
                for rr = 1:size(test_list, 1)
                    base_edge_test_id = test_list(rr,:);
                    base_edge_test = base_edge(:, base_edge_test_id);
                    for sg = 1:2^(ss-1)
                        base_edge_test_sign = repmat([1,sign_matrix(sg,:)], [num_spec, 1]) .* base_edge_test;
                        if isequal(sum(base_edge_test_sign,2), zeros(num_spec, 1))
                            comp_reject = 0;
                            break
                        end
                    end
                    if comp_reject == 0
                        break
                    end
                end
                if comp_reject == 0
                    break
                end
            end
        end


        Solution = {};
        Index = {};

        
        if size(stoi_origin,2) == rank(stoi_origin) & comp_reject
            Solution = {};
            Index = {};
            cum_num_filtered_by_pos_sol = cum_num_filtered_by_pos_sol + 1;
            cum_num_filtered_by_num_comp = cum_num_filtered_by_num_comp + 1;
        else
            % The below code (line 201) checks the necessary condition for
            % having WR after translation (Theorem 3.1) in the paper.
            [pos_sol, pos_sol_TF] = CRN_positive_solution_optim(stoi_origin);
            
            if pos_sol_TF == 0 & comp_reject
                Solution = {};
                Index = {};
                cum_num_filtered_by_pos_sol = cum_num_filtered_by_pos_sol + 1;
                cum_num_filtered_by_num_comp = cum_num_filtered_by_num_comp + 1;
            else
                wr_trans_TF = -cum_num_trans_WR;
                dz_trans_TF = -cum_num_trans_ZD;
                wrdz_trans_TF = -cum_num_trans_WR_ZD;
                
                % Perform network translation (Step 2 & 3)
                
                [Solution, Index, cum_num_cand, cum_num_cand_filtered_by_def, cum_num_trans_WR, cum_num_trans_ZD, cum_num_trans_WR_ZD] ...
                    = CRN_translation_v1(sources, products, max_order, cum_num_cand, cum_num_cand_filtered_by_def, cum_num_trans_WR, cum_num_trans_ZD, cum_num_trans_WR_ZD);  % Merging reactions
                
                wr_trans_TF = wr_trans_TF+cum_num_trans_WR;
                dz_trans_TF = dz_trans_TF+cum_num_trans_ZD;
                wrdz_trans_TF = wrdz_trans_TF+cum_num_trans_WR_ZD;
                
                % Solution contains the source and product complex
                % matrices of found translated CRNs with WR and ZD.
                % Index contains the permutation of reactions between the
                % original CRN and the found CRNs.
                
            end
        end
        cum_num_WR_ZD = cum_num_WR_ZD + numel(Solution)/2;
    end
    
end
timespent0 = toc;
save(['TOWARDZ_results.mat']);
end



