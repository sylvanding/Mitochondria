% 模拟线粒体图像点坐标并保存为 CSV 文件

clear;

num_iterations = 1;

for k = 1:num_iterations
    % --- 参数定义 (来自 call_mito.m) ---
    px = 160; %nm
    NA = 1.4; %Numerica aperture
    dens = 1; %Labelling density
    FOV = 6; %Field of view in microns
    % nframes = 500; %number of frames
    anti_length = 20; %Length of antibody in nm; double it if you are using primary and secondary
    mito_diam = 0.5; %Circular diameter of mitochondria (probably should not change)
    ep_dens = 15.49; %1 alpha/beta mitochondria epitope per 14.49 nm slong length of microtubule
    % k2 = 0.0135; %Chance of moving from on state to triplet state ms^-1
    % k3 = 0.333e-04; %Chance of moving from triplet state to on state ms^-1
    % k4 = 0.0011; %Chance of moving from on state to bleached state ms^-1
    % frame_rate = 20; %ms
    % zeroth_int = 500; %Desired photon count of zeroth order
    % first_int = 1500; %Desired photon count of first order
    % Ibgp = 15; %Bacground noise of zeroth order (photons*100)
    % Ibg = 3 * Ibgp; %Background noise of first order (photons*100)
    % RN = 100; %readout noise (photons*100)
    mito_length = 1.0; %um
    heterogeneity = 1; %Heterogeneity of mitochondria sizes; scales from 0 to 1
    numChains = 12; %number of mitochondria
    target_point_count = 16384; % target number of points
    % dissociation = 0;
    % --- 参数定义结束 ---

    % --- 参数定义 (新增参数) ---
    output_path = 'mito_points';

    if ~exist(output_path, 'dir')
        mkdir(output_path);
    end

    % --- 参数定义结束 ---

    % --- 模拟生成线粒体结构 ---
    [mito, mito_edge, mito_viewable, tempx, tempy, zpos, spec_hetr, mito_label_map] = sim_mito_3D_2(dens, px, numChains, anti_length, mito_diam, ep_dens, FOV, mito_length, heterogeneity);
    % --- 模拟生成结束 ---

    % --- 生成用于点坐标提取的中间表示 fin_gt2 ---
    % in_mit_fact = 1; % 此行已移除，因为 inner 点不再生成
    edge_mit = 1;
    mito_edge3 = zeros(size(mito_edge)); % Preallocate mito_edge3
    fin_gt2 = zeros(size(mito_edge)); % Preallocate fin_gt2

    for n = 1:size(mito_edge, 3)
        % inner = imfill(mito_edge(:, :, n)); % inner 点相关代码已移除
        % inner = inner .* (rand(size(mito_edge, 1), size(mito_edge, 2)) > 0.975); % inner 点相关代码已移除
        mito_edge3(:, :, n) = mito_edge(:, :, n) .* (rand(size(mito_edge, 1), size(mito_edge, 2)) > 0.55);
        % inner(inner > 0) = in_mit_fact; % inner 点相关代码已移除
        mito_edget = mito_edge3(:, :, n) * edge_mit;
        % backgr = rand(size(mito_edge3, 1), size(mito_edge3, 2)) > 0.998;
        % fin_gt2(:, :, n) = mito_edget + inner + backgr;
        fin_gt2(:, :, n) = mito_edget; % 修改：只使用 edge 点
    end

    % --- fin_gt2 生成结束 ---

    % --- 从 fin_gt2 提取点坐标 ---
    vals = unique(fin_gt2);
    % gt = []; % 此变量未使用，已移除
    GT_list = [];
    % chance = []; % 此变量未使用，已移除
    total_points = 0; % 计算总点数以预分配 GT_list

    for n = 1:length(vals)
        num_points = length(find(fin_gt2 == vals(n))) * double(vals(n));
        total_points = total_points + num_points;
    end

    GT_list = zeros(total_points, 1); % 预分配 GT_list
    current_pos = 1;

    for n = 1:length(vals)
        indices_val = find(fin_gt2 == vals(n));
        num_repeats = vals(n);

        if num_repeats > 0 && ~isempty(indices_val)
            points_to_add = repmat(indices_val, [num_repeats, 1]);
            GT_list(current_pos:current_pos + length(points_to_add) - 1) = points_to_add;
            current_pos = current_pos + length(points_to_add);
        end

    end

    [x_idx, y_idx, z_idx] = ind2sub([size(fin_gt2, 1), size(fin_gt2, 2), size(fin_gt2, 3)], GT_list);

    % --- 新增：从 mito_label_map 获取每个点的线粒体ID ---
    % 确保 mito_label_map 和 fin_gt2 (从中获得索引) 的维度是匹配的
    % 或者说，确保 sim_mito_3D_2 输出的 mito_label_map 的维度
    % 与 fin_gt2 的维度一致。
    % fin_gt2 的大小由 mito_edge 决定，而 mito_edge 来自 sim_mito_3D_2
    
    mito_ids = zeros(length(x_idx), 1); % 预分配线粒体ID数组
    if ~isempty(x_idx) % 只有在有点的时候才查找ID
        for i = 1:length(x_idx)
            % 使用体素索引 x_idx, y_idx, z_idx 从 mito_label_map 获取ID
            % 假设 mito_label_map 的维度顺序与 ind2sub 的输出一致：
            % (dim1_idx corresponds to x_idx, dim2_idx to y_idx, dim3_idx to z_idx)
            mito_ids(i) = mito_label_map(x_idx(i), y_idx(i), z_idx(i));
        end
    end
    % --- 线粒体ID获取结束 ---

    % 将体素索引转换为nm单位坐标
    x_nm = 16 * x_idx;
    y_nm = 16 * y_idx;
    z_nm = 16 * z_idx;

    % % --- 新增：采样点到指定数量 ---
    % current_num_points = length(x_nm);
    % fprintf('采样前的点数: %d\n', current_num_points);

    % if current_num_points > target_point_count && target_point_count > 0
    %     sample_indices = randperm(current_num_points, target_point_count);
    %     x_nm = x_nm(sample_indices);
    %     y_nm = y_nm(sample_indices);
    %     z_nm = z_nm(sample_indices);
    %     mito_ids = mito_ids(sample_indices); % 同步采样 mito_ids
    %     fprintf('采样后的点数: %d\n', length(x_nm));
    % elseif target_point_count <= 0 && current_num_points > 0
    %     fprintf('目标点数设置为 %d (<=0)，不进行采样，保留所有点。\n', target_point_count);
    %     fprintf('采样后的点数: %d\n', current_num_points);
    %     % 如果不采样，x_nm, y_nm, z_nm, mito_ids 保持不变
    % elseif current_num_points == 0
    %     fprintf('采样前的点数为0，无需采样。\n');
    %     fprintf('采样后的点数: 0\n');
    %     % x_nm, y_nm, z_nm, mito_ids 此时应为空
    % else % current_num_points <= target_point_count (and target_point_count > 0)
    %     fprintf('采样后的点数: %d (未进行采样，因点数 %d <= 目标点数 %d 或目标点数无效)\n', current_num_points, current_num_points, target_point_count);
    %     % 如果不采样，x_nm, y_nm, z_nm, mito_ids 保持不变
    % end
    % % --- 采样结束 ---

    % x = 16 * (x + rand(size(x)) - .5); y = 16 * (y + rand(size(y)) - .5); z = 16 * (z + rand(size(z)) - .5); % 转换为 nm 并添加随机偏移
    % x = x + dissociation * randn(length(x), 1);
    % y = y + dissociation * randn(length(x), 1);
    % z = z + dissociation * randn(length(x), 1);
    % z0 = mean(z);
    % --- 点坐标提取结束 ---

    % --- 生成闪烁和最终坐标 TR0 (包含不确定性) ---
    % 此部分代码已根据您的要求移除 (原行号 85-143)
    % --- TR0 生成结束 ---

    % --- 保存点坐标数据到 CSV 文件 ---
    savefile = sprintf('%s/%d.csv', output_path, k); % 定义输出文件名

    if ~isempty(x_nm) % 使用采样后的 x_nm 检查是否为空
        % mito_ids 已经在此处正确生成并采样
        points_table = table(x_nm, y_nm, z_nm, mito_ids, 'VariableNames', {'x [nm]', 'y [nm]', 'z [nm]', 'MitoID'});
        writetable(points_table, savefile);
        fprintf('已生成并保存文件: %s (包含 %d 个点)\n', savefile, length(x_nm));
    else
        % 如果没有点，创建一个包含表头的空CSV文件
        empty_coords = double.empty(0, 1);
        empty_ids = double.empty(0,1); % 为空表也定义ID列
        points_table = table(empty_coords, empty_coords, empty_coords, empty_ids, 'VariableNames', {'x [nm]', 'y [nm]', 'z [nm]', 'MitoID'});
        writetable(points_table, savefile);
        fprintf('已生成并保存文件: %s (包含 0 个点)\n', savefile);
    end

    % --- 保存结束 ---

    clearvars -except num_iterations;

end

clear;

disp('所有模拟完成！');
