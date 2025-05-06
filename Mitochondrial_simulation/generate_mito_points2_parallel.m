% 模拟线粒体图像点坐标并保存为 CSV 文件 (并行版本)

function generate_mito_points2_parallel()
    % 模拟线粒体图像点坐标并保存为 CSV 文件 (并行版本)

    clear; % 清除函数工作区或基础工作区（取决于调用方式）

    num_iterations = 4; % 可以根据需要修改迭代次数

    % --- 定义所有将传递给工作函数的参数 ---
    sim_params.px = 160; %nm
    sim_params.NA = 1.4; %Numerica aperture
    sim_params.dens = 1; %Labelling density
    sim_params.FOV = 6; %Field of view in microns
    sim_params.anti_length = 20; %Length of antibody in nm; double it if you are using primary and secondary
    sim_params.mito_diam = 0.5; %Circular diameter of mitochondria (probably should not change)
    sim_params.ep_dens = 15.49; %1 alpha/beta mitochondria epitope per 14.49 nm slong length of microtubule
    % sim_params.k2 = 0.0135; %Chance of moving from on state to triplet state ms^-1
    % sim_params.k3 = 0.333e-04; %Chance of moving from triplet state to on state ms^-1
    % sim_params.k4 = 0.0011; %Chance of moving from on state to bleached state ms^-1
    % sim_params.frame_rate = 20; %ms
    % sim_params.zeroth_int = 500; %Desired photon count of zeroth order
    % sim_params.first_int = 1500; %Desired photon count of first order
    % sim_params.Ibgp = 15; %Bacground noise of zeroth order (photons*100)
    % sim_params.Ibg = 3 * sim_params.Ibgp; %Background noise of first order (photons*100)
    % sim_params.RN = 100; %readout noise (photons*100)
    sim_params.mito_length = 1.0; %um
    sim_params.heterogeneity = 1; %Heterogeneity of mitochondria sizes; scales from 0 to 1
    sim_params.numChains = 12; %number of mitochondria
    % sim_params.target_point_count = 16384; % target number of points (注释掉，与原脚本一致)
    % sim_params.dissociation = 0; % (注释掉，与原脚本一致)

    sim_params.output_path = 'mito_points'; % 输出路径

    % --- 创建输出目录 (在循环开始前执行一次) ---
    if ~exist(sim_params.output_path, 'dir')
        mkdir(sim_params.output_path);
    end

    % --- 并行循环 ---
    % 先检查是否有现有的并行池，如果有则删除
    existingPool = gcp('nocreate');

    if ~isempty(existingPool)
        delete(existingPool);
        disp('已删除现有并行池');
    end

    % 确保 MATLAB Parallel Computing Toolbox 可用
    % 可以显式启动并行池:
    pool = parpool(4);

    try

        parfor k_iter = 1:num_iterations
            process_single_simulation(k_iter, sim_params);
        end

    catch ME
        disp('并行循环中发生错误:');
        disp(ME.message);
        disp('错误堆栈跟踪:');
        disp(ME.stack);

    end

    if ~isempty(pool)
        delete(pool);
        disp('已关闭并行池。');
    end

    disp('所有模拟完成 (并行版本)！');

    clear; % 再次清除函数工作区或基础工作区

end

% --- 工作函数 (本地函数) ---
function process_single_simulation(k, params)
    % 此函数包含单次迭代的核心逻辑

    % 从结构体中解包参数
    px = params.px;
    NA = params.NA;
    dens = params.dens;
    FOV = params.FOV;
    anti_length = params.anti_length;
    mito_diam = params.mito_diam;
    ep_dens = params.ep_dens;
    mito_length = params.mito_length;
    heterogeneity = params.heterogeneity;
    numChains = params.numChains;
    output_path = params.output_path;
    % target_point_count = params.target_point_count; % 如果使用采样
    % dissociation = params.dissociation; % 如果使用

    % --- 模拟生成线粒体结构 ---
    % sim_mito_3D_2 函数必须在 MATLAB 路径上，或者作为本地/嵌套函数定义
    [~, mito_edge, ~, ~, ~, ~, ~, mito_label_map] = sim_mito_3D_2(dens, px, numChains, anti_length, mito_diam, ep_dens, FOV, mito_length, heterogeneity);

    % --- 生成用于点坐标提取的中间表示 fin_gt2 ---
    edge_mit = 1;
    mito_edge3 = zeros(size(mito_edge)); % Preallocate mito_edge3
    fin_gt2 = zeros(size(mito_edge)); % Preallocate fin_gt2

    for n_loop_gt2 = 1:size(mito_edge, 3) % 重命名循环变量 n
        mito_edge3(:, :, n_loop_gt2) = mito_edge(:, :, n_loop_gt2) .* (rand(size(mito_edge, 1), size(mito_edge, 2)) > 0.55);
        mito_edget = mito_edge3(:, :, n_loop_gt2) * edge_mit;
        fin_gt2(:, :, n_loop_gt2) = mito_edget; % 修改：只使用 edge 点
    end

    % --- 从 fin_gt2 提取点坐标 ---
    vals = unique(fin_gt2);
    total_points = 0;

    for n_val_loop = 1:length(vals)
        val_current = vals(n_val_loop);

        if val_current == 0 % 跳过值为0的情况，因为它不贡献点数
            continue;
        end

        num_points_for_val = length(find(fin_gt2 == val_current)) * double(val_current);
        total_points = total_points + num_points_for_val;
    end

    GT_list = zeros(total_points, 1); % 预分配 GT_list
    current_pos = 1;

    for n_val_loop = 1:length(vals)
        val_current = vals(n_val_loop);

        if val_current == 0 % 再次跳过值为0的情况
            continue;
        end

        indices_val = find(fin_gt2 == val_current);
        num_repeats = val_current;

        if num_repeats > 0 && ~isempty(indices_val)
            points_to_add = repmat(indices_val, [double(num_repeats), 1]);
            len_points_to_add = length(points_to_add);

            if (current_pos + len_points_to_add - 1) <= total_points % 确保不超出预分配的 GT_list 大小
                GT_list(current_pos:(current_pos + len_points_to_add - 1)) = points_to_add;
                current_pos = current_pos + len_points_to_add;
            else
                % fprintf('警告 (迭代 %d): GT_list 空间不足，可能总点数计算有误或 repmat 结果超出预期。\n', k);
                % 理论上，如果 total_points 计算正确，这里不应发生
                % 如果发生，可能需要截断 points_to_add 或重新检查 total_points 计算
                actual_len_to_add = total_points - current_pos + 1;

                if actual_len_to_add > 0 && actual_len_to_add <= len_points_to_add
                    GT_list(current_pos:total_points) = points_to_add(1:actual_len_to_add);
                    current_pos = total_points + 1;
                end

            end

        end

    end

    % 如果 GT_list 由于某种原因没有完全填满 (例如，所有 val_current 都是0)
    % 或者 current_pos 没有达到 total_points + 1，则修剪 GT_list
    if current_pos <= total_points
        GT_list = GT_list(1:current_pos - 1);
    end

    if isempty(GT_list)
        x_idx = [];
        y_idx = [];
        z_idx = [];
    else
        [x_idx, y_idx, z_idx] = ind2sub([size(fin_gt2, 1), size(fin_gt2, 2), size(fin_gt2, 3)], GT_list);
    end

    % --- 新增：从 mito_label_map 获取每个点的线粒体ID ---
    mito_ids = zeros(length(x_idx), 1); % 预分配线粒体ID数组

    if ~isempty(x_idx) % 只有在有点的时候才查找ID

        for i_coord = 1:length(x_idx) % 重命名循环变量 i
            mito_ids(i_coord) = mito_label_map(x_idx(i_coord), y_idx(i_coord), z_idx(i_coord));
        end

    end

    % 将体素索引转换为nm单位坐标
    x_nm = 16 * x_idx;
    y_nm = 16 * y_idx;
    z_nm = 16 * z_idx;

    % --- 采样点到指定数量 (注释部分保持不变) ---
    % current_num_points = length(x_nm);
    % fprintf('迭代 %d: 采样前点数: %d\n', k, current_num_points);
    % if current_num_points > target_point_count && target_point_count > 0
    %     sample_indices = randperm(current_num_points, target_point_count);
    %     x_nm = x_nm(sample_indices);
    %     y_nm = y_nm(sample_indices);
    %     z_nm = z_nm(sample_indices);
    %     mito_ids = mito_ids(sample_indices); % 同步采样 mito_ids
    %     fprintf('迭代 %d: 采样后点数: %d\n', k, length(x_nm));
    % elseif target_point_count <= 0 && current_num_points > 0
    %     fprintf('迭代 %d: 目标点数设置为 %d (<=0)，不进行采样，保留所有点。\n', k, target_point_count);
    %     fprintf('迭代 %d: 采样后点数: %d\n', k, current_num_points);
    % elseif current_num_points == 0
    %     fprintf('迭代 %d: 采样前的点数为0，无需采样。\n', k);
    %     fprintf('迭代 %d: 采样后的点数: 0\n', k);
    % else % current_num_points <= target_point_count (and target_point_count > 0)
    %     fprintf('迭代 %d: 采样后点数: %d (未进行采样，因点数 %d <= 目标点数 %d 或目标点数无效)\n', k, current_num_points, current_num_points, target_point_count);
    % end

    % --- 保存点坐标数据到 CSV 文件 ---
    savefile = sprintf('%s/%d.csv', output_path, k); % 定义输出文件名

    if ~isempty(x_nm) % 使用 x_nm 检查是否为空
        points_table = table(x_nm, y_nm, z_nm, mito_ids, 'VariableNames', {'x [nm]', 'y [nm]', 'z [nm]', 'MitoID'});
        writetable(points_table, savefile);
        fprintf('已生成并保存文件 (迭代 %d): %s (包含 %d 个点)\n', k, savefile, length(x_nm));
    else
        empty_coords = double.empty(0, 1);
        empty_ids = double.empty(0, 1); % 为空表也定义ID列
        points_table = table(empty_coords, empty_coords, empty_coords, empty_ids, 'VariableNames', {'x [nm]', 'y [nm]', 'z [nm]', 'MitoID'});
        writetable(points_table, savefile);
        fprintf('已生成并保存文件 (迭代 %d): %s (包含 0 个点)\n', k, savefile);
    end

    % parfor 循环会自动处理工作区清理，此处无需 clearvars
end
