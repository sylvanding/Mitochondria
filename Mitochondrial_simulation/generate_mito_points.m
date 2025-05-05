% 模拟线粒体图像点坐标并保存为 CSV 文件

num_iterations = 1;

for k = 1:num_iterations
    % --- 参数定义 (来自 call_mito.m) ---
    px = 160; %nm
    NA = 1.4; %Numerica aperture
    dens = 1; %Labelling density
    FOV = 12; %Field of view in microns
    nframes = 500; %number of frames
    anti_length = 20; %Length of antibody in nm; double it if you are using primary and secondary
    mito_diam = 0.5; %Circular diameter of mitochondria (probably should not change)
    ep_dens = 15.49; %1 alpha/beta mitochondria epitope per 14.49 nm slong length of microtubule
    k2 = 0.0135; %Chance of moving from on state to triplet state ms^-1
    k3 = 0.333e-04; %Chance of moving from triplet state to on state ms^-1
    k4 = 0.0011; %Chance of moving from on state to bleached state ms^-1
    frame_rate = 20; %ms
    zeroth_int = 500; %Desired photon count of zeroth order
    first_int = 1500; %Desired photon count of first order
    Ibgp = 15; %Bacground noise of zeroth order (photons*100)
    Ibg = 3 * Ibgp; %Background noise of first order (photons*100)
    RN = 100; %readout noise (photons*100)
    mito_length = 1.0; %um
    heterogeneity = 1; %Heterogeneity of mitochondria sizes; scales from 0 to 1
    numChains = 3; %number of mitochondria
    dissociation = 0;
    % --- 参数定义结束 ---

    % --- 参数定义 (新增参数) ---
    output_path = 'mito_points';
    if ~exist(output_path, 'dir')
        mkdir(output_path);
    end
    % --- 参数定义结束 ---

    % --- 模拟生成线粒体结构 ---
    [mito, mito_edge, mito_viewable, tempx, tempy, zpos, spec_hetr] = sim_mito_3D(dens, px, numChains, anti_length, mito_diam, ep_dens, FOV, mito_length, heterogeneity);
    % --- 模拟生成结束 ---

    % --- 生成用于点坐标提取的中间表示 fin_gt2 ---
    in_mit_fact = 1;
    edge_mit = 1;
    mito_edge3 = zeros(size(mito_edge)); % Preallocate mito_edge3
    fin_gt2 = zeros(size(mito_edge)); % Preallocate fin_gt2

    for n = 1:size(mito_edge, 3)
        inner = imfill(mito_edge(:, :, n));
        inner = inner .* (rand(size(mito_edge, 1), size(mito_edge, 2)) > 0.975);
        mito_edge3(:, :, n) = mito_edge(:, :, n) .* (rand(size(mito_edge, 1), size(mito_edge, 2)) > 0.55);
        inner(inner > 0) = in_mit_fact;
        mito_edget = mito_edge3(:, :, n) * edge_mit;
        % backgr = rand(size(mito_edge3, 1), size(mito_edge3, 2)) > 0.998;
        % fin_gt2(:, :, n) = mito_edget + inner + backgr;
        fin_gt2(:, :, n) = mito_edget + inner;
    end

    % --- fin_gt2 生成结束 ---

    % --- 从 fin_gt2 提取点坐标 ---
    vals = unique(fin_gt2);
    gt = [];
    GT_list = [];
    chance = [];
    total_points = 0; % 计算总点数以预分配 GT_list

    for n = 1:length(vals)
        num_points = length(find(fin_gt2 == vals(n))) * double(vals(n));
        chance(n) = num_points;
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

    [x, y, z] = ind2sub([size(fin_gt2, 1), size(fin_gt2, 2), size(fin_gt2, 3)], GT_list);

    x = 16 * x; y = 16 * y; z = 16 * z; % 转换为 nm

    % x = 16 * (x + rand(size(x)) - .5); y = 16 * (y + rand(size(y)) - .5); z = 16 * (z + rand(size(z)) - .5); % 转换为 nm 并添加随机偏移
    % x = x + dissociation * randn(length(x), 1);
    % y = y + dissociation * randn(length(x), 1);
    % z = z + dissociation * randn(length(x), 1);
    % z0 = mean(z);
    % --- 点坐标提取结束 ---

    % --- 生成闪烁和最终坐标 TR0 (包含不确定性) ---
    nblinks = 30000; % 定义闪烁次数
    order = randperm(length(y)); % 随机排列点
    mu = 3; % 平均闪烁持续时间 (帧数)
    blinks = round(exprnd(mu, nblinks, 1)); % 生成指数分布的闪烁持续时间
    blinks = blinks(blinks > 0); % 移除 0 持续时间的闪烁
    actual_blinks = sum(blinks); % 实际总闪烁帧数

    indices = zeros(actual_blinks, 1); % 预分配 indices
    current_pos = 1;
    blink_count = 1;

    for n = 1:length(blinks)

        if blink_count <= length(order) % 确保不会超出 order 的范围
            num_blinks_for_point = blinks(n);
            indices(current_pos:current_pos + num_blinks_for_point - 1) = repmat(order(blink_count), num_blinks_for_point, 1);
            current_pos = current_pos + num_blinks_for_point;
            blink_count = blink_count + 1;
        else
            % 如果 order 中的点用完了，就停止生成 indices
            indices = indices(1:current_pos - 1); % 裁剪未使用的部分
            break;
        end

    end

    xn = x(indices); % 根据闪烁索引选择 x 坐标
    yn = y(indices); % 根据闪烁索引选择 y 坐标
    zn = z(indices); % 根据闪烁索引选择 z 坐标

    TR0gt = (1:length(indices))'; % Ground truth 帧 ID
    TR0gt(:, 2) = ones(length(indices), 1); % Ground truth (通常是粒子 ID 或类别，这里全为 1)
    TR0gt(:, 3:5) = [xn, yn, zn]; % Ground truth 坐标 (nm)

    TR0 = TR0gt; % 复制一份作为带有不确定性的坐标的基础

    % lat_unc = (15.96 + 10.37 * randn(length(TR0gt(:, 3)), 1)); % 计算横向不确定度 (nm)
    % ax_unc = 2.5641 * lat_unc; % 计算轴向不确定度 (nm)

    % 添加高斯噪声模拟定位不确定性
    % TR0(:, 3) = TR0gt(:, 3) + lat_unc .* randn(length(TR0gt(:, 3)), 1); % 添加 x 不确定性
    % TR0(:, 4) = TR0gt(:, 4) + lat_unc .* randn(length(TR0gt(:, 3)), 1); % 添加 y 不确定性
    % TR0(:, 5) = TR0gt(:, 5) + ax_unc .* randn(length(TR0gt(:, 3)), 1); % 添加 z 不确定性
    % --- TR0 生成结束 ---

    % --- 保存 TR0 数据到 CSV 文件 ---
    A6 = {'id', 'frame', 'x [nm]', 'y [nm]', 'z [nm]'}; % 定义 CSV 文件头
    savefile = sprintf('%s/%d.csv', output_path, k); % 定义输出文件名，格式为 k.csv
    writecell(A6, savefile); % 写入文件头
    dlmwrite(savefile, TR0, 'delimiter', ',', '-append'); % 追加写入数据
    % --- 保存结束 ---

    fprintf('已生成并保存文件: %s\n', savefile);
    clear;

end

disp('所有模拟完成！');
