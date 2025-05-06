function [mito, mito_edge, mito_viewable, tempx, tempy, zpos, spec_hetr, mito_label_map] = sim_mito_3D_2(dens, px, numChains, anti_length, mito_diam, ep_dens, FOV, mito_length, heterogeneity)
    %numChains=200;
    siz2 = (FOV * 1000) / px;
    mag = 10;
    sizee = siz2 * mag; %image size
    sizee_ax = 3 * mag;
    sizee_ax2 = round(3 / ((px / 1000) / mag));
    size3 = (FOV * 1000) / sizee; %pixel size
    %NA=1.4;
    anti_length = (anti_length / px) * mag;
    %ep_dens=14.89; %1 epitope/16 nm
    ep_dens = 1;
    mito_diam = (mito_diam * 1000) / size3;
    sz = 2000 / px;
    %SRzoom=8;

    persisLength = 600 / size3;
    mito_viewable = zeros(sizee, sizee, sizee_ax2);
    % mito_edge=zeros(sizee); % 旧的 mito_edge 初始化，可以保留或修改其用途
    mito_edge = zeros(sizee, sizee, sizee_ax2); % 确保 mito_edge 是3D的，如果之前不是
    mito_label_map = zeros(sizee, sizee, sizee_ax2); % <--- 新增
    mito = zeros(sizee);
    %PSF = genPSFparam(sz,px,NA,lambda);
    tempx = [];
    tempy = [];
    tempz = [];
    ind1 = 0; ind2 = 0; ind3 = 0;

    for chin = 1:numChains
        blank = zeros(sizee, sizee, sizee_ax2);

        %numVerts = (mag)*lognrnd(log((((mito_length*1000)/size3))),heterogeneity);
        numVerts = ((mito_length * 1000) / size3) * 2 * rand * heterogeneity + ((mito_length * 1000) / size3) * (1 - heterogeneity);
        mic = [.125 * sizee + .75 * sizee * rand, .125 * sizee + .75 * sizee * rand, (.0125 * sizee_ax) * rand + .4938 * sizee_ax2];
        vec_ang = 2 * pi() * rand;
        phi_ang = 0; %0.0313*pi()*rand-.0157*pi();
        radius = (mito_diam / 2);
        radius = radius + randn * radius / 10;
        radius2 = (mito_diam / 2);
        radius2 = radius2 + randn * radius2 / 10;

        for ind = 1:(32 / size3):numVerts
            radius = (32 / size3) * (radius / 60) * randn + radius;
            radius2 = (radius2 / 30) * randn * (32 / size3) + radius2;

            if radius > (mito_diam)
                radius = mito_diam;
            end

            if radius < (mito_diam / 4)
                radius = mito_diam / 4;
            end

            if radius2 > (mito_diam)
                radius2 = mito_diam;
            end

            if radius2 < (mito_diam / 4)
                radius2 = mito_diam / 4;
            end

            if radius < (radius2 / 2)
                radius = radius2 / 2;
            end

            if radius2 < (radius / 2)
                radius2 = radius / 2;
            end

            angle = sqrt((1 / (persisLength * 2))) * (rand - .5);
            phi = sqrt((1 / (persisLength * 16))) * (rand - .5);
            vec_ang = vec_ang + (32 / size3) * angle;
            phi_ang = 0; %phi_ang+(32/size3)*phi;

            rand_disp = (32 / size3) * ep_dens; %*(1+.3*randn(1));
            x = cos(vec_ang);
            y = sin(vec_ang);
            z = sin(phi);
            fact = x ^ 2 + y ^ 2 + z ^ 2;
            rand_disp = rand_disp / fact;
            mic_new = mic(end, :) + rand_disp * [cos(vec_ang), sin(vec_ang), sin(phi_ang)];
            mic = [mic; mic_new];

            if mic_new(1) < 0
                % ind1=ind1+1 % 这些 ind 变量似乎没有在后续使用，可以考虑移除或正确使用
                break;
            elseif mic_new(1) > sizee
                % ind1=ind1+1
                break;
            else
                x_pos = mic_new(1);
            end

            if mic_new(2) < 0
                % ind2=ind2+1
                break;
            elseif mic_new(2) > sizee
                % ind2=ind2+1
                break;
            else
                y_pos = mic_new(2);
            end

            if mic_new(3) < 0
                % ind3=ind3+1
                break;
            elseif mic_new(3) > sizee_ax2
                % ind3=ind3+1
                break;
            else
                z_pos = mic_new(3);
            end

            shiftx = 0; %(mito_diam*(rand-.5)*sin(vec_ang));
            shifty = 0; %(mito_diam*(rand-.5)*cos(vec_ang));
            shiftz = 0;

            if (floor(sizee * x_pos / sizee) + 1 + shiftx) < 1
                shiftx = 0;
            end

            if (floor(sizee * y_pos / sizee) + 1 + shifty) < 1
                shifty = 0;
            end

            if (floor(sizee * z_pos / sizee_ax) + 1 + shiftz) < 1 % 应该是 sizee_ax2 ?
                shiftz = 0;
            end

            shift2x = ((shiftx / (sizee / siz2)));
            shift2y = ((shifty / (sizee / siz2)));
            shift2z = ((shiftz / (sizee / siz2))); % 应该是 sizee_ax2 ?

            if (floor(siz2 * x_pos / 512) + 1 + shift2x) < 1 % 512 是什么？FOV/px ?
                shift2x = 0;
            end

            if (floor(siz2 * y_pos / 512) + 1 + shift2y) < 1
                shift2y = 0;
            end

            if (floor(siz2 * z_pos / 512) + 1 + shift2z) < 1
                shift2z = 0;
            end

            lateralanti_shift = (rand * anti_length);
            xanti = 2 * (rand - .5) * lateralanti_shift;
            yanti = sign(rand - .5) * sqrt(lateralanti_shift ^ 2 - xanti ^ 2);
            zanti = anti_length * (2 * ((1 - lateralanti_shift / anti_length) - .5));

            if isnan(zanti)
                zanti = 0;
            end

            if (floor(sizee * x_pos / sizee + shiftx + (xanti)) + 1) < 1
                xanti = 0;
            end

            if (floor(sizee * y_pos / sizee + shifty + (yanti)) + 1) < 1
                yanti = 0;
            end

            if (floor(sizee * z_pos / sizee + shiftz + (zanti)) + 1) < 1 % 应该是 sizee_ax2 ?
                zanti = 0;
            end

            xanti2 = ((xanti / (sizee / siz2)));
            yanti2 = ((yanti / (sizee / siz2)));
            zanti2 = ((zanti / (sizee_ax / siz2))); % 应该是 sizee_ax2 ?

            if (floor(siz2 * x_pos / sizee + shift2x + xanti2) + 1) < 1
                xanti2 = 0;
            end

            if (floor(siz2 * y_pos / sizee + shift2y + yanti2) + 1) < 1
                yanti2 = 0;
            end

            if (floor(siz2 * y_pos / sizee_ax + shift2y + yanti2) + 1) < 1 % 应该是 sizee_ax2 ?
                zanti2 = 0;
            end

            % 确保 mape 是3D的
            mape = crop_sphere(sizee, sizee, sizee_ax2, floor(sizee * x_pos / sizee + shiftx + (xanti)) + 1, floor(sizee * y_pos / sizee + shifty + (yanti)) + 1, floor(sizee_ax * z_pos / sizee_ax + shiftz + (zanti)) + 1, radius, radius2); % z_pos/sizee_ax 应该对应 sizee_ax2 的尺度

            blank = blank + mape; % blank 应该是3D的
            % mito 的更新也应该是3D的
            mito(floor(sizee * x_pos / sizee + shiftx + (xanti)) + 1, floor(sizee * y_pos / sizee + shifty + (yanti)) + 1, floor(sizee * z_pos / sizee + shiftz + (zanti)) + 1) = rand < dens; % z_pos/sizee 应该对应 sizee_ax2 的尺度
            % mito2 的更新也应该是3D的
            mito2(floor(siz2 * x_pos / sizee + shift2x + xanti2) + 1, floor(siz2 * y_pos / sizee + shift2y + yanti2) + 1, floor(siz2 * z_pos / sizee + shift2z + zanti2) + 1) = rand < dens; % z_pos/sizee 应该对应 sizee_ax2 的尺度
            tempx = [tempx; siz2 * x_pos / sizee + shift2x + xanti2];
            tempy = [tempy; siz2 * y_pos / sizee + shift2y + yanti2];
            tempz = [tempz; siz2 * z_pos / 512 + shift2z + zanti2]; % z_pos/512 ?
            % mito_viewable 的索引也需要是3D
            if mito_viewable(floor(x_pos) + 1, floor(y_pos) + 1, floor(z_pos) + 1) > 0 % z_pos 应该对应 sizee_ax2 的尺度
                break;
            end

        end

        %blank=conv2(blank,mape,'same')>0; % conv2 是2D卷积，如果blank是3D，这里需要调整
        blank = blank > 0;
        blank = logical((blank - mito_viewable) > 0); % 确保 blank 和 mito_viewable 维度一致
        mito_viewable = mito_viewable + blank; % 确保 blank 和 mito_viewable 维度一致

        current_mito_edge_mask = zeros(size(blank)); % <--- 新增

        for n = 1:size(blank, 3)
            current_mito_edge_mask(:, :, n) = edge(blank(:, :, n), 'canny'); % <--- 修改
        end

        % 根据重叠策略选择以下之一：
        % 策略1: 新ID覆盖旧ID
        mito_label_map(current_mito_edge_mask > 0) = chin; % <--- 新增

        % 策略2: 保留旧ID，只标记未被标记的点
        % unlabeled_edge_points = (current_mito_edge_mask > 0) & (mito_label_map == 0);
        % mito_label_map(unlabeled_edge_points) = chin;

        % 更新总的、不区分ID的 mito_edge (如果还需要的话)
        % 确保 current_mito_edge_mask 和 mito_edge 维度匹配
        mito_edge = mito_edge | (current_mito_edge_mask > 0); % <--- 修改/新增 (确保 mito_edge 是3D的)

        %final(:,:,chin)=mic;
    end

    %mito_diffraction=conv2(mito2,PSF,'same'); % PSF 未定义， mito2 维度?
    zpos = 0.01 * randn(length(tempx), 1);
    spec_hetr = randn(length(tempx), 1);
end
