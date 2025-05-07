import numpy as np
import pandas as pd
import re
import os
import pathlib
import math
import multiprocessing
import h5py


def normalize_pc_pair(input, gt=None, return_params=False):
    """
    对点云数据进行归一化,支持单个点云或点云对的归一化

    Args:
        input: 输入点云 (b, n, 3)
        gt: 目标点云 (b, m, 3), 可选
        return_params: 是否返回归一化参数

    Returns:
        normalized_input: 归一化后的输入点云
        normalized_gt: 归一化后的目标点云(如果提供)
        params: 归一化参数(如果return_params=True)
    """
    # 计算中心点和最远距离
    centroid = np.mean(input, axis=1, keepdims=True)  # (b, 1, 3)
    input_centered = input - centroid
    furthest_distances = np.amax(
        np.abs(input_centered), axis=(1, 2), keepdims=True
    )  # (b, 1, 1)
    furthest_distances = np.repeat(furthest_distances, 3, axis=2)

    # 归一化input
    normalized_input = input_centered / furthest_distances

    # 如果提供了gt,使用相同参数归一化
    normalized_gt = None
    if gt is not None:
        gt_centered = gt - centroid
        normalized_gt = gt_centered / furthest_distances

    if return_params:
        params = {"centroid": centroid, "scale": furthest_distances}
        return normalized_input, normalized_gt, params

    return normalized_input, normalized_gt


def pc_random_sample(pc_data, n_points_per_pc):
    """
    对点云数据进行随机采样，支持输入为 numpy 数组或列表

    Args:
        pc_data: 输入点云数据，可以是形状为 (num_point_clouds, num_points, 3) 的 numpy 数组，
                也可以是包含多个点云的列表，每个点云形状为 (num_points, 3)
        n_points_per_pc: 每个点云的目标点数

    Returns:
        sampled_pc_data: 采样后的点云数据，形状为 (num_point_clouds, n_points_per_pc, 3) 的 numpy 数组
    """
    # 检查输入类型并获取点云数量
    if isinstance(pc_data, list):
        num_point_clouds = len(pc_data)
    else:  # 假设是 numpy 数组
        num_point_clouds = pc_data.shape[0]

    # 初始化输出数组
    sampled_pc_data = np.zeros(
        (num_point_clouds, n_points_per_pc, 3), dtype=np.float32
    )

    for i in range(num_point_clouds):
        # 获取当前点云的点数
        if isinstance(pc_data, list):
            current_points = pc_data[i].shape[0] if pc_data[i] is not None else 0
            current_pc = pc_data[i] if pc_data[i] is not None else np.zeros((0, 3))
        else:
            current_points = pc_data.shape[1]
            current_pc = pc_data[i]

        if current_points == 0:
            continue  # 如果没有点，直接跳过，保持全零数据

        if current_points < n_points_per_pc:
            # 如果点数少于目标点数，进行重复采样
            repeat_factor = (
                n_points_per_pc + current_points - 1
            ) // current_points  # 向上取整
            remainder = n_points_per_pc % current_points
            if remainder == 0:
                sampled_data = np.tile(current_pc, (repeat_factor, 1))
            else:
                full_repeats = np.tile(current_pc, (repeat_factor - 1, 1))
                extra_indices = np.random.choice(
                    current_points, remainder, replace=False
                )
                extra_points = current_pc[extra_indices]
                sampled_data = np.vstack([full_repeats, extra_points])
        elif current_points > n_points_per_pc:
            # 如果点数多于目标点数，进行随机采样
            sampled_indices = np.random.choice(
                current_points, n_points_per_pc, replace=False
            )
            sampled_data = current_pc[sampled_indices]
        else:
            # 如果点数等于目标点数，不做任何操作
            sampled_data = current_pc

        sampled_pc_data[i] = sampled_data

    return sampled_pc_data


def read_mito_data(filename):
    """读取mitochondrial source data文件"""
    mito_data = pd.read_csv(filename)
    mito_data = mito_data.to_numpy()
    mito_data = mito_data[:, 0:3]  # 只保留x, y, z坐标

    current_points = len(mito_data)
    if current_points == 0:
        return None  # 如果没有点，直接返回None

    return mito_data


def filter_and_pad_point_cloud(pc_data, z_threshold=(-25, 100)):
    """
    过滤掉点云中 Z 坐标大于阈值的点，并复制点以恢复到原始点数。

    Args:
        pc_data (np.ndarray): 输入的点云数据，形状为 (num_point_clouds, N, 3)。
        z_threshold (float): Z 坐标的过滤阈值。默认为 0.25。

    Returns:
        np.ndarray: 处理后的点云数据，形状与输入相同 (num_point_clouds, N, 3)。
        对于每个点云，Z > z_threshold 的点被移除，
        然后通过随机复制剩余点的方式将点数恢复到 N。
        如果过滤后没有剩余点，则返回一个由零填充的点云。
    """
    num_point_clouds, original_num_points, _ = pc_data.shape
    processed_pcs = []

    for i in range(num_point_clouds):
        single_pc = pc_data[i]

        # 1. 过滤点
        mask = (single_pc[:, 2] >= z_threshold[0]) & (single_pc[:, 2] <= z_threshold[1])
        filtered_pc = single_pc[mask]

        num_filtered_points = filtered_pc.shape[0]

        # 2. 复制点以恢复到 original_num_points
        if num_filtered_points == original_num_points:
            # 没有点被过滤
            processed_pcs.append(filtered_pc)
        elif num_filtered_points == 0:
            # 所有点都被过滤掉了，用零填充
            processed_pcs.append(
                np.zeros((original_num_points, 3), dtype=pc_data.dtype)
            )
        else:
            # 点数不足，需要复制
            points_to_add = original_num_points - num_filtered_points
            # 从过滤后的点中随机选择索引（允许重复）
            indices_to_duplicate = np.random.choice(
                num_filtered_points, points_to_add, replace=True
            )
            duplicated_points = filtered_pc[indices_to_duplicate]
            # 合并过滤后的点和复制的点
            final_pc = np.vstack((filtered_pc, duplicated_points))
            processed_pcs.append(final_pc)

    return np.stack(processed_pcs, axis=0)


def main():
    # 参数设置
    n_points_per_pc = 16384  # 每个点云的点数

    # 下采样参数
    downsample_rate = 8
    n_points_per_pc_sparse = n_points_per_pc // downsample_rate  # 2048

    # 输出目录和 H5 文件名
    output_file_name = "mito_pc_%s_%s" % (
        n_points_per_pc,
        n_points_per_pc_sparse,
    )
    base_output_dir = pathlib.Path("./outputs")
    base_output_dir.mkdir(parents=True, exist_ok=True)
    # 定义聚合 H5 文件的完整路径
    aggregate_h5_filepath = base_output_dir / f"{output_file_name}.h5"

    # 获取所有微管中轴线的 txt 文件
    input_dir = pathlib.Path("./Mitochondrial_simulation/mito_points")
    txt_files = list(input_dir.glob("*.csv"))

    # 设置进程数
    # num_processes = multiprocessing.cpu_count()
    num_processes = 4
    print(f"使用 {num_processes} 个进程进行并行处理...")

    # 准备传递给 process_file 的参数列表
    args_list = [txt_file for txt_file in txt_files]

    # 创建进程池并执行任务，收集结果（现在是点云数组或 None）
    gt_data = []
    processed_file_names = []  # 存储成功处理的文件名，以便追踪
    with multiprocessing.Pool(processes=num_processes) as pool:
        results = pool.map(read_mito_data, args_list)
        # 过滤掉 None 的结果 (处理失败或未生成点的文件)
        for i, result in enumerate(results):
            if result is not None:
                gt_data.append(result)
                processed_file_names.append(txt_files[i].name)  # 记录成功处理的文件名

    # --- 将所有收集到的点云保存到一个 H5 文件中 ---
    if gt_data:
        print(f"\n成功处理了 {len(gt_data)} 个文件，准备聚合保存...")
        # gt_data 采样
        gt_data = pc_random_sample(gt_data, n_points_per_pc)

        # 将列表中的 NumPy 数组堆叠成一个大的 NumPy 数组
        # 形状将是 (num_files, n_points_per_pc, 3)
        gt_data = np.stack(gt_data, axis=0)
        
        print(f"gt_data.shape: {gt_data.shape}")

        # filter
        # gt_data = filter_and_pad_point_cloud(gt_data, z_threshold=(-25, 100))

        # sample and normalize
        input_data = pc_random_sample(gt_data, n_points_per_pc_sparse)
        gt_data, input_data, norm_params = normalize_pc_pair(gt_data, input_data, True)

        print(f"np.max(input_data): {np.max(input_data, axis=(0, 1))}")
        print(f"np.min(input_data): {np.min(input_data, axis=(0, 1))}")
        print(f"np.max(gt_data): {np.max(gt_data, axis=(0, 1))}")
        print(f"np.min(gt_data): {np.min(gt_data, axis=(0, 1))}")

        try:
            with h5py.File(aggregate_h5_filepath, "w") as f:
                # 保存聚合的点云数据
                f.create_dataset("input_data", data=input_data.astype(np.float32))
                f.create_dataset("gt_data", data=gt_data.astype(np.float32))
                # Store normalization parameters
                f.create_dataset(
                    "norm_params/centroid",
                    data=norm_params["centroid"].astype(np.float32),
                )
                f.create_dataset(
                    "norm_params/scale", data=norm_params["scale"].astype(np.float32)
                )

            print(f"\n所有点云已成功聚合保存到: {aggregate_h5_filepath}")
            print(f"H5 文件包含 'gt_data' 数据集，形状为: {gt_data.shape}")
            print(f"'source_files' 数据集，包含 {len(processed_file_names)} 个源文件。")

        except Exception as save_e:
            print(f"\n聚合保存到 H5 时出错: {save_e}")
    else:
        print("\n未能成功处理任何文件，未生成聚合 H5 文件。")

    print("所有文件处理完成。")


if __name__ == "__main__":
    main()
