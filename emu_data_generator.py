# emu_data_generator.py
import numpy as np
from scipy.special import comb

from scripts.src.core.common.config import CoreConstants as natural_abundance


def generate_input_emu_data(input_emu_dict, input_metabolite_dict=None, custom_labeling=None):
    """
    生成输入EMU数据字典，为每个EMU生成同位素分布

    参数:
    input_emu_dict: EMU方程分析返回的输入EMU需求，包含每个EMU的碳原子数和结构
    input_metabolite_dict: 输入代谢物的同位素标记信息，格式为：
                          {metabolite_name: [{"ratio_list": [0,1,..], "abundance": 0.X}, ...]}
                          例如: {"AcCoA": [{"ratio_list": [0,1], "abundance": 0.25}, ...]}
    natural_abundance: 自然丰度的13C比例 (默认1%)
    custom_labeling: 自定义代谢物标记字典 {emu_name: distribution_array}
                     例如: {'AcCoA__11': np.array([0.5, 0.25, 0.25])}
                     其中分布数组的长度等于碳原子数+1，表示[M+0, M+1, ..., M+n]

    返回:
    input_emu_data_dict: 输入EMU数据字典，格式为 {emu_name: np.array([M+0, M+1, ..., M+n])}
    """

    input_emu_data_dict = {}
    natural_dist_cache = {}

    # 处理每个代谢物的EMU
    for metabolite, emu_value in input_emu_dict.items():
        # 确保emu_value是一个列表
        emu_items = emu_value if isinstance(emu_value, list) else [emu_value]

        for emu in emu_items:
            try:
                emu_name = emu.full_name
                carbon_count = emu.emu_carbon_num
            except AttributeError as e:
                raise AttributeError(f"EMU对象访问属性错误: {e}")

            # 检查是否有自定义标记
            if custom_labeling and emu_name in custom_labeling:
                dist = custom_labeling[emu_name]
                # 验证自定义分布
                if not isinstance(dist, np.ndarray):
                    dist = np.array(dist)

                # 验证分布数组长度
                if len(dist) != carbon_count + 1:
                    raise ValueError(f"EMU {emu_name} 的自定义分布长度({len(dist)})与碳原子数({carbon_count})不匹配，"
                                     f"应为 {carbon_count+1} (M+0到M+{carbon_count})")

                # 验证分布总和为1
                if not np.isclose(np.sum(dist), 1.0, atol=1e-5):
                    print(f"警告: EMU {emu_name} 的分布总和为 {np.sum(dist)}，已自动归一化")
                    dist = dist / np.sum(dist)

                # 验证分布值域在 [0,1] 之间
                if np.any(dist < 0) or np.any(dist > 1):
                    raise ValueError(f"EMU {emu_name} 的分布包含无效值 (需在0-1范围内): {dist}")

                input_emu_data_dict[emu_name] = dist.copy()
            else:
                # 检查是否有输入代谢物标记信息
                metabolite_base_name = emu.metabolite_name.split('_')[0]  # 处理像OAC_input这样的特殊情况
                if input_metabolite_dict and metabolite_base_name in input_metabolite_dict:
                    # 从输入代谢物标记信息中提取该代谢物的标记情况
                    labeling_info = input_metabolite_dict[metabolite_base_name]
                    # 对于每个碳位置，根据EMU的选定碳原子构建分布
                    carbon_atoms = [i for i, is_selected in enumerate(emu.selected_carbon_list) if is_selected]

                    # 改进：使用更高效的卷积计算
                    # 初始化分布为所有可能的同位素组合 - 采用零数组表示空状态
                    complete_distribution = np.zeros(carbon_count + 1)

                    # 遍历所有标记情况
                    for item in labeling_info:
                        abundance = item["abundance"]
                        ratio_list = item["ratio_list"]

                        # 创建当前标记情况的分布数组
                        current_dist = np.zeros(carbon_count + 1)

                        # 计算选定碳原子中标记的数量
                        labeled_count = sum(1 for i in carbon_atoms if i < len(ratio_list) and ratio_list[i] == 1)

                        # 设置相应位置的丰度值
                        current_dist[labeled_count] = abundance

                        # 将当前分布加入到总分布中
                        complete_distribution += current_dist

                    # 确保分布总和为1（处理舍入误差）
                    if not np.isclose(np.sum(complete_distribution), 1.0, atol=1e-5):
                        if np.sum(complete_distribution) > 0:  # 防止除以零
                            complete_distribution = complete_distribution / np.sum(complete_distribution)
                        else:
                            raise ValueError(f"EMU {emu_name} 的分布总和为零，无法归一化")

                    input_emu_data_dict[emu_name] = complete_distribution
                else:
                    # 使用自然丰度计算
                    if carbon_count in natural_dist_cache:
                        # 使用缓存的自然丰度分布
                        input_emu_data_dict[emu_name] = natural_dist_cache[carbon_count].copy()
                    else:
                        # 为单碳EMU计算分布
                        if carbon_count == 1:
                            dist = np.array([1 - natural_abundance, natural_abundance])
                        else:
                            # 使用卷积计算多碳EMU的分布
                            dist = np.array([1.0])  # 初始化为M+0的概率为1
                            single_carbon_dist = np.array([1 - natural_abundance, natural_abundance])

                            # 对每个碳位置进行卷积
                            for _ in range(carbon_count):
                                dist = np.convolve(dist, single_carbon_dist)[:carbon_count + 1]

                        # 缓存计算结果
                        natural_dist_cache[carbon_count] = dist.copy()
                        input_emu_data_dict[emu_name] = dist.copy()

    # 最终验证所有分布
    for emu_name, dist in input_emu_data_dict.items():
        if not validate_distribution(dist):
            raise ValueError(f"生成的EMU {emu_name} 分布无效: {dist}")

    return input_emu_data_dict


def validate_distribution(distribution):
    """
    验证同位素分布是否有效

    参数:
    distribution: 同位素分布数组 [M+0, M+1, ..., M+n]

    返回:
    bool: 分布是否有效 (非负数值且总和为1)
    """
    if not isinstance(distribution, np.ndarray):
        return False
    return np.all(distribution >= 0) and np.isclose(np.sum(distribution), 1.0, atol=1e-5)
