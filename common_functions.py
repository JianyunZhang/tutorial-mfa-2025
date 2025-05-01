import numpy as np
from scripts.src.core.common.config import CoreConstants

def generate_input_emu_data(input_metabolite_dict, input_emu_dict, natural_abundance=CoreConstants.natural_c13_ratio):
    """
    根据输入的同位素标记底物信息和预先计算的EMU列表生成输入EMU数据字典
    现在支持自然丰度的考虑

    Parameters:
    -----------
    input_metabolite_dict: Dict - 底物的同位素标记信息，指定每种底物的同位素标记模式及其比例
    input_emu_dict: Dict - 预先计算的EMU字典，包含需要计算MID的EMU列表
    natural_abundance: float - 自然界中碳13同位素的丰度比例，默认使用核心常量中定义的值

    Returns:
    --------
    input_emu_data_dict: Dict - 以 {emu_name: np.array([M+0, M+1, ..., M+n])} 格式的输入EMU数据字典
    """

    # 初始化结果字典
    input_emu_data_dict = {}

    # 处理输入代谢物字典
    if input_metabolite_dict is None:
        input_metabolite_dict = {}

    # 遍历input_emu_dict中的每个EMU
    for emu_name, emu_obj in input_emu_dict.items():
        # 获取代谢物名称和选定的碳原子列表
        metabolite_name = emu_obj.metabolite_name
        selected_carbon_list = emu_obj.selected_carbon_list

        # 计算选定碳原子的数量
        selected_carbons_count = sum(selected_carbon_list)

        # 初始化MID为全零数组
        mid = np.zeros(selected_carbons_count + 1)

        # 创建标记信息列表，无论代谢物是否在input_metabolite_dict中
        labeling_info_list = []

        # 如果代谢物在input_metabolite_dict中有标记信息，使用提供的信息
        if metabolite_name in input_metabolite_dict:
            labeling_info_list = input_metabolite_dict[metabolite_name]
        else:
            # 如果没有特定的标记信息，创建使用自然丰度的默认标记信息
            # 创建一个包含所有碳原子的自然丰度标记，abundance为1.0
            natural_ratio_list = [natural_abundance] * len(selected_carbon_list)
            labeling_info_list = [{'ratio_list': natural_ratio_list, 'abundance': 1.0}]

        # 统一处理所有标记信息
        for info in labeling_info_list:
            # 获取目前指定的标记反应底物的同位素标记MID丰度情况,以及所占比例
            ratio_list = info.get('ratio_list', [])
            abundance = info.get('abundance', 1.0)

            # 选择与selected_carbon_list位置对应的原子标记
            selected_ratio = []
            for i, sel in enumerate(selected_carbon_list):
                if sel == 1:  # 如果此碳原子被选择
                    if i < len(ratio_list):  # 确保索引在范围内
                        selected_ratio.append(ratio_list[i])
                    else:
                        # 超出范围的碳原子使用自然丰度
                        selected_ratio.append(natural_abundance)

            # 计算此标记模式下的MID（通过卷积）
            pattern_mid = np.array([1.0])  # 初始化为单位向量 [1]

            for ratio in selected_ratio:
                atom_mid = np.array([1.0-ratio, ratio])  # 单个碳原子的MID
                # 使用卷积计算整个EMU的MID
                pattern_mid = np.convolve(pattern_mid, atom_mid)

            # 将此标记模式的MID按比例添加到最终MID中
            mid += pattern_mid * abundance

        # 将卷积计算得到的该emu对应的MID添加到结果字典中
        input_emu_data_dict[emu_name] = mid

    return input_emu_data_dict


def prepare_mid_data_list(graph_results):
    """
    Prepares MID data list from graph_results for base_prediction_function
    """
    # Extract EMU index information dictionary
    emu_name_index_size_dict = graph_results[2]

    # Find maximum index
    max_index = max([index for _, (index, _) in emu_name_index_size_dict.items()], default=-1)

    # Create index to size mapping
    index_to_size = {index: size for _, (index, size) in emu_name_index_size_dict.items()}

    # Initialize list with appropriate distribution arrays
    complete_predicted_mid_data_list = []
    while len(complete_predicted_mid_data_list) <= max_index:
        size = index_to_size.get(len(complete_predicted_mid_data_list), 2)
        initial_values = np.ones(size) / size
        complete_predicted_mid_data_list.append(initial_values)

    # Copy input EMU data to correct positions
    input_emu_data_list = graph_results[0]
    for i, data in enumerate(input_emu_data_list):
        if i < len(complete_predicted_mid_data_list):
            complete_predicted_mid_data_list[i] = data

    return complete_predicted_mid_data_list


def print_mid_results(target_emu_name_list, target_indices, complete_predicted_mid_data_list):
    """
    打印 MID 预测结果，以小数形式显示（而非百分比）
    """
    for target_name, target_idx in zip(target_emu_name_list, target_indices):
        if target_idx < len(complete_predicted_mid_data_list):
            mid = complete_predicted_mid_data_list[target_idx]
            print(f"\n{target_name}的预测MID:")
            for i, value in enumerate(mid):
                print(f"M+{i}: {value:.4f}")    # 改为小数形式，增加精度到4位

