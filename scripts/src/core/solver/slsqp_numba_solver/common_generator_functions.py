


# This function constructs the EMU graph for MFA using pure Python.
# Parameters:
# emu_mid_equation_dict: Dictionary of EMU MID equations.
# flux_name_index_dict: Dictionary mapping flux names to indices.
# input_emu_data_dict: Dictionary of input EMU data.
# target_emu_name_list: List of target EMU names.
# all_target_metabolite_emu_name_list: List of all target metabolite EMU names.
# Returns: A tuple containing input EMU data list, operation list, complete EMU name index size dictionary, target EMU index dictionary, and all target EMU index dictionary.
# 一个元组，包含输入EMU数据列表、操作列表、完整的EMU名称索引大小字典、目标EMU索引字典和所有目标EMU索引字典。
# 这段代码定义了一个名为 emu_graph_constructor_optimized_pure_python 的函数，
# 它的主要任务是生成与代谢流分析相关的矩阵更新列表，并构造关于基本代谢单元（EMU）的图形和索引。
# 具体来说，它处理了代谢网络中EMU方程的构建与优化，并返回一组用于进一步计算和分析的数据结构。
# 传入形参：
# emu_mid_equation_dict：包含EMU方程的字典，每个代谢物的同位素分布方程。
# flux_name_index_dict：包含反应名称与反应索引的映射，用于查找每个反应的索引。
# input_emu_data_dict：输入的EMU数据字典，存储了与各代谢物相关的详细信息。
# target_emu_name_list：目标EMU名称列表，通常是需要分析的代谢物。
# all_target_metabolite_emu_name_list：所有目标代谢物的EMU名称列表，用于生成相关的矩阵方程。
# 总结：
# emu_graph_constructor_optimized_pure_python 函数的主要任务是构建代谢流分析中的EMU图，并生成用于计算代谢流的矩阵方程。
# 它通过遍历代谢物方程、生成矩阵更新列表、更新代谢物的索引，最终返回一系列用于后续代谢流优化和分析的数据结构。

def emu_graph_constructor_optimized_pure_python(
        emu_mid_equation_dict, flux_name_index_dict, input_emu_data_dict, target_emu_name_list,
        all_target_metabolite_emu_name_list):
    """
        构建用于代谢通量分析的EMU (基本代谢单位) 图结构。
        该函数生成数据结构，通过矩阵运算来计算代谢网络中的同位素分布。

        参数：
            emu_mid_equation_dict: 按碳原子数组织的EMU方程字典。
            flux_name_index_dict: 反应ID到通量数组中索引的映射。
            input_emu_data_dict: 包含同位素分布的输入EMU数据字典。
            target_emu_name_list: 目标EMU名称列表。
            all_target_metabolite_emu_name_list: 所有目标代谢物EMU名称列表。

        返回：
            input_emu_data_list: 输入EMU数据列表。
            operation_list: 矩阵更新操作列表。
            complete_emu_name_index_size_dict: EMU名称到 (索引,大小) 元组的映射。
            target_emu_index_dict: 目标EMU名称到其索引的映射。
            all_target_emu_index_dict: 所有目标代谢物EMU名称到其索引的映射。

        示例输出结构：
            all_emu_list = [np_array_for_EMU1, np_array_for_EMU2, ...]
            operation_list = [
                (carbon_num, matrix_a_dim, matrix_b_col, matrix_A_update_list, matrix_B_update_list,
                 matrix_Y_update_list),
                ...
            ]
            target_emu_index_dict = {target_emu_name: target_emu_index, ...}
    """

    # 这段文档字符串提供了该函数的输出数据结构的示意。
    # all_emu_list：所有EMU的NumPy数组表示。
    # operation_list：每个操作的更新列表，包含矩阵的维度、更新索引和更新的数据。
    # target_emu_index_dict：目标EMU名称到索引的映射。
    """

    all_emu_list = [np_array_for_EMU1, np_array_for_EMU2, ... ,]
    operation_list = [
        (carbon_num, matrix_a_dim, matrix_b_col, matrix_A_update_list, matrix_B_update_list,
        matrix_Y_update_list),
        ...]
    target_emu_index_dict = {
        target_emu_name: target_emu_index,
        ...
    }

    :param all_target_metabolite_emu_name_list:
    :param emu_mid_equation_dict:
    :param flux_name_index_dict:
    :param input_emu_data_dict:
    :param target_emu_name_list:
    :return:
    """

    # 内部函数：生成矩阵更新列表
    def generate_updating_list_of_matrix(_matrix_flux_location_dict, _flux_name_index_dict):
        """
                根据反应通量位置生成矩阵更新列表。

                参数：
                    _matrix_flux_location_dict: 存储矩阵位置和反应系数的字典。
                    _flux_name_index_dict: 反应ID到流量数组索引的映射。

                返回：
                    matrix_updates_list: 更新列表，每个更新包含行、列和通量更新信息。
        """
        matrix_updates_list = []
        # 遍历矩阵位置及其对应的反应数据
        for (row_index, col_index), reaction_dict in _matrix_flux_location_dict.items():
            flux_value_update_list = []
            # 收集该位置每个反应的通量索引和系数
            for reaction_id, coefficient in reaction_dict.items():
                flux_value_update_list.append((_flux_name_index_dict[reaction_id], coefficient))
            # 将更新元组添加到列表中
            matrix_updates_list.append((row_index, col_index, flux_value_update_list))
        return matrix_updates_list

    # 初始化用于处理的EMU方程字典
    emu_matrix_equation_carbon_num_dict = emu_mid_equation_dict
    # 字典，用于存储EMU名称到其索引和大小的映射
    complete_emu_name_index_size_dict = {}
    # 列表，用于存储输入EMU数据以供后续使用
    input_emu_data_list = []
    # 从输入EMU数据中填充数据列表和索引-大小字典
    for index, (emu_name, emu_data) in enumerate(input_emu_data_dict.items()):
        complete_emu_name_index_size_dict[emu_name] = (index, len(emu_data))
        input_emu_data_list.append(emu_data)

    # 列表，用于存储每个碳数层的矩阵更新操作
    operation_list = []
    # 处理每个碳数层及其关联的EMU方程
    for carbon_num, (
            this_layer_emu_dict_list, input_and_lower_layer_emu_dict_list, matrix_a_flux_location_dict,
            matrix_b_flux_location_dict, matrix_a_dim, matrix_b_col) in emu_matrix_equation_carbon_num_dict.items():
        # 跳过没有EMU的层
        if len(this_layer_emu_dict_list) == 0:
            continue

        # 使用内部函数生成矩阵A和B的更新列表
        matrix_a_update_list = generate_updating_list_of_matrix(matrix_a_flux_location_dict, flux_name_index_dict)
        matrix_b_update_list = generate_updating_list_of_matrix(matrix_b_flux_location_dict, flux_name_index_dict)

        # 生成矩阵Y的更新列表，表示对输入或低层EMU的依赖
        matrix_y_update_list = []
        for emu_name, input_and_lower_layer_emu in input_and_lower_layer_emu_dict_list.items():
            # 对于EMU的组合，获取每个组成EMU的索引
            if isinstance(input_and_lower_layer_emu, list):
                emu_index_list = [
                    complete_emu_name_index_size_dict[emu.full_name][0] for emu in input_and_lower_layer_emu]
             # 对于单个EMU，直接获取其索引
            else:
                emu_index_list = [complete_emu_name_index_size_dict[emu_name][0]]
            matrix_y_update_list.append(emu_index_list)

        # 为当前层的EMU分配索引
        this_layer_matrix_x_emu_index_list = []
        start_index = len(complete_emu_name_index_size_dict)
        for index, emu_obj in enumerate(this_layer_emu_dict_list):
            complete_index = start_index + index
            # 存储EMU的索引和大小（carbon_num + 1 反映同位素分布的长度）
            # current_complete_index = len(complete_emu_name_index_size_dict)
            # this_layer_matrix_x_emu_index_list.append(current_complete_index)
            complete_emu_name_index_size_dict[emu_obj.full_name] = (start_index + index, carbon_num + 1)
            this_layer_matrix_x_emu_index_list.append(complete_index)
        # 将此碳数层的操作详情添加到操作列表
        operation_list.append((
            carbon_num, matrix_a_dim, matrix_b_col, matrix_a_update_list, matrix_b_update_list, matrix_y_update_list,
            this_layer_matrix_x_emu_index_list))

    # 创建目标EMU名称到其索引的映射字典
    target_emu_index_dict = {
        target_emu_name: complete_emu_name_index_size_dict[target_emu_name][0]
        for target_emu_name in target_emu_name_list}


    # 创建所有目标代谢物EMU名称到其索引的映射字典
    all_target_emu_index_dict = {
        target_emu_name: complete_emu_name_index_size_dict[target_emu_name][0]
        for target_emu_name in all_target_metabolite_emu_name_list}

    # 返回所有构建的数据结构：
    # input_emu_data_list：包含所有输入EMU的数据。
    # operation_list：包含矩阵更新所需的操作数据。
    # complete_emu_name_index_size_dict：代谢物名称到索引的映射。
    # target_emu_index_dict 和 all_target_emu_index_dict：分别为目标EMUs的索引映射。
    return input_emu_data_list, operation_list, complete_emu_name_index_size_dict, target_emu_index_dict, \
        all_target_emu_index_dict