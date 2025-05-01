from ...common.packages import np, np_float_type
from ...common.functions import np_log_eps



# 该函数对给定的 EMU 数据列表和索引列表进行卷积操作，返回卷积结果。
# complete_emu_data_list: 完整的 EMU 数据列表。
# emu_index_list: EMU 索引列表。
# 返回值: 卷积结果数组。
def emu_index_list_conv(complete_emu_data_list, emu_index_list):
    # 添加参数检查
    if not emu_index_list:
        raise ValueError("emu_index_list cannot be empty")

    # 检查所有索引是否都在范围内
    for idx in emu_index_list:
        if idx >= len(complete_emu_data_list):
            raise ValueError(
                f"Index {idx} out of range for complete_emu_data_list with length {len(complete_emu_data_list)}")

    # 执行卷积操作
    result_array = complete_emu_data_list[emu_index_list[0]]
    for input_emu_index in emu_index_list[1:]:
        input_array = complete_emu_data_list[input_emu_index]
        result_array = np.convolve(result_array, input_array, mode='full')
    return result_array



# 这个函数 construct_and_update_matrix 用于根据提供的代谢流向量（flux_vector）和更新列表（update_list），构建并更新一个矩阵。
# 矩阵的每个元素会根据代谢流向量中的相应值以及更新系数进行计算。
# flux_vector：代谢流向量，表示每个反应的流速（通常是一个数组）。
# row_num：矩阵的行数。
# col_num：矩阵的列数。
# update_list：一个列表，包含了矩阵元素更新的信息。每个元素由 (row_index, col_index, flux_value_update_list) 组成，其中：
# row_index 和 col_index 指定矩阵的行和列位置。
# flux_value_update_list 包含了与该矩阵元素更新相关的代谢流索引和系数。
def construct_and_update_matrix(flux_vector, row_num, col_num, update_list):
    # matrix.fill(0)
    # 初始化矩阵
    # matrix：使用 np.zeros() 初始化一个全零矩阵，大小为 (row_num, col_num)，即指定的行数和列数。
    # dtype=np_float_type：指定矩阵的数据类型为 np_float_type，即浮动类型。这通常是 np.float32 或 np.float64，具体取决于代码的上下文。
    matrix = np.zeros((row_num, col_num), dtype = np_float_type)
    # 遍历并更新矩阵元素 update_list，该列表包含了矩阵元素的更新信息。
    # update_list：遍历提供的更新列表，
    # 每个元素是一个三元组 (row_index, col_index, flux_value_update_list)，表示要更新矩阵中的某个位置的值。
    for (row_index, col_index, flux_value_update_list) in update_list:
        # 计算矩阵元素的更新值
        element_update_value = 0
        # flux_value_update_list：每个 (flux_index, coefficient) 对表示如何根据代谢流向量中的某些值和系数来更新矩阵元素。
        # flux_index：指代在 flux_vector 中的索引位置，表示代谢流向量中的某一反应。
        # coefficient：与该反应流相乘的系数。
        for (flux_index, coefficient) in flux_value_update_list:
            # 更新计算：通过将 flux_vector[flux_index]（代谢流向量中的某个元素）与 coefficient（系数）相乘
            # 然后累加到 element_update_value，最终计算出该矩阵元素的更新值。
            element_update_value += flux_vector[flux_index] * coefficient
        # flux_index_array, coefficient_array
        # gathered_flux_vector = gather(flux_vector, flux_index_array)
        # element_update_value = sum(gathered_flux_vector * coefficient_array)
        # 更新矩阵中的元素
        # matrix[row_index, col_index]：根据 row_index 和 col_index 更新矩阵中的指定元素。
        # += element_update_value：将计算出的 element_update_value 加到该位置已有的值上，进行增量更新。
        matrix[row_index, col_index] += element_update_value
    # rol_col_pair_index_array = [(row_index, col_index), (row_index, col_index), ...]
    # element_update_value_array = [element_update_value, element_update_value, ...]
    # matrix = scatter_nd(rol_col_pair_index_array, element_update_value_array, (rol_num, col_num))
    # 返回矩阵：返回更新后的矩阵 matrix，矩阵的每个元素已经根据代谢流向量和系数进行了相应的更新。
    return matrix

# 该函数使用给定的操作列表和流量向量，预测并更新完整的 MID 数据列表。
# flux_vector: 流量向量。
# complete_predicted_mid_data_list: 完整的预测 MID 数据列表。
# operation_list: 操作列表，包含矩阵更新和 EMU 索引列表。
# 返回值: 无返回值，直接更新 complete_predicted_mid_data_list。
# 该函数 base_prediction_function 用于基于已知的反应流（flux）向量和一系列操作（操作列表），来计算代谢网络中代谢物的预测数据
# 它通过构建和更新矩阵，进行求解并最终更新预测的中间数据列表。以下是代码逐行的详细注释：
# flux_vector: 这是一个代谢流向量，表示每个反应的流速。
# complete_predicted_mid_data_list: 存储预测的中间代谢数据的列表。
# operation_list: 该列表包含一系列操作，每个操作包含了与代谢流相关的矩阵信息，后续会用于矩阵的更新和计算。
# 功能：base_prediction_function 主要用于基于输入的代谢流（flux_vector）和代谢网络的操作列表（operation_list）
# 通过构建矩阵方程，求解代谢物的预测数据并更新预测列表。
# 步骤：
# 构建和更新矩阵A和矩阵B：根据输入的代谢流和操作列表。
# 初始化并更新矩阵Y：根据给定的EMU索引和转换函数更新矩阵Y。
# 求解线性方程组：使用 np.linalg.solve() 求解代谢流的预测结果。
# 更新预测数据列表：将计算结果更新到 complete_predicted_mid_data_list 中。
# 该函数是代谢流预测模型的核心，基于线性方程组的求解机制，通过矩阵运算实现了代谢物预测的自动化计算和更新。
def base_prediction_function(flux_vector, complete_predicted_mid_data_list, operation_list):
    # 循环遍历操作列表，每次迭代都会更新一个矩阵并计算预测的中间数据。
    # 遍历 operation_list 中的每个操作，操作列表包含多个操作，每个操作对应于不同的碳数，以及与代谢流相关的矩阵信息（矩阵A、矩阵B、矩阵Y的更新列表等）。
    for (
            carbon_num, matrix_a_dim, matrix_b_col, matrix_a_update_list, matrix_b_update_list,
            matrix_y_update_list, this_layer_matrix_x_emu_index_list) in operation_list:

        # 矩阵A和矩阵B的构建与更新
        # matrix_a 和 matrix_b：调用 construct_and_update_matrix 函数来构建和更新矩阵A和矩阵B。
        # 这个函数利用传入的flux_vector（代谢流向量）、矩阵的维度以及更新列表来构建和更新矩阵。
        # matrix_a 的维度是 (matrix_a_dim, matrix_a_dim)，即为方阵。
        # matrix_b 的维度是 (matrix_a_dim, matrix_b_col)，矩阵B的列数由 matrix_b_col 给定。
        matrix_a = construct_and_update_matrix(flux_vector, matrix_a_dim, matrix_a_dim, matrix_a_update_list)
        matrix_b = construct_and_update_matrix(flux_vector, matrix_a_dim, matrix_b_col, matrix_b_update_list)

        # matrix_y.fill(0)
        matrix_y = np.zeros((matrix_b_col, carbon_num + 1))
        for row_index, emu_index_list in enumerate(matrix_y_update_list):
            # np.copyto(matrix_y[row_index], emu_index_list_conv(complete_predicted_mid_data_list, emu_index_list))
            matrix_y[row_index] = emu_index_list_conv(complete_predicted_mid_data_list, emu_index_list)

        matrix_x = np.linalg.solve(matrix_a, matrix_b @ matrix_y)
        for row, row_emu_index in zip(matrix_x, this_layer_matrix_x_emu_index_list):
            # complete_predicted_mid_data_list.append(row)
            # np.copyto(complete_predicted_mid_data_list[row_emu_index], row)
            complete_predicted_mid_data_list[row_emu_index] = row

    # return complete_predicted_mid_data_list