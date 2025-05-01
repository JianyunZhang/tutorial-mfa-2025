import numpy as np
from scripts.src.core.common.config import CoreConstants
from scripts.src.core.model.common_model_process_functions import model_preprocess
from scripts.src.core.model.emu_analyzer_functions import emu_equation_analyzer
from scripts.src.core.model.emu_analyzer_functions import emu_matrix_equation_generator
from scripts.src.core.model.emu_analyzer_functions import new_emu_matrix_equation_generator
from scripts.src.core.solver.slsqp_numba_solver.common_generator_functions import (
    emu_graph_constructor_optimized_pure_python
)
from scripts.src.core.solver.slsqp_numba_solver.python_initializer_optimizer_functions import (
    base_prediction_function
)
from common_functions import generate_input_emu_data, prepare_mid_data_list, print_mid_results


def run_mfa_pipeline(
        metabolite_reaction_dict,
        input_metabolite_name_set,
        complete_metabolite_dim_dict,
        target_metabolite_name_list,
        flux_name_index_dict,
        flux_vector,
        #custom_labeling=None,
        input_metabolite_dict=None,
        verbose=True
):
    """
    Generic MFA pipeline that runs the entire analysis process

    Parameters:
    -----------
    metabolite_reaction_dict: Dict - Metabolic network definition
    input_metabolite_name_set: Set - Set of input metabolite names
    complete_metabolite_dim_dict: Dict - Carbon atom counts for each metabolite
    target_metabolite_names: List - Names of target metabolites
    flux_name_index_dict: Dict - Maps flux names to indices
    flux_vector: ndarray - Flux values
    custom_labeling: Dict - Custom isotope labeling patterns
    natural_abundance: float - Natural abundance of 13C
    verbose: bool - Whether to print detailed output

    Returns:
    --------
    Dict containing all results including MID predictions
    """

    """运行代谢通量分析(MFA)流程的函数"""
    if verbose:
        print("启动代谢通量分析...")

    # Step 0: 预处理模型
    if verbose:
        print("\n0. 模型预处理...")

    # symmetrical_metabolite_set = {
    #     'Suc', 'Fum'
    # }
    symmetrical_metabolite_set = {}
    balance_excluded_metabolite_set = {
        'CO2'
    }

    # 调用model_preprocess对代谢模型进行预处理
    processed_model = model_preprocess(
        metabolite_reaction_dict,
        input_metabolite_name_set,
        complete_metabolite_dim_dict,
        flux_name_index_dict,
        symmetrical_metabolite_set=symmetrical_metabolite_set,
        emu_excluded_metabolite_set=balance_excluded_metabolite_set,
        balance_excluded_metabolite_set=balance_excluded_metabolite_set,
    )

    # 从预处理结果获取更新后的数据
    metabolite_reaction_dict = processed_model['metabolite_reaction_dict']
    input_metabolite_name_set = processed_model['input_metabolite_name_set']
    complete_metabolite_dim_dict = processed_model['complete_metabolite_dim_dict']
    flux_name_index_dict = processed_model['flux_name_index_dict']

    if verbose:
        print(f"  - 代谢物数量: {len(metabolite_reaction_dict)}")
        print(f"  - 输入代谢物: {len(input_metabolite_name_set)}")
        print(f"  - 目标代谢物: {len(target_metabolite_name_list)}")
        print(f"  - 反应通量数: {len(flux_name_index_dict)}")

    # Step 1: 构建EMU方程 EMU equation analysis
    if verbose:
        print("\n1. Analyzing EMU equations...")

    emu_equations, input_emu_dict = emu_equation_analyzer(
        metabolite_reaction_dict,
        target_metabolite_name_list,
        input_metabolite_name_set,
        complete_metabolite_dim_dict
    )

    if verbose:
    #     print(f"EMU equations generated: {len(emu_equations)}")
    #     # Count equations by type
    #     # Safely print information about EMU equations
    #     try:
    #         # Extract statistics about EMU equations
    #         metabolites_involved = set()
    #         for eq in emu_equations:
    #             if hasattr(eq, 'target_emu') and hasattr(eq.target_emu, 'metabolite_name'):
    #                 metabolites_involved.add(eq.target_emu.metabolite_name)
    #
    #         print(f"  - Total metabolites involved: {len(metabolites_involved)}")
            print(f"  - Total EMU equations: {len(emu_equations)}")
            print(f"  - EMU equations data: {emu_equations}")
            print(f"  - Input EMU dictionary: {input_emu_dict}")
    #     except Exception as e:
    #         print(f"  - Could not analyze EMU equation details: {e}")

    # Step 2: Generate input EMU data
    if verbose:
        print("\n2. Preparing input EMU data...")

    # 首先根据 emu_equations 函数输出的判断结果 input_emu_dict，得到该反应方程全过程的输入 emu
    # 然后根据给定的同位素标记的底物 input_metabolite_dict，计算获得对应的 input_emu_dict 每种代谢物所对应的同位素标记的含量
    # 反应结果以 input_emu_data_dict 返回
    input_emu_data_dict = generate_input_emu_data(
        input_metabolite_dict,
        input_emu_dict,
        natural_abundance=0
    )

    if verbose:
    #     print(f"Generated {len(input_emu_data_dict)} input EMU data items")
    #     print(f"Input EMU data dictionary: {input_emu_data_dict}" )
    #     # Count by metabolite
    #     metabolite_counts = {}
    #     for emu_name in input_emu_data_dict:
    #         metab = emu_name.split('__')[0]
    #         metabolite_counts[metab] = metabolite_counts.get(metab, 0) + 1
    #
    #     for metab, count in metabolite_counts.items():
    #         print(f"  - {metab}: {count} EMUs")
    #
    #     # Count custom vs natural abundance
    #     # custom_count = sum(1 for emu in input_emu_data_dict if custom_labeling and emu in custom_labeling)
    #     # print(f"  - Custom labeled EMUs: {custom_count}")
    #     # print(f"  - Natural abundance EMUs: {len(input_emu_data_dict) - custom_count}")
        pass

    # Step 3: Build matrix equations
    if verbose:
        print("\n3. Building EMU matrix equations...")

    emu_matrix_equations = emu_matrix_equation_generator(emu_equations, input_emu_data_dict)

    if verbose:
        print(f"Matrix equations generated: {len(emu_matrix_equations)}")
        print(f"Matrix equations data: {emu_matrix_equations}")
    #     # Safely analyze matrix sizes
    #     try:
    #         size_counts = {}
    #         # 检查 emu_matrix_equations 是否为字典
    #         if isinstance(emu_matrix_equations, dict):
    #             for carbon_num, matrix_data in emu_matrix_equations.items():
    #                 # 提取矩阵维度信息 (位置索引4是 matrix_a_dim)
    #                 if len(matrix_data) >= 5:
    #                     matrix_size = matrix_data[4]  # matrix_a_dim
    #                     size_counts[matrix_size] = size_counts.get(matrix_size, 0) + 1
    #                     print(f"  - Carbon {carbon_num}: Matrix size {matrix_size}x{matrix_size}")
    #         # 暂时不用这个
    #         # for matrix_eq in emu_matrix_equations:
    #         #     # Check if matrix_eq is a dictionary and has 'A' key
    #         #     if isinstance(matrix_eq, dict) and 'A' in matrix_eq:
    #         #         A_matrix = matrix_eq['A']
    #         #         if hasattr(A_matrix, 'shape'):
    #         #             size = A_matrix.shape[0]
    #         #             size_counts[size] = size_counts.get(size, 0) + 1
    #
    #         if size_counts:
    #             # 处理矩阵尺寸统计...
    #             print("  Matrix sizes:")
    #             for size, count in sorted(size_counts.items()):
    #                 print(f"  - {size}×{size}: {count} matrices")
    #
    #             max_size = max(size_counts.keys()) if size_counts else 0
    #             print(f"  - Largest matrix: {max_size}×{max_size}")
    #         else:
    #             print("  - No valid matrix dimensions found")
    #     except Exception as e:
    #         print(f"  - Error analyzing matrix equations: {e}")

    # Step 4: Build EMU graph
    if verbose:
        print("\n4. Constructing EMU graph...")

    # 为目标生成完整的EMU名称
    # Generate full EMU names for targets
    target_emu_name_list = []
    all_target_metabolite_emu_name_list = []

    for target in target_metabolite_name_list:
        carbon_count = complete_metabolite_dim_dict[target]
        full_emu_name = f"{target}__{'1' * carbon_count}"
        target_emu_name_list.append(full_emu_name)
        all_target_metabolite_emu_name_list.append(full_emu_name)

    graph_results = emu_graph_constructor_optimized_pure_python(
        emu_matrix_equations,
        flux_name_index_dict,
        input_emu_data_dict,
        target_emu_name_list,
        all_target_metabolite_emu_name_list
    )

    if verbose:
        print(f"EMU graph constructed successfully")
        # Extract statistics from graph_results
        input_emu_data_list = graph_results[0]
        operation_list = graph_results[1]
        emu_name_index_size_dict = graph_results[2]
        target_emu_index_dict = graph_results[3]

        print(f"  - Input EMUs: {len(input_emu_data_list)}")
        print(f"  - Input EMUs data: {input_emu_data_list}")
        print(f"  - Graph operations: {len(operation_list)}")
        print(f"  - Operation list: {operation_list}")
        print(f"  - Total EMU nodes: {len(emu_name_index_size_dict)}")
        print(f"  - Target EMUs: {len(target_emu_index_dict)}")


        # Count operation types
        op_types = {}
        for op in operation_list:
            op_type = op[0] if isinstance(op, tuple) else "unknown"
            op_types[op_type] = op_types.get(op_type, 0) + 1

        print("  Operation types:")
        for op_type, count in op_types.items():
            print(f"  - {op_type}: {count}")

    # Step 5: Run prediction
    if verbose:
        print("\n5. Running MFA prediction...")

    complete_predicted_mid_data_list = prepare_mid_data_list(graph_results)
    operation_list = graph_results[1]

    # Execute prediction
    base_prediction_function(
        flux_vector,
        complete_predicted_mid_data_list,
        operation_list
    )

    # Collect target indices
    target_indices = []
    target_emu_index_dict = graph_results[3]

    for target_name in target_emu_name_list:
        if target_name in target_emu_index_dict:
            target_indices.append(target_emu_index_dict[target_name])

    # Print results if verbose
    if verbose:
        print_mid_results(target_emu_name_list, target_indices, complete_predicted_mid_data_list)

    # Return results
    return {
        'emu_equations': emu_equations,
        'input_emu_dict': input_emu_dict,
        'input_emu_data_dict': input_emu_data_dict,
        'graph_results': graph_results,
        'target_emu_name_list': target_emu_name_list,
        'target_indices': target_indices,
        'predicted_mid_data': complete_predicted_mid_data_list,
    }

