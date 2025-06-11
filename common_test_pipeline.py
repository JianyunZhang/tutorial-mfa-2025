import numpy as np
import scipy.linalg # Added for QR decomposition and SVD
from scripts.src.core.common.config import CoreConstants
from scripts.src.core.model.common_model_process_functions import model_preprocess, compart_all_metabolites, metabolite_carbon_number_verification
from scripts.src.core.model.emu_analyzer_functions import emu_equation_analyzer
from scripts.src.core.model.emu_analyzer_functions import emu_matrix_equation_generator
from scripts.src.core.solver.slsqp_numba_solver.common_generator_functions import (
    emu_graph_constructor_optimized_pure_python
)
from scripts.src.core.solver.slsqp_numba_solver.python_initializer_optimizer_functions import (
    base_prediction_function
)
from scripts.src.core.solver.solver_construction_functions.common_construct_functions import (
    mixing_equation_constructor
)
from common_functions import generate_input_emu_data, prepare_mid_data_list, print_mid_results
from scripts.src.core.data.data_class import MIDData
from scipy.optimize import minimize
from scipy.optimize import linprog

from scripts.data.hct116_cultured_cell_line.data_metabolite_to_standard_name_dict import data_metabolite_to_standard_name_dict
from scripts.data.hct116_cultured_cell_line.specific_data_parameters import SpecificParameters
from scripts.src.core.solver.slsqp_numba_solver.python_initializer_optimizer_functions import (
    mix_mid_data_list,
    solver_objective_func
)
from scripts.src.core.model.common_model_process_functions import flux_balance_equation_generator
from scripts.data.common_functions import common_data_loader
from scripts.src.common.config import DataType
from typing import List, Tuple, Dict, Any, Optional, Union, Callable



def load_hct116_experimental_data(experiment_name="HCT116_WQ2101", condition="ctrl", index="average"):
    """
    加载HCT116细胞系的实验数据

    Parameters:
    -----------
    experiment_name: str - 实验名称，默认为"HCT116_WQ2101"
    condition: str - 实验条件，默认为"ctrl"
    index: int - 实验重复序号，默认为1

    Returns:
    --------
    experimental_mid_data_obj_dict: Dict - 实验测量的MID数据对象字典
    model_metabolite_to_standard_name_dict: Dict - 模型代谢物到标准名称的映射
    """
    from scripts.data.common_functions import common_data_loader
    from scripts.src.common.config import DataType
    from scripts.data.hct116_cultured_cell_line.data_metabolite_to_standard_name_dict import data_metabolite_to_standard_name_dict

    # 用common_data_loader加载数据，自动初始化并填充complete_dataset
    data_wrap_obj, keyword = common_data_loader(DataType.hct116_cultured_cell_line)

    # param_dict 里的 experiment_name 必须和 sheet_name_dict 里的 experiments 字段一致
    param_dict = {
        keyword.experiments: experiment_name,
        keyword.condition: condition,
        keyword.index: index
    }

    # 获取实验数据（MFAData对象）
    mfa_data_obj = data_wrap_obj.return_dataset(param_dict)

    # 返回MIDData对象字典和标准名映射
    return mfa_data_obj.experimental_mid_data_obj_dict, data_metabolite_to_standard_name_dict


# def diagnose_and_fix_mixing_equation_problem(
#         experimental_mid_data_obj_dict,
#         model_metabolite_to_standard_name_dict,
#         user_metabolite_to_standard_name_dict,
#         all_metabolite_name_list,
#         model_target_metabolite_compartment_dict,
#         verbose=True):
#     """
#     诊断并修复混合方程构建问题
#
#     Parameters:
#     -----------
#     experimental_mid_data_obj_dict: 实验数据字典
#     model_metabolite_to_standard_name_dict: 模型代谢物映射字典
#     user_metabolite_to_standard_name_dict: 用户定义的映射字典
#     all_metabolite_name_list: 模型中的所有代谢物名称
#     model_target_metabolite_compartment_dict: 模型代谢物组分字典
#     verbose: 是否显示详细信息
#
#     Returns:
#     --------
#     updated_model_metabolite_to_standard_name_dict: 更新后的映射字典
#     manual_mix_equation_dict: 手动构建的混合方程字典
#     """
#     # 初始化返回值
#     manual_mix_equation_dict = {}
#
#     if verbose:
#         print("\n开始诊断混合方程构建问题...")
#         print(f"实验数据中的代谢物: {len(experimental_mid_data_obj_dict)} 个")
#         print(
#             f"模型映射字典中的条目: {len(model_metabolite_to_standard_name_dict) if model_metabolite_to_standard_name_dict else 0} 个")
#         print(
#             f"用户定义的映射条目: {len(user_metabolite_to_standard_name_dict) if user_metabolite_to_standard_name_dict else 0} 个")
#
#     # 确保映射字典存在
#     if model_metabolite_to_standard_name_dict is None:
#         model_metabolite_to_standard_name_dict = {}
#
#     # 1. 合并用户定义的映射
#     if user_metabolite_to_standard_name_dict:
#         if verbose:
#             print("\n合并用户定义的映射...")
#
#         for exp_name, model_name in user_metabolite_to_standard_name_dict.items():
#             if exp_name in model_metabolite_to_standard_name_dict:
#                 if model_metabolite_to_standard_name_dict[exp_name] != model_name:
#                     if verbose:
#                         print(
#                             f"  更新映射: {exp_name} 从 {model_metabolite_to_standard_name_dict[exp_name]} 到 {model_name}")
#                     model_metabolite_to_standard_name_dict[exp_name] = model_name
#                 else:
#                     if verbose:
#                         print(f"  保持已有映射: {exp_name} → {model_name}")
#             else:
#                 model_metabolite_to_standard_name_dict[exp_name] = model_name
#                 if verbose:
#                     print(f"  添加新映射: {exp_name} → {model_name}")
#
#     # 2. 检查实验数据中的代谢物是否都有映射
#     unmapped_exp_metabolites = []
#     for exp_metabolite in experimental_mid_data_obj_dict.keys():
#         if exp_metabolite not in model_metabolite_to_standard_name_dict:
#             unmapped_exp_metabolites.append(exp_metabolite)
#
#     if unmapped_exp_metabolites and verbose:
#         print(f"\n发现 {len(unmapped_exp_metabolites)} 个未映射的实验代谢物:")
#         for metab in unmapped_exp_metabolites:
#             print(f"  - {metab}")
#
#     # 3. 尝试自动映射未映射的代谢物
#     if unmapped_exp_metabolites:
#         if verbose:
#             print("\n尝试自动映射未映射的代谢物...")
#
#         for exp_metabolite in unmapped_exp_metabolites:
#             # 使用名称相似度进行自动映射
#             best_match = None
#             highest_similarity = 0
#
#             exp_name_simple = exp_metabolite.lower().replace('-', '').replace('_', '').replace(' ', '')
#
#             for model_metabolite in all_metabolite_name_list:
#                 model_name_simple = model_metabolite.lower().replace('-', '').replace('_', '').replace(' ', '')
#
#                 # 完全包含关系
#                 if exp_name_simple in model_name_simple or model_name_simple in exp_name_simple:
#                     similarity = len(set(exp_name_simple) & set(model_name_simple)) / len(
#                         set(exp_name_simple) | set(model_name_simple))
#                     if similarity > highest_similarity:
#                         highest_similarity = similarity
#                         best_match = model_metabolite
#
#             if best_match and highest_similarity > 0.5:  # 设置相似度阈值
#                 model_metabolite_to_standard_name_dict[exp_metabolite] = best_match
#                 if verbose:
#                     print(f"  自动映射: {exp_metabolite} → {best_match} (相似度: {highest_similarity:.2f})")
#
#     # 4. 检查映射后的代谢物是否在模型中存在
#     invalid_mappings = []
#     for exp_metabolite, model_metabolite in model_metabolite_to_standard_name_dict.items():
#         if model_metabolite not in all_metabolite_name_list:
#             invalid_mappings.append((exp_metabolite, model_metabolite))
#
#     if invalid_mappings and verbose:
#         print(f"\n发现 {len(invalid_mappings)} 个无效映射 (模型中不存在目标代谢物):")
#         for exp_metab, model_metab in invalid_mappings:
#             print(f"  - {exp_metab} → {model_metab}")
#
#     # 5. 手动构建混合方程
#     if verbose:
#         print("\n准备手动构建混合方程...")
#
#     # 遍历所有有效映射
#     for exp_metabolite in experimental_mid_data_obj_dict.keys():
#         if exp_metabolite in model_metabolite_to_standard_name_dict:
#             model_metabolite = model_metabolite_to_standard_name_dict[exp_metabolite]
#             if model_metabolite in all_metabolite_name_list:
#                 # 查找此代谢物在模型组分中的位置
#                 compartments_found = []
#                 for comp in model_target_metabolite_compartment_dict.get('default', {}):
#                     if model_metabolite in model_target_metabolite_compartment_dict['default'][comp]:
#                         compartments_found.append(comp)
#
#                 if compartments_found:
#                     # 创建混合方程 (简单情况: 直接映射到模型代谢物)
#                     manual_mix_equation_dict[exp_metabolite] = {model_metabolite: 1.0}
#                     if verbose:
#                         print(
#                             f"  创建混合方程: {exp_metabolite} → {model_metabolite} (组分: {', '.join(compartments_found)})")
#                 else:
#                     if verbose:
#                         print(f"  警告: {model_metabolite} 在模型组分字典中未找到")
#
#     if verbose:
#         print(f"\n手动构建了 {len(manual_mix_equation_dict)} 个混合方程")
#         print(f"最终映射字典包含 {len(model_metabolite_to_standard_name_dict)} 个条目")
#
#     return model_metabolite_to_standard_name_dict, manual_mix_equation_dict

def convert_emu_dict_to_bal_dict(metabolite_reaction_dict_for_emu):
    """
    将 {metabolite: [Reaction, ...]} 结构的 EMU 字典，转换为 {metabolite: (dict, dict)} 结构的通量平衡字典。
    用于 flux_balance_equation_generator。
    """
    from collections import defaultdict
    # 这里用 DefaultDict((DefaultDict(0), DefaultDict(0))) 结构
    bal_dict = {}
    for metabolite in metabolite_reaction_dict_for_emu:
        bal_dict[metabolite] = (defaultdict(float), defaultdict(float))
    # 遍历所有反应，填充底物和产物
    for metabolite, reaction_list in metabolite_reaction_dict_for_emu.items():
        for reaction in reaction_list:
            # 处理底物
            for substrate_node in reaction.substrate_list:
                if substrate_node.name in bal_dict:
                    bal_dict[substrate_node.name][0][reaction.reaction_id] += substrate_node.coefficient
            # 处理产物
            for product_node in reaction.product_list:
                if product_node.name in bal_dict:
                    bal_dict[product_node.name][1][reaction.reaction_id] += product_node.coefficient
    return bal_dict


def convert_loss_operation_list_for_squared_loss(loss_operation_list):
    """
    将 loss_operation_list 中的 experimental_mid_data 字段从 MIDData 对象转换为 numpy 数组。
    """
    new_loss_operation_list = []
    for item in loss_operation_list:
        # 解包四元组
        mixed_metabolite_name, predicted_mid_index, experimental_mid_data, valid_index_array = item
        # 如果 experimental_mid_data 是 MIDData 对象，则取其 mid_list
        if hasattr(experimental_mid_data, 'mid_list'):
            exp_mid_array = np.array(experimental_mid_data.mid_list)
        else:
            exp_mid_array = np.array(experimental_mid_data)
        # 重新组装四元组
        new_loss_operation_list.append(
            (mixed_metabolite_name, predicted_mid_index, exp_mid_array, valid_index_array)
        )
    return new_loss_operation_list


# 在 SLSQP scipy.minimize 调用前添加参数验证
def validate_optimization_parameters(
        initial_flux_vector: np.ndarray,
        bounds: Optional[List[Tuple[Optional[float], Optional[float]]]],
        constraints: Optional[List[Dict[str, Any]]]
) -> Tuple[np.ndarray, Optional[List[Tuple[Optional[float], Optional[float]]]], Optional[List[Dict[str, Any]]]]:
    """
    验证并修正优化参数，确保参数符合优化算法的要求
    包含约束矩阵奇异性检测和处理

    Args:
        initial_flux_vector: 初始通量向量
        bounds: 边界条件
        constraints: 约束条件

    Returns:
        修正后的参数元组 (initial_flux_vector, bounds, constraints)
    """
    # 确保初始通量向量是numpy数组
    if not isinstance(initial_flux_vector, np.ndarray):
        initial_flux_vector = np.array(initial_flux_vector, dtype=float)

    # 检查初始向量是否包含无效值
    if np.any(np.isnan(initial_flux_vector)) or np.any(np.isinf(initial_flux_vector)):
        raise ValueError("初始通量向量包含NaN或Inf值")

    # 检查边界条件
    if bounds is not None:
        # 确保边界长度与初始向量匹配
        if len(bounds) != len(initial_flux_vector):
            raise ValueError(f"边界条件长度({len(bounds)})与初始向量长度({len(initial_flux_vector)})不匹配")

        # 确保每个边界都是有效的
        valid_bounds = []
        for i, (lb, ub) in enumerate(bounds):
            if lb is not None and ub is not None and lb > ub:
                raise ValueError(f"第{i}个变量的下界({lb})大于上界({ub})")

            # 处理None值
            if lb is None:
                lb = -np.inf
            if ub is None:
                ub = np.inf

            valid_bounds.append((lb, ub))
        bounds = valid_bounds

    # 检查约束条件并处理奇异矩阵问题
    if constraints is not None:
        processed_constraints = []
        
        for i, constraint in enumerate(constraints):
            # 检查约束类型
            if 'type' not in constraint:
                raise ValueError(f"第{i}个约束缺少'type'字段")

            # 检查约束函数
            if ('fun' not in constraint) or (not callable(constraint['fun'])):
                raise ValueError(f"第{i}个约束缺少有效的'fun'字段")

            # 对于等式约束，进行奇异矩阵检测和处理
            if constraint['type'] == 'eq' and 'jac' in constraint and callable(constraint['jac']):
                try:
                    # 获取雅可比矩阵（约束矩阵）
                    jac_func = constraint['jac']
                    A_eq = jac_func(initial_flux_vector)
                    
                    if isinstance(A_eq, np.ndarray) and A_eq.ndim == 2:
                        # 使用SVD检测矩阵秩
                        try:
                            U_svd, s_svd, Vh_svd = scipy.linalg.svd(A_eq)
                            tol_svd = np.finfo(float).eps * max(A_eq.shape) * (s_svd[0] if len(s_svd) > 0 else 1)
                            rank_A_eq = np.sum(s_svd > tol_svd)
                            
                            print(f"约束 {i}: 矩阵形状 {A_eq.shape}, 秩 {rank_A_eq}")
                            
                            # 如果矩阵是奇异的（秩缺陷）
                            if rank_A_eq < A_eq.shape[0]:
                                print(f"检测到奇异约束矩阵 {i}: 秩为 {rank_A_eq}, 行数为 {A_eq.shape[0]}")
                                
                                # 使用QR分解选择线性无关的约束
                                _Q_at, _R_at, P_at_cols = scipy.linalg.qr(A_eq.T, pivoting=True, mode='economic')
                                independent_row_indices = sorted(P_at_cols[:rank_A_eq])
                                
                                A_eq_processed = A_eq[independent_row_indices, :]
                                print(f"移除冗余约束后，矩阵形状从 {A_eq.shape} 变为 {A_eq_processed.shape}")
                                
                                # 获取约束函数的右侧向量
                                b_eq = constraint['fun'](initial_flux_vector) + A_eq @ initial_flux_vector
                                b_eq_processed = b_eq[independent_row_indices]
                                
                                # 创建新的约束函数
                                def new_constraint_fun(x):
                                    return A_eq_processed @ x - b_eq_processed
                                
                                def new_constraint_jac(x):
                                    return A_eq_processed
                                
                                # 更新约束
                                processed_constraint = constraint.copy()
                                processed_constraint['fun'] = new_constraint_fun
                                processed_constraint['jac'] = new_constraint_jac
                                processed_constraints.append(processed_constraint)
                                
                            else:
                                # 矩阵满秩，直接使用原约束
                                processed_constraints.append(constraint)
                                
                        except np.linalg.LinAlgError as la_err:
                            print(f"约束 {i} SVD分解失败: {la_err}, 使用原约束")
                            processed_constraints.append(constraint)
                            
                    else:
                        # 非矩阵约束，直接使用
                        processed_constraints.append(constraint)
                        
                except Exception as e:
                    print(f"处理约束 {i} 时出错: {e}, 使用原约束")
                    processed_constraints.append(constraint)
            else:
                # 非等式约束或无雅可比矩阵，直接使用
                processed_constraints.append(constraint)
                
        constraints = processed_constraints if processed_constraints else constraints

    return initial_flux_vector, bounds, constraints


def run_optimization(
        objective_function: Callable,
        initial_flux_vector: np.ndarray,
        bounds: Optional[List[Tuple[float, float]]] = None,
        constraints: Optional[List[Dict[str, Any]]] = None,
        method: str = 'SLSQP',
        verbose: bool = False
) -> Any:
    """
    执行优化算法，包含完善的错误处理和回退策略
    增强了约束矩阵奇异性检测和处理能力

    Args:
        objective_function: 目标函数
        initial_flux_vector: 初始通量向量
        bounds: 边界条件
        constraints: 约束条件
        method: 优化方法
        verbose: 是否显示详细信息

    Returns:
        优化结果
    """
    try:
        # 验证参数并处理奇异矩阵问题
        initial_flux_vector, bounds, constraints = validate_optimization_parameters(
            initial_flux_vector, bounds, constraints)

        # 确保初始值在边界范围内
        if bounds is not None:
            for i in range(len(initial_flux_vector)):
                lower = bounds[i][0]
                upper = bounds[i][1]

                if lower != -np.inf and initial_flux_vector[i] < lower:
                    initial_flux_vector[i] = lower + 1e-8  # 添加小偏移避免边界问题
                if upper != np.inf and initial_flux_vector[i] > upper:
                    initial_flux_vector[i] = upper - 1e-8  # 添加小偏移避免边界问题

        if verbose:
            print(f"优化参数验证通过：变量数量={len(initial_flux_vector)}, 边界数量={len(bounds) if bounds else 0}, 约束数量={len(constraints) if constraints else 0}")

        # 检查约束的条件数
        if constraints and method == 'SLSQP':
            for i, constraint in enumerate(constraints):
                if constraint['type'] == 'eq' and 'jac' in constraint:
                    try:
                        A_eq = constraint['jac'](initial_flux_vector)
                        if isinstance(A_eq, np.ndarray) and A_eq.ndim == 2:
                            cond_num = np.linalg.cond(A_eq)
                            if cond_num > 1e12:
                                print(f"警告: 约束 {i} 的条件数很大 ({cond_num:.2e}), 可能导致数值不稳定")
                    except:
                        pass

        # 运行优化
        if method == 'SLSQP':
            opt_result = minimize(
                objective_function,
                initial_flux_vector,
                method='SLSQP',
                bounds=bounds,
                constraints=constraints,
                options={
                    'disp': verbose,
                    'maxiter': 2000,
                    'ftol': 1e-10,
                    'eps': 1e-9
                }
            )
        else:
            opt_result = minimize(
                objective_function,
                initial_flux_vector,
                method=method,
                bounds=bounds,
                constraints=constraints if method.upper() in ['SLSQP', 'COBYLA', 'TRUST-CONSTR'] else None,
                options={
                    'disp': verbose,
                    'maxiter': 2000,
                    'ftol': 1e-10
                }
            )

        # 检查优化结果并提供详细信息
        if not opt_result.success and verbose:
            print(f"SLSQP优化失败: {opt_result.message}")
            if hasattr(opt_result, 'nit'):
                print(f"迭代次数: {opt_result.nit}")

        return opt_result

    except ValueError as e:
        print(f"优化参数验证失败或优化过程中出现错误: {e}")
        
        # 如果是奇异矩阵相关错误，尝试去除约束
        if "singular" in str(e).lower() or "Singular matrix" in str(e):
            print("检测到奇异矩阵错误，尝试使用仅边界约束的优化...")
            try:
                opt_result = minimize(
                    objective_function,
                    initial_flux_vector,
                    method='L-BFGS-B',
                    bounds=bounds,
                    options={
                        'disp': verbose,
                        'maxiter': 2000,
                        'ftol': 1e-10
                    }
                )
                return opt_result
            except Exception as e2:
                print(f"L-BFGS-B方法也失败: {e2}")
        
        # 其他错误，尝试L-BFGS-B方法
        print("尝试使用L-BFGS-B方法...")
        try:
            opt_result = minimize(
                objective_function,
                initial_flux_vector,
                method='L-BFGS-B',
                bounds=bounds,
                options={
                    'disp': verbose,
                    'maxiter': 2000,
                    'ftol': 1e-10
                }
            )
            return opt_result

        except Exception as e2:
            print(f"L-BFGS-B优化也失败了: {e2}")
            # 最后尝试无约束优化
            try:
                print("尝试使用Nelder-Mead方法（无约束）...")
                opt_result = minimize(
                    objective_function,
                    initial_flux_vector,
                    method='Nelder-Mead',
                    options={
                        'disp': verbose,
                        'maxiter': 3000
                    }
                )
                return opt_result
            except Exception as e3:
                print(f"所有优化方法都失败了: {e3}")
                # 创建一个失败的结果对象
                return type('Result', (), {
                    'success': False,
                    'message': f"无法完成优化: {e3}",
                    'x': initial_flux_vector,
                    'fun': float('inf')
                })()

    except Exception as e:
        print(f"优化过程中出现未知错误: {e}")

        # 检查是否是奇异矩阵错误
        if "singular" in str(e).lower() or "Singular matrix" in str(e):
            print("检测到奇异矩阵错误，尝试无约束优化...")
            try:
                opt_result = minimize(
                    objective_function,
                    initial_flux_vector,
                    method='Nelder-Mead',
                    options={
                        'disp': verbose,
                        'maxiter': 3000
                    }
                )
                return opt_result
            except Exception as e2:
                print(f"无约束优化也失败: {e2}")

        # 创建一个失败的结果对象
        return type('Result', (), {
            'success': False,
            'message': str(e),
            'x': initial_flux_vector,
            'fun': float('inf')
        })()


def run_mfa_pipeline(
        metabolite_reaction_dict,
        input_metabolite_name_set,
        complete_metabolite_dim_dict,
        target_metabolite_name_list,
        flux_name_index_dict,
        flux_vector,
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
    balance_excluded_metabolite_set = {'CO2'}

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
    metabolite_reaction_dict = processed_model['metabolite_reaction_dict'] # 这里对应flux balance constranct
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
        carbon_count = complete_metabolite_dim_dict[target] # 提取选中代谢物的碳原子数量 (从给定输入的 complete_metabolite_dim_dict 中获取)
        full_emu_name = f"{target}__{'1' * carbon_count}"   # 构造完整的目标代谢物 EMU 名称, 如 'Glu__11111'
        target_emu_name_list.append(full_emu_name)  # 把构造好的目标代谢物 EMU 名称添加进入空白列表 target_emu_name_list
        all_target_metabolite_emu_name_list.append(full_emu_name) # 把构造好的目标代谢物 EMU 名称添加进入空白列表 all_target_metabolite_emu_name_list
    # 元素格式：每个元素都是格式为 "{代谢物名}__{'1' * 碳原子数量}" 的字符串
    # 示例：如果代谢物 "F" 有 3 个碳原子，则其 EMU 名称为 "F__111"
    # target_emu_name_list: 用于最终结果输出和预测，标识需要计算质量同位素分布(MID)的目标 EMU
    # all_target_metabolite_emu_name_list: 用于 EMU 图构建过程，确保图中包含所有需要的目标代谢物 EMU 节点

    graph_results = emu_graph_constructor_optimized_pure_python(
        emu_matrix_equations,   # 先前一步生成的 EMU 矩阵方程 (包含了 A*X=B*Y 四个矩阵的数据, 以及矩阵维度信息)
        flux_name_index_dict,   # 反应名称-对应索引序号的映射字典, 这也是人为提前设置输入的参数
        input_emu_data_dict,    # 这是根据人为输入的同位素标记底物, 与先前 emu_equation_analyzer() 步骤算出的 input_emu_dict 相结合, 通过 generate_input_emu_data() 函数计算得出的 input_emu_data_dict 每个底物所占的同位素标记的比例
        target_emu_name_list,   # 目标代谢物 EMU 名称列表, 由刚刚一步生成
        all_target_metabolite_emu_name_list # 目标代谢物 EMU 名称列表, 由刚刚一步生成
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
    # operation_list 包含 碳原子数量, A矩阵维度, B矩阵列数
    # A矩阵更新列表, B矩阵更新列表 (和 emu_matrix_equation_generator 生成的 matrix_a_flux_location_dict 没什么区别)
    # Y矩阵更新列表 matrix_y_update_list (目前碳原子数量下所需的反应底物 EMU, 在列表中对标的位置, 以确定反应底物的碳原子系数)
    # X矩阵位置列表 this_layer_matrix_x_emu_index_list (存储了当前碳原子层数下, 每个底物在全部输入代谢物 EMU 列表中的索引位置)

    # Execute prediction
    base_prediction_function(
        flux_vector,    # 手工输入的反应名称对应编号
        complete_predicted_mid_data_list,   # 预测的同位素分布数据列表
        operation_list  # 上面构造的操作列表, 来源于 graph_results[1]
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

def run_mfa_pipeline_with_experimental_data(
        metabolite_reaction_dict,  # 代谢物反应字典
        input_metabolite_name_set, # 输入代谢物集合
        complete_metabolite_dim_dict, # 所有代谢物碳原子数量字典
        target_metabolite_name_list, # 目标代谢物名称列表
        flux_name_index_dict,      # 反应名称索引字典
        flux_vector,               # 反应系数向量

        user_metabolite_to_standard_name_dict, # 用户定义的代谢物到标准名称的映射
        flux_balance_reaction_dict,
        experimental_mid_data_obj_dict, # 实验测量的MID数据对象字典
        model_metabolite_to_standard_name_dict, # 模型代谢物到标准名称的映射
        specific_flux_range_dict=None, # 特定反应系数范围字典
        common_mix_ratio_range=(0.05, 0.95),  # 通用混合比例范围
        mix_ratio_multiplier=1.0,    # 混合比例乘数

        input_metabolite_dict=None,   # 输入代谢物字典
        verbose=True                  # 是否输出详细信息
):
    """
    集成实验数据的代谢通量分析(MFA)流程

    此函数扩展了基本的MFA流程，处理实验数据与模型数据之间的不匹配问题。
    它使用mixing_equation_constructor函数来构建混合方程，解决同位素质谱测量数据
    与代谢网络模型中代谢物节点可能不匹配的问题。

    Parameters:
    -----------
    metabolite_reaction_dict: Dict - 代谢网络模型定义
    input_metabolite_name_set: Set - 输入代谢物名称集合
    complete_metabolite_dim_dict: Dict - 每个代谢物的碳原子数量
    flux_name_index_dict: Dict - 通量名称到索引的映射
    flux_vector: ndarray - 流量值
    experimental_mid_data_obj_dict: Dict - 实验测量的MID数据对象字典
    model_metabolite_to_standard_name_dict: Dict - 模型代谢物到标准名称的映射
    specific_flux_range_dict: Dict - 特定通量范围字典，默认为None
    common_mix_ratio_range: tuple - 混合比例范围，默认为(0.0, 1.0)
    mix_ratio_multiplier: float - 混合比例乘数，默认为1.0
    input_metabolite_dict: Dict - 输入代谢物同位素标记字典，默认为None
    verbose: bool - 是否输出详细信息，默认为True

    Returns:
    --------
    Dict - 包含所有结果的字典，包括MID预测
    """
    """集成实验数据的代谢流量分析(MFA)流程，包含线性规划初始化和SLSQP优化"""

    # 步骤1: 运行基础MFA计算以获取模型结构
    if verbose:
        print("第1步: 运行基础MFA计算...")

    # 获取所有代谢物名称列表
    all_metabolite_name_list = list(complete_metabolite_dim_dict.keys())

    # 执行 target_metabolite_name_list
    # 提取目标代谢物名称，将实验数据中的代谢物名称进行提取并添加到 target_metabolite_name_list 中
    # 运行基于模拟代谢物的 MFA 分析步骤以前，需要确保 target_metabolite_name_list 中包含自行设定的模型及实验中存在的代谢物名称
    for exp_metabolite in experimental_mid_data_obj_dict:
        if exp_metabolite in user_metabolite_to_standard_name_dict:
            model_metabolite = user_metabolite_to_standard_name_dict[exp_metabolite]
            if model_metabolite in all_metabolite_name_list and model_metabolite not in target_metabolite_name_list:
                target_metabolite_name_list.append(model_metabolite)

    if not target_metabolite_name_list:
        print("警告: 没有找到匹配的目标代谢物!")
        target_metabolite_name_list = list(set(all_metabolite_name_list) - input_metabolite_name_set)

    # 运行基本的 MFA 分析流程获取模型结构
    base_results = run_mfa_pipeline(
        metabolite_reaction_dict=metabolite_reaction_dict,
        input_metabolite_name_set=input_metabolite_name_set,
        complete_metabolite_dim_dict=complete_metabolite_dim_dict,
        target_metabolite_name_list=target_metabolite_name_list,
        flux_name_index_dict=flux_name_index_dict,
        flux_vector=flux_vector,
        input_metabolite_dict=input_metabolite_dict,
        verbose=True
    )

    # 步骤2: 处理代谢物映射关系
    if verbose:
        print("\n第2步: 处理代谢物映射关系...")

    # 处理细胞组分和标准名称映射 (细胞质和线粒体)
    model_compartment_set = {'c', 'm'}

    # 从模型中提取所有代谢物信息，生成线粒体和细胞质里面的代谢物
    complete_tissue_compartment_metabolite_dict, metabolite_bare_metabolite_name_dict = compart_all_metabolites(
        all_metabolite_name_list,
        model_compartment_set
    )

    # 验证代谢物的碳原子数一致性
    bare_metabolite_dim_dict = metabolite_carbon_number_verification(
        complete_metabolite_dim_dict,
        metabolite_bare_metabolite_name_dict
    )

    # 创建模型目标代谢物组分字典
    # 使用'default'作为组织名（单组织模型）
    model_target_metabolite_compartment_dict = {'default': {}}

    # 为每个细胞组分初始化空集合
    for comp in model_compartment_set:
        model_target_metabolite_compartment_dict['default'][comp] = set()

    # 将代谢物添加到相应的组分集合
    for metabolite in target_metabolite_name_list:
        if metabolite in complete_tissue_compartment_metabolite_dict:
            # 将代谢物添加到所属的组分集合
            for comp in complete_tissue_compartment_metabolite_dict[metabolite]:
                model_target_metabolite_compartment_dict['default'][comp].add(metabolite)
        else:
            # 如果没有组分信息，默认为细胞质
            model_target_metabolite_compartment_dict['default']['c'].add(metabolite)

    if verbose:
        print(f"\n模型组分结构:")
        for tissue, compartments in model_target_metabolite_compartment_dict.items():
            print(f"  组织 {tissue}:")
            for comp, metabolites in compartments.items():
                print(f"    细胞器 {comp}: {len(metabolites)} 个代谢物")
                if len(metabolites) > 0:
                    print(f"      代谢物: {list(metabolites)[:5]}{'...' if len(metabolites) > 5 else ''}")



    # 步骤3: 构建混合方程
    if verbose:
        print("\n第3步: 构建混合方程...")

    # 获取基本的 MFA 分析流程得到的模型结构预测结果
    graph_results = base_results['graph_results']
    target_emu_name_list = base_results['target_emu_name_list']
    target_indices = base_results['target_indices']
    complete_predicted_mid_data_list = base_results['predicted_mid_data']
    emu_name_index_dict = base_results['graph_results'][3]

    # 在run_mfa_pipeline_with_experimental_data函数中
    print("\n开始构建混合方程...")
    print(f"  - 实验代谢物数量: {len(experimental_mid_data_obj_dict)}")
    print(f"  - 模型到标准名称映射数量: {len(model_metabolite_to_standard_name_dict)}")

    metabolite_bare_metabolite_name_dict = dict(zip(model_metabolite_to_standard_name_dict.values(), model_metabolite_to_standard_name_dict.keys()))

    # 构建混合方程 - 处理模型代谢物与实验代谢物之间的映射关系
    # 这里面要生成 target_metabolite_name_list,然后用在run_mfa_pipeline里面
    nested_mix_equation_dict, mix_ratio_name_index_dict, mix_ratio_balance_list, updated_specific_flux_range_dict = mixing_equation_constructor(
        experimental_mid_data_obj_dict=experimental_mid_data_obj_dict,
        model_target_metabolite_compartment_dict=model_target_metabolite_compartment_dict,
        model_metabolite_to_standard_name_dict=model_metabolite_to_standard_name_dict,
        metabolite_bare_metabolite_name_dict=metabolite_bare_metabolite_name_dict,
        specific_flux_range_dict=specific_flux_range_dict,
        common_mix_ratio_range=common_mix_ratio_range,
        mix_ratio_multiplier=mix_ratio_multiplier
        # list_of_case_name=None # None(默认值) 适用于单组实验数据情况，函数直接处理传入的 experimental_mid_data_obj_dict
    )

    # 调用函数后
    print(f"\n混合方程构建结果:")
    print(f"  - 混合方程数量: {len(nested_mix_equation_dict)}")
    print(f"  - 混合比率数量: {len(mix_ratio_name_index_dict)}")
    print(f"  - 平衡约束数量: {len(mix_ratio_balance_list)}")

    # 步骤4: 使用线性规划生成初始反应流量的系数值
    if verbose:
        print("\n第4步: 使用线性规划生成初始通量值...")

    # 构建质量平衡约束条件 s·v = 0
    # 设置反应系数约束
    flux_bounds = []
    for flux_name, idx in flux_name_index_dict.items():
        if specific_flux_range_dict and flux_name in specific_flux_range_dict:
            min_val, max_val = specific_flux_range_dict[flux_name]
        else:
            min_val, max_val = 0.1, 1000.0  # 约束默认范围
        flux_bounds.append((min_val, max_val))

    # 构建线性规划目标函数（使用随机目标函数，避免解的偏向性）,这个是要随机的，每次都要随机生成，有正有负，所以要-0.5，建议用 np random去生成
    c = np.random.random(len(flux_name_index_dict)) - 0.5

    # 从代谢网络构建约束条件
    # 获取代谢反应平衡矩阵（如果有的话），这里需要用到完整的代谢物反应字典和通量名称索引字典，不能只用每个代谢物节点流入的反应
    # processed_model = model_preprocess(
    #     metabolite_reaction_dict=flux_balance_reaction_dict,
    #     input_metabolite_name_set=input_metabolite_name_set,
    #     complete_metabolite_dim_dict=complete_metabolite_dim_dict,
    #     flux_name_index_dict=flux_name_index_dict,
    #     #reaction_list=,
    #     #symmetrical_metabolite_set={''},
    #     #added_input_metabolite_set=,
    #     emu_excluded_metabolite_set={'CO2'},
    #     balance_excluded_metabolite_set={'CO2'}
    #     #target_metabolite_name_list=,
    #     #composite_reaction_list=
    # )


    # metabolite_reaction_dict_for_bal = convert_emu_dict_to_bal_dict(processed_model['metabolite_reaction_dict'])
    # 然后 bal_dict 可直接传给 flux_balance_equation_generator
    
    # 尝试生成反应流量平衡约束条件
    # try:
    # TODO:
    # del metabolite_reaction_dict_for_bal['Glu']

    # 这里需要重点修改
    # metabolite: (flux_dict_as_substrate, flux_dict_as_product)
    metabolite_reaction_dict_for_bal = {
        'OAC': ({'v1': 1.0, 'v7': 1.0}, {'v8': 1.0, 'v6': 1.0}),
        'Cit': ({'v2': 1.0}, {'v1': 1.0}),
        'AKG': ({'v3': 1.0, 'v4': 1.0}, {'v2': 1.0}),
        'Suc': ({'v5': 1.0}, {'v4': 1.0}),
        'Fum': ({'v6': 1.0}, {'v5': 1.0, 'v7': 1.0})
    }

    flux_balance_matrix, flux_balance_right_side_vector = flux_balance_equation_generator(
        metabolite_reaction_dict=metabolite_reaction_dict_for_bal,  # 代谢物反应字典, 这里必须是processed_model的
        composite_reaction_dict={},    # 在这个范例一般来说没有复合反应
        flux_name_index_dict=flux_name_index_dict
    )
    A_eq = flux_balance_matrix
    b_eq = flux_balance_right_side_vector

    # 运行线性规划
    result = linprog(
        c=c,
        A_eq=A_eq if A_eq is not None and len(A_eq) > 0 else None,
        b_eq=b_eq if b_eq is not None and len(b_eq) > 0 else None,
        bounds=flux_bounds,
        method='simplex'
    )

    if result.success:
        initial_flux_vector = result.x
        if verbose:
            print("线性规划生成初始反应系数成功")
            print("初始反应流量:")
            for flux_name, idx in flux_name_index_dict.items():
                print(f"  {flux_name}: {initial_flux_vector[idx]:.4f}")
    else:
        if verbose:
            print(f"线性规划求解失败: {result.message}，使用输入的反应系数作为初始值")
        initial_flux_vector = flux_vector


    # 步骤5: 准备SLSQP优化
    if verbose:
        print("\n第5步: 准备SLSQP优化...")

    # 调试信息
    print(f"混合方程字典内容：{nested_mix_equation_dict}")
    print(f"代谢物维度字典中代谢物数量: {len(complete_metabolite_dim_dict)}")
    print(f"流量名称索引字典中通量数量: {len(flux_name_index_dict)}")
    print(f"EMU名称索引字典中EMU数量: {len(emu_name_index_dict)}")




    # 准备混合操作列表
    mix_operation_list = []
    for mixed_metabolite, nested_mix_equation in nested_mix_equation_dict.items():
        mix_op = []
        print(f"处理混合代谢物: {mixed_metabolite}, 混合方程: {nested_mix_equation}")

        # 处理简化的直接映射（字符串到字符串）
        if isinstance(nested_mix_equation, str):
            model_metabolite = nested_mix_equation
            print(f"  检查简化映射组分: {model_metabolite}")
            
            # 检查是否是流量名称
            flux_index = flux_name_index_dict.get(model_metabolite)
            if flux_index is None:
                print(f"  - {model_metabolite} 不在流量索引字典中")
                continue

            # 检查碳原子数
            carbon_num = complete_metabolite_dim_dict.get(model_metabolite)
            if carbon_num is None:
                print(f"  - {model_metabolite} 不在代谢物维度字典中")
                continue

            # 构造EMU名称并查找索引
            emu_name = f"{model_metabolite}__{'1' * carbon_num}"
            print(f"  - EMU名称: {emu_name}")

            emu_index = emu_name_index_dict.get(emu_name)
            if emu_index is None:
                print(f"  - {emu_name} 不在EMU索引字典中")
                # 尝试查找类似的EMU
                similar_emus = [name for name in emu_name_index_dict.keys() if model_metabolite in name]
                if similar_emus:
                    print(f"  - 找到类似EMU: {similar_emus}")
                    emu_name = similar_emus[0]
                    emu_index = emu_name_index_dict.get(emu_name)
                if emu_index is None:
                    continue

            print(f"  - 添加: ({flux_index}, {emu_index})")
            mix_op.append((flux_index, emu_index))
            
        # 处理复杂的嵌套字典结构
        elif isinstance(nested_mix_equation, dict):
            for model_metabolite, coeff in nested_mix_equation.items():
                print(f"  检查嵌套组分: {model_metabolite}, 系数: {coeff}")
                if coeff > 0:
                    # 检查是否是流量名称
                    flux_index = flux_name_index_dict.get(model_metabolite)
                    if flux_index is None:
                        print(f"  - {model_metabolite} 不在流量索引字典中")
                        continue

                    # 检查碳原子数
                    carbon_num = complete_metabolite_dim_dict.get(model_metabolite)
                    if carbon_num is None:
                        print(f"  - {model_metabolite} 不在代谢物维度字典中")
                        continue

                    # 构造EMU名称并查找索引
                    emu_name = f"{model_metabolite}__{'1' * carbon_num}"
                    print(f"  - EMU名称: {emu_name}")

                    emu_index = emu_name_index_dict.get(emu_name)
                    if emu_index is None:
                        print(f"  - {emu_name} 不在EMU索引字典中")
                        # 尝试查找类似的EMU
                        similar_emus = [name for name in emu_name_index_dict.keys() if model_metabolite in name]
                        if similar_emus:
                            print(f"  - 找到类似EMU: {similar_emus}")
                            emu_name = similar_emus[0]
                            emu_index = emu_name_index_dict.get(emu_name)
                        if emu_index is None:
                            continue

                    print(f"  - 添加: ({flux_index}, {emu_index})")
                    mix_op.append((flux_index, emu_index))

        if mix_op:
            print(f"  成功添加混合操作: {mix_op}")
            # 对于简化的映射，我们创建一个特殊的混合操作格式
            if isinstance(nested_mix_equation, str):
                # 简化映射不需要混合比例索引，直接使用 -1 或者省略
                mix_operation_list.append((mixed_metabolite, mix_op))
            else:
                # 复杂映射需要混合比例索引
                mix_ratio_index = mix_ratio_name_index_dict.get(mixed_metabolite, 0)
                mix_operation_list.append((mixed_metabolite, mix_ratio_index, mix_op))
        else:
            print(f"  警告: {mixed_metabolite}的混合操作列表为空")

    print(f"最终混合操作列表长度: {len(mix_operation_list)}")



    # 步骤5中创建损失操作列表时，需要完整构建四元组
    loss_operation_list = []

    # 处理混合操作情况
    if len(mix_operation_list)==0:
        for mix_item in mix_operation_list:
            # 处理不同格式的混合操作
            if len(mix_item) == 2:
                # 简化格式: (mixed_metabolite, mix_op)
                mixed_metabolite, mix_op = mix_item
                mix_ratio_index = None  # 简化映射不需要混合比例索引
            elif len(mix_item) == 3:
                # 完整格式: (mixed_metabolite, mix_ratio_index, mix_op)
                mixed_metabolite, mix_ratio_index, mix_op = mix_item
            else:
                print(f"警告: 无法识别的混合操作格式: {mix_item}")
                continue
            
            if mixed_metabolite in experimental_mid_data_obj_dict:
                # 获取实验MID数据
                experimental_mid_data = experimental_mid_data_obj_dict[mixed_metabolite]

                # 确定预测索引
                if mix_ratio_index is not None:
                    # 有混合比例的情况，预测索引是混合后的索引
                    predicted_mid_index = mix_ratio_index
                else:
                    # 简化映射的情况，直接使用EMU索引
                    if mix_op and len(mix_op) > 0:
                        _, emu_index = mix_op[0]  # 取第一个EMU索引
                        predicted_mid_index = emu_index
                    else:
                        print(f"警告: {mixed_metabolite} 的混合操作为空")
                        continue

                # 创建有效索引数组 - 这里使用实验数据中的所有位置
                # 假设experimental_mid_data.mid_list包含MID数据
                if hasattr(experimental_mid_data, 'mid_list') and experimental_mid_data.mid_list is not None:
                    valid_index_array = np.arange(len(experimental_mid_data.mid_list), dtype=np.int32)
                elif hasattr(experimental_mid_data, 'data_vector') and experimental_mid_data.data_vector is not None:
                    valid_index_array = np.arange(len(experimental_mid_data.data_vector), dtype=np.int32)
                else:
                    # 如果无法获取MID列表长度，使用默认值
                    valid_index_array = np.array([0, 1, 2, 3, 4], dtype=np.int32)  # 默认使用前5个质量同位素

                # 添加完整的四元组
                loss_operation_list.append(
                    [mixed_metabolite, predicted_mid_index, experimental_mid_data, valid_index_array])
                if verbose:
                    print(
                        f"添加完整损失操作: {mixed_metabolite}, 预测索引: {predicted_mid_index}, 有效索引数组: {valid_index_array}")
    else:
        # 直接使用模型代谢物预测结果的情况
        print("没有混合操作，尝试直接使用模型预测结果...")
        for exp_name, model_name in matched_metabolites:
            if exp_name in experimental_mid_data_obj_dict:
                # 获取实验MID数据
                experimental_mid_data = experimental_mid_data_obj_dict[exp_name]

                # 获取模型代谢物EMU名称
                carbon_num = complete_metabolite_dim_dict.get(model_name)
                if carbon_num is not None:
                    emu_name = f"{model_name}__{'1' * carbon_num}"

                    # 获取EMU索引
                    if emu_name in emu_name_index_dict:
                        # 处理可能的列表或单值情况
                        predicted_mid_index = emu_name_index_dict[emu_name]
                        if isinstance(predicted_mid_index, (list, tuple)):
                            predicted_mid_index = predicted_mid_index[0]

                        # 创建有效索引数组
                        if hasattr(experimental_mid_data, 'mid_list') and experimental_mid_data.mid_list is not None:
                            valid_index_array = np.arange(len(experimental_mid_data.mid_list), dtype=np.int32)
                        elif hasattr(experimental_mid_data, 'data_vector') and experimental_mid_data.data_vector is not None:
                            valid_index_array = np.arange(len(experimental_mid_data.data_vector), dtype=np.int32)
                        else:
                            valid_index_array = np.array([0, 1, 2, 3, 4], dtype=np.int32)

                        # 添加完整的四元组
                        loss_operation_list.append(
                            [exp_name, predicted_mid_index, experimental_mid_data, valid_index_array])
                        if verbose:
                            print(
                                f"添加完整损失操作: {exp_name}, 预测索引: {predicted_mid_index}, 有效索引数组: {valid_index_array}")

    if verbose:
        print(f"\n操作列表构建结果:")
        print(f"  - 混合操作数量: {len(mix_operation_list)}")
        print(f"  - 损失操作数量: {len(loss_operation_list)}")

    # 步骤6: 运行SLSQP优化
    if verbose:
        print("\n第6步: 运行SLSQP优化...")

    # 获取基础操作列表
    operation_list = base_results['graph_results'][1]
    emu_name_index_dict = base_results['graph_results'][3]

    # 将 loss_operation_list 中的 experimental_mid_data 字段从 MIDData 对象转换为 numpy 数组。
    # 满足 solver_objective_func 的要求
    new_loss_operation_list=convert_loss_operation_list_for_squared_loss(loss_operation_list)

    # 定义目标函数
    def objective_function(x):
        # 创建当前预测数据的副本
        curr_predicted_mid_data_list = [arr.copy() for arr in complete_predicted_mid_data_list]

        # 计算损失
        loss = solver_objective_func(
            flux_vector=x,
            complete_predicted_mid_data_list=curr_predicted_mid_data_list,
            operation_list=operation_list,
            mix_operation_list=mix_operation_list,
            mix_ratio_multiplier=mix_ratio_multiplier,
            loss_operation_list=new_loss_operation_list,
            loss_code=1, # 使用均方差损失
            optimal_cross_entropy=0.0# 熵损失参数
        )
        return loss

    initial_obj = objective_function(initial_flux_vector)
    print(f"初始目标函数值: {initial_obj:.4f}")
    opt_result = run_optimization(
        objective_function=objective_function,
        initial_flux_vector=initial_flux_vector,
        bounds=flux_bounds,
        constraints=[{
            'type': 'eq',
            'fun': lambda x: A_eq @ x - b_eq,
            'jac': lambda x: A_eq
        }],
        method='SLSQP',  # 默认使用SLSQP方法
        verbose=True  # 启动详细输出
    )
    opt_result.fun


    # 设置通量约束边界，这个可以用初始解的bonds
    bounds = []
    for flux_name, idx in flux_name_index_dict.items():
        if specific_flux_range_dict and flux_name in specific_flux_range_dict:
            min_val, max_val = specific_flux_range_dict[flux_name]
        else:
            min_val, max_val = 0.1, 1000.0
        bounds.append((min_val, max_val))

    # 构建线性等式约束条件
    # A_eq @ x - b_eq
    # 还要提供 equality constraints ( jacob constrains ) 默认等于 A_eq

    # 定义线性等式约束 (A_eq * x = b_eq)
    A_eq = flux_balance_matrix
    b_eq = flux_balance_right_side_vector

    # 定义线性不等式约束 (如果需要的话，格式为 A_ub * x <= b_ub)
    # 示例：假设有一些流量总和约束
    # A_ub = np.array([[1, 1, 0, ...], [0, 1, 1, ...]])  # 根据实际需求定义
    # b_ub = np.array([upper_bound1, upper_bound2, ...])  # 对应的上界

    # 创建约束条件列表
    constraints = []
    if (isinstance(A_eq, np.ndarray) and isinstance(b_eq, np.ndarray) and
       A_eq.ndim == 2 and b_eq.ndim == 1 and
       A_eq.shape[0] == b_eq.shape[0] and
       A_eq.shape[0] > 0 and
       (A_eq.shape[1] == len(initial_flux_vector) if len(initial_flux_vector) > 0 else A_eq.shape[1] == 0)):
        constraints.append({
            'type': 'eq',
            'fun': lambda x: A_eq @ x - b_eq,
            'jac': lambda x: A_eq
        })
        if verbose:
            print(f"添加了 {A_eq.shape[0]} 个线性等式约束 (Ax-b=0). A_eq shape: {A_eq.shape}, b_eq shape: {b_eq.shape}, initial_flux_vector length: {len(initial_flux_vector)}")
    elif verbose:
        print(f"未添加线性等式约束。A_eq is None: {A_eq is None}, b_eq is None: {b_eq is None}")
        if isinstance(A_eq, np.ndarray): print(f"  A_eq.ndim: {A_eq.ndim}, A_eq.shape: {A_eq.shape}")
        if isinstance(b_eq, np.ndarray): print(f"  b_eq.ndim: {b_eq.ndim}, b_eq.shape: {b_eq.shape}")
        if isinstance(initial_flux_vector, np.ndarray): print(f"  initial_flux_vector length: {len(initial_flux_vector)}")

    # 步骤8: 运行优化
    if verbose:
        print("\n开始运行SLSQP优化流程...")

    # 验证参数 (仍然保留这一步，确保参数格式正确)
    # initial_flux_vector, bounds, constraints = validate_optimization_parameters(
    #     initial_flux_vector, bounds, constraints)

    # 检查 initial_flux_vector 是否为空
    if len(initial_flux_vector) == 0:
        if verbose:
            print("优化中止：初始通量向量为空 (无优化变量)。")
        # 创建一个模拟的失败优化结果对象
        opt_result = type('OptimizationResult', (), {
            'x': np.array([]),
            'fun': float('inf'),
            'success': False,
            'message': "无优化变量，优化未执行。"
        })()
    else:
        # initial_flux_vector, bounds, constraints = validate_optimization_parameters(
        # initial_flux_vector, bounds, constraints) # Validation is inside run_optimization
        opt_result = run_optimization(
            objective_function=objective_function,
            initial_flux_vector=initial_flux_vector,
            bounds=bounds,
            constraints=constraints,
            method='SLSQP',  # 默认使用SLSQP方法
            verbose=True # 启动详细输出
        )

    # 记录优化结果
    optimized_flux_vector = opt_result.x
    objective_value = opt_result.fun
    optimization_success = opt_result.success
    optimization_message = opt_result.message if hasattr(opt_result, 'message') else "无优化信息"

    if verbose:
        print(f"优化结果: {'成功' if optimization_success else '失败'}")
        print(f"目标函数值: {objective_value}")
        print(f"优化信息: {optimization_message}")


    if opt_result.success:
        optimized_flux_vector = opt_result.x
        if verbose:
            print(f"优化成功! 最终损失值: {opt_result.fun:.6f}")
            print("\n优化后的通量值:")
            for flux_name, idx in flux_name_index_dict.items():
                print(f"  {flux_name}: {optimized_flux_vector[idx]:.4f}")
    else:
        if verbose:
            print(f"优化失败: {opt_result.message}")
        optimized_flux_vector = initial_flux_vector

    # 步骤7: 使用优化后的流量系数值重新计算预测结果
    if verbose:
        print("\n第7步: 计算最终预测结果...")

    # 创建最终预测数据的副本
    final_predicted_mid_data_list = [arr.copy() for arr in complete_predicted_mid_data_list]

    # 使用优化后的通量值重新计算预测
    base_prediction_function(
        optimized_flux_vector,
        final_predicted_mid_data_list,
        operation_list
    )

    # 处理混合代谢物
    if mix_operation_list:
        mix_mid_data_list(
            optimized_flux_vector,
            final_predicted_mid_data_list,
            mix_operation_list,
            mix_ratio_multiplier
        )

    # 步骤8: 计算预测与实验数据的比较
    if verbose:
        print("\n第8步: 计算预测与实验数据比较...")

    comparison_data = []
    for exp_metabolite, exp_data in experimental_mid_data_obj_dict.items():
        if exp_metabolite in model_metabolite_to_standard_name_dict:
            model_metabolite = model_metabolite_to_standard_name_dict[exp_metabolite]
            if model_metabolite in bare_metabolite_dim_dict:
                emu_name = f"{model_metabolite}__{'1' * bare_metabolite_dim_dict[model_metabolite]}"

                if emu_name in emu_name_index_dict:
                    # 修复：emu_name_index_dict可能返回列表或元组，确保只取第一个索引
                    if isinstance(emu_name_index_dict[emu_name], (list, tuple)):
                        predicted_mid_index = emu_name_index_dict[emu_name][0]
                    else:
                        predicted_mid_index = emu_name_index_dict[emu_name]

                    predicted_mid = final_predicted_mid_data_list[predicted_mid_index]

                    # 确保exp_data有mid_list属性并且不为空
                    if hasattr(exp_data, 'mid_list') and exp_data.mid_list is not None:
                        exp_mid_data = np.array(exp_data.mid_list)

                        # 确保比较长度一致
                        min_length = min(len(predicted_mid), len(exp_mid_data))

                        # 计算RMSE
                        rmse = np.sqrt(np.mean((predicted_mid[:min_length] - exp_mid_data[:min_length]) ** 2))

                        comparison_data.append({
                            'experimental_metabolite': exp_metabolite,
                            'model_metabolite': model_metabolite,
                            'predicted_mid': predicted_mid,
                            'experimental_mid': exp_mid_data,
                            'RMSE': rmse
                        })

    if verbose and comparison_data:
        print("\n代谢物MID预测与实验比较:")
        for data in comparison_data:
            print(f"\n{data['experimental_metabolite']} (模型: {data['model_metabolite']}):")
            print(f"  RMSE: {data['RMSE']:.4f}")
            exp_mid = data['experimental_mid']
            pred_mid = data['predicted_mid']
            print("  预测MID vs 实验MID:")
            for i in range(min(len(exp_mid), len(pred_mid))):
                print(f"  M+{i}: {pred_mid[i]:.4f} vs {exp_mid[i]:.4f}")

    # 步骤9: 计算流量控制系数和敏感性分析
    if verbose:
        print("\n第9步: 流量控制分析...")

    # 计算流量控制系数
    flux_control_coefficients = {}
    # 对每个流量计算控制系数
    for flux_name, idx in flux_name_index_dict.items():
        # 跳过流量为零的反应
        if abs(optimized_flux_vector[idx]) < 1e-6:
            continue

        # 将当前流量增加1%
        perturbed_flux = optimized_flux_vector.copy()
        perturbed_flux[idx] *= 1.01

        # 计算扰动后的预测MID
        perturbed_mid_data_list = [arr.copy() for arr in complete_predicted_mid_data_list]
        base_prediction_function(
            perturbed_flux,
            perturbed_mid_data_list,
            operation_list
        )

        # 如果有混合操作，处理混合代谢物
        if mix_operation_list:
            mix_mid_data_list(
                perturbed_flux,
                perturbed_mid_data_list,
                mix_operation_list,
                mix_ratio_multiplier
            )

        # 计算扰动对目标代谢物MID的影响
        coefficients = {}
        for exp_metabolite in experimental_mid_data_obj_dict:
            if exp_metabolite in model_metabolite_to_standard_name_dict:
                model_metabolite = model_metabolite_to_standard_name_dict[exp_metabolite]
                if model_metabolite in bare_metabolite_dim_dict:
                    emu_name = f"{model_metabolite}__{'1' * bare_metabolite_dim_dict[model_metabolite]}"

                    if emu_name in emu_name_index_dict:
                        # 修复：emu_name_index_dict可能返回列表或元组，确保只取第一个索引
                        if isinstance(emu_name_index_dict[emu_name], (list, tuple)):
                            mid_idx = emu_name_index_dict[emu_name][0]
                        else:
                            mid_idx = emu_name_index_dict[emu_name]

                        original_mid = final_predicted_mid_data_list[mid_idx]
                        perturbed_mid = perturbed_mid_data_list[mid_idx]

                        # 计算相对变化
                        relative_change = np.sum(np.abs(perturbed_mid - original_mid)) / np.sum(original_mid)
                        coefficients[exp_metabolite] = relative_change / 0.01  # 归一化为1%变化的影响

        flux_control_coefficients[flux_name] = coefficients

    if verbose and flux_control_coefficients:
        print("\n流量控制系数 (对1%流量变化的敏感性):")
        for flux_name, coeffs in flux_control_coefficients.items():
            print(f"\n  {flux_name}:")
            for metab, coeff in coeffs.items():
                print(f"    {metab}: {coeff:.4f}")


    # 返回完整结果
    return {
        'base_results': base_results,
        'nested_mix_equation_dict': nested_mix_equation_dict,
        'mix_ratio_name_index_dict': mix_ratio_name_index_dict,
        'initial_flux_vector': initial_flux_vector,
        'optimized_flux_vector': optimized_flux_vector,
        'target_emu_name_list': target_emu_name_list,
        'target_indices': target_indices,
        'predicted_mid_data': final_predicted_mid_data_list,
        'comparison_data': comparison_data,
        'flux_control_coefficients': flux_control_coefficients,
        'optimization_result': opt_result
    }