import numpy as np
import scipy.linalg
from scripts.src.core.common.config import CoreConstants
from scripts.src.core.model.model_class import Reaction, Node
from scripts.src.core.model.common_model_process_functions import model_preprocess
from scripts.src.core.model.emu_analyzer_functions import emu_equation_analyzer, emu_matrix_equation_generator
from scripts.src.core.solver.slsqp_numba_solver.common_generator_functions import (
    emu_graph_constructor_optimized_pure_python
)
from scripts.src.core.solver.slsqp_numba_solver.python_initializer_optimizer_functions import (
    base_prediction_function, solver_objective_func, mix_mid_data_list
)
from scripts.src.core.solver.solver_construction_functions.common_construct_functions import (
    mixing_equation_constructor
)
from common_functions import generate_input_emu_data, prepare_mid_data_list, print_mid_results
from scripts.src.core.data.data_class import MIDData
from scipy.optimize import minimize, linprog
from scripts.src.core.model.common_model_process_functions import flux_balance_equation_generator
from scripts.data.common_functions import common_data_loader
from scripts.src.common.config import DataType
from scripts.data.hct116_cultured_cell_line.data_metabolite_to_standard_name_dict import data_metabolite_to_standard_name_dict
from typing import List, Tuple, Dict, Any, Optional, Union, Callable


def setup_tca_network_from_paper():
    """
    根据论文DOI:10.1016/j.ymben.2006.09.001设置TCA循环代谢网络
    基于论文中的实验数据和代谢路径构建TCA网络
    """
    # 根据论文，构建TCA循环的关键反应
    # 反应1: OAC + AcCoA → Cit (碳原子转换: abcd + ef → dcbfea)
    v1 = Reaction(id='v1',
                 sub=[Node(name='OAC', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd']),
                      Node(name='AcCoA', coefficient=1, carbon_composition_string=['e', 'f'])],
                 pro=[Node(name='Cit', coefficient=1, carbon_composition_string=['d', 'c', 'b', 'f', 'e', 'a'])],
                 reverse=False)

    # 反应2: Cit → AKG + CO2 (碳原子转换: abcdef → abcde + f)
    v2 = Reaction(id='v2',
                 sub=[Node(name='Cit', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd', 'e', 'f'])],
                 pro=[Node(name='AKG', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd', 'e']),
                      Node(name='CO2', coefficient=1, carbon_composition_string=['f'])],
                 reverse=False)

    # 反应3: AKG → Glu
    v3 = Reaction(id='v3',
                 sub=[Node(name='AKG', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd', 'e'])],
                 pro=[Node(name='Glu', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd', 'e'])],
                 reverse=False)

    # 反应4: AKG → Suc + CO2
    v4 = Reaction(id='v4',
                 sub=[Node(name='AKG', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd', 'e'])],
                 pro=[Node(name='Suc', coefficient=1, carbon_composition_string=['b', 'c', 'd', 'e']),
                      Node(name='CO2', coefficient=1, carbon_composition_string=['a'])],
                 reverse=False)

    # 反应5: Suc → Fum (考虑碳原子对称性)
    v5 = Reaction(id='v5',
                 sub=[Node(name='Suc', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd'])],
                 pro=[Node(name='Fum', coefficient=0.5, carbon_composition_string=['a', 'b', 'c', 'd']),
                      Node(name='Fum', coefficient=0.5, carbon_composition_string=['d', 'c', 'b', 'a'])],
                 reverse=False)

    # 反应6: Fum → OAC
    v6 = Reaction(id='v6',
                 sub=[Node(name='Fum', coefficient=0.5, carbon_composition_string=['a', 'b', 'c', 'd']),
                      Node(name='Fum', coefficient=0.5, carbon_composition_string=['d', 'c', 'b', 'a'])],
                 pro=[Node(name='OAC', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd'])],
                 reverse=False)

    # 反应7: OAC → Fum
    v7 = Reaction(id='v7',
                 sub=[Node(name='OAC', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd'])],
                 pro=[Node(name='Fum', coefficient=0.5, carbon_composition_string=['a', 'b', 'c', 'd']),
                      Node(name='Fum', coefficient=0.5, carbon_composition_string=['d', 'c', 'b', 'a'])],
                 reverse=False)

    # 反应8: Asp → OAC
    v8 = Reaction(id='v8',
                 sub=[Node(name='Asp', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd'])],
                 pro=[Node(name='OAC', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd'])],
                 reverse=False)

    # 构建代谢物反应字典
    metabolite_reaction_dict = {
        # 输入代谢物
        'Asp': [],
        'AcCoA': [],
        # TCA循环代谢物
        'Cit': [v1],
        'AKG': [v2],
        'Suc': [v4],
        'Fum': [v5, v7],
        'OAC': [v6, v8],
        # 输出代谢物
        'Glu': [v3],
        'CO2': [v2, v4]
    }

    # 流量平衡反应字典
    flux_balance_reaction_dict = {
        'Asp': [v8],
        'AcCoA': [v1],
        'Cit': [v1, v2],
        'AKG': [v2, v3, v4],
        'Suc': [v4, v5],
        'Fum': [v5, v7, v6],
        'OAC': [v6, v8, v1, v7],
        'Glu': [v3],
        'CO2': [v2, v4]
    }

    # 代谢物碳原子数
    complete_metabolite_dim_dict = {
        'Asp': 4,
        'AcCoA': 2,
        'Cit': 6,
        'AKG': 5,
        'Suc': 4,
        'Fum': 4,
        'OAC': 4,
        'Glu': 5,
        'CO2': 1
    }

    # 流量名称索引字典
    flux_name_index_dict = {
        'v1': 0,
        'v2': 1,
        'v3': 2,
        'v4': 3,
        'v5': 4,
        'v6': 5,
        'v7': 6,
        'v8': 7
    }

    return metabolite_reaction_dict, flux_balance_reaction_dict, complete_metabolite_dim_dict, flux_name_index_dict


def create_simulated_experimental_data():
    """
    根据论文数据创建模拟的质谱实验数据
    基于TCA循环的13C标记实验数据
    """
    experimental_mid_data_obj_dict = {}
    
    # 模拟Glu的质谱实验数据 (基于论文中的TCA循环13C标记实验)
    # 创建模拟的MID数据 - 基于[1,2-13C2]乙酰CoA标记实验
    glu_mid_data = MIDData(
        raw_name='glutamate',
        raw_data_list=[0.3464, 0.2695, 0.2708, 0.0807, 0.0286, 0.0039],  # M+0到M+5的分布
        compartment_list=['c'],  # 细胞质区室
        tissue_list=['default']  # 默认组织
    )
    glu_mid_data.mid_list = [0.3464, 0.2695, 0.2708, 0.0807, 0.0286, 0.0039]  # M+0到M+5的分布
    glu_mid_data.data_vector = np.array(glu_mid_data.mid_list)
    glu_mid_data.name = 'glutamate'
    
    experimental_mid_data_obj_dict['glutamate'] = glu_mid_data
    
    # 如果需要更多代谢物的实验数据，可以在这里添加
    # 例如琥珀酸、柠檬酸等的MID数据
    
    return experimental_mid_data_obj_dict


def run_optimization_with_constraints(
        objective_function: Callable,
        initial_flux_vector: np.ndarray,
        flux_bounds: List[Tuple[float, float]],
        flux_balance_matrix: Optional[np.ndarray] = None,
        flux_balance_rhs: Optional[np.ndarray] = None,
        verbose: bool = True
) -> Any:
    """
    执行带约束的优化
    """
    constraints = []
    
    # 添加质量平衡约束
    if flux_balance_matrix is not None and flux_balance_rhs is not None:
        if flux_balance_matrix.shape[0] > 0:
            constraints.append({
                'type': 'eq',
                'fun': lambda x: flux_balance_matrix @ x - flux_balance_rhs,
                'jac': lambda x: flux_balance_matrix
            })
    
    try:
        result = minimize(
            objective_function,
            initial_flux_vector,
            method='SLSQP',
            bounds=flux_bounds,
            constraints=constraints,
            options={
                'disp': verbose,
                'maxiter': 2000,
                'ftol': 1e-10
            }
        )
        return result
    except Exception as e:
        if verbose:
            print(f"优化失败: {e}")
        # 返回失败的结果对象
        return type('Result', (), {
            'success': False,
            'message': str(e),
            'x': initial_flux_vector,
            'fun': float('inf')
        })()

#TODO: model_preprocess 处理模型

#TODO：加载实验数据（模拟实验数据即可），包括target_metabolites MID与input metabolites ratio list

#TODO: EMU matrix与EMU graph构建(operation_list)

#TODO: loss_operation_list构建

#TODO: base_prediction_function构建

#TODO: flux_balance和bounds构建（包括v8的特殊值）

#TODO: initial_vector生成（利用flux_balance和bounds）

#TODO: 执行optimization（initial_vector, base_prediction_function, flux_balance, bounds），包括初始流量、初始obj、优化后的流量、优化后obj

def common_test_pipeline():
    """
    基于论文DOI:10.1016/j.ymben.2006.09.001的TCA实验数据执行完整的MFA分析流程
    """
    print("启动基于论文TCA实验数据的代谢通量分析(MFA)流程...\n")
    
    # TODO 1: model_preprocess 处理模型
    print("步骤1: 模型预处理...")
    
    # 设置TCA网络
    metabolite_reaction_dict, flux_balance_reaction_dict, complete_metabolite_dim_dict, flux_name_index_dict = setup_tca_network_from_paper()
    
    # 设置输入代谢物
    input_metabolite_name_set = {'AcCoA', 'Asp'}
    
    # 预处理模型
    symmetrical_metabolite_set = {'Suc', 'Fum'}
    balance_excluded_metabolite_set = {'CO2'}
    
    processed_model = model_preprocess(
        metabolite_reaction_dict,
        input_metabolite_name_set,
        complete_metabolite_dim_dict,
        flux_name_index_dict,
        symmetrical_metabolite_set=symmetrical_metabolite_set,
        emu_excluded_metabolite_set=balance_excluded_metabolite_set,
        balance_excluded_metabolite_set=balance_excluded_metabolite_set,
    )
    
    # 更新预处理后的数据
    metabolite_reaction_dict = processed_model['metabolite_reaction_dict']
    input_metabolite_name_set = processed_model['input_metabolite_name_set']
    complete_metabolite_dim_dict = processed_model['complete_metabolite_dim_dict']
    flux_name_index_dict = processed_model['flux_name_index_dict']
    
    print(f"  - 代谢物数量: {len(metabolite_reaction_dict)}")
    print(f"  - 输入代谢物: {input_metabolite_name_set}")
    print(f"  - 反应流量数: {len(flux_name_index_dict)}")
    
    # TODO 2: 加载实验数据（模拟实验数据），包括target_metabolites MID与input metabolites ratio list
    print("\n步骤2: 加载实验数据...")
    
    # 创建模拟质谱实验数据
    experimental_mid_data_obj_dict = create_simulated_experimental_data()
    
    # 设置目标代谢物
    target_metabolite_name_list = ['Glu']
    
    # 设置输入代谢物标记模式 - 基于论文的[1,2-13C2]乙酰CoA标记实验
    input_metabolite_dict = {
        "AcCoA": [
            {
                "ratio_list": [0, 1],  # 第2个碳原子标记
                "abundance": 0.25,  # 25%的分子有这种标记
            },
            {
                "ratio_list": [1, 1],  # 两个碳原子都标记
                "abundance": 0.25,  # 25%的分子有这种标记
            },
            {
                "ratio_list": [0, 0],  # 没有标记
                "abundance": 0.5,  # 50%的分子没有标记
            }
        ]
    }
    
    print(f"  - 实验数据包含 {len(experimental_mid_data_obj_dict)} 个代谢物")
    print(f"  - 目标代谢物: {target_metabolite_name_list}")
    
    # TODO 3: EMU matrix与EMU graph构建(operation_list)
    print("\n步骤3: EMU方程分析和矩阵方程构建...")
    
    # EMU方程分析
    emu_equations, input_emu_dict = emu_equation_analyzer(
        metabolite_reaction_dict,
        target_metabolite_name_list,
        input_metabolite_name_set,
        complete_metabolite_dim_dict
    )
    
    print(f"  - EMU方程数量: {len(emu_equations)}")
    
    # 生成输入EMU数据
    input_emu_data_dict = generate_input_emu_data(
        input_metabolite_dict,
        input_emu_dict,
        natural_abundance=0.01
    )
    
    # 构建EMU矩阵方程
    emu_matrix_equations = emu_matrix_equation_generator(emu_equations, input_emu_data_dict)
    
    print(f"  - 矩阵方程数量: {len(emu_matrix_equations)}")
    
    # 构建EMU图
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
    
    print(f"  - EMU图构建完成")
    
    # TODO 4: loss_operation_list构建
    print("\n步骤4: 构建损失操作列表...")
    
    # 获取EMU图结果
    input_emu_data_list = graph_results[0]
    operation_list = graph_results[1]
    emu_name_index_dict = graph_results[3]
    
    # 构建损失操作列表
    loss_operation_list = []
    user_metabolite_to_standard_name_dict = {'glutamate': 'Glu'}
    
    for exp_metabolite, exp_data in experimental_mid_data_obj_dict.items():
        if exp_metabolite in user_metabolite_to_standard_name_dict:
            model_metabolite = user_metabolite_to_standard_name_dict[exp_metabolite]
            carbon_num = complete_metabolite_dim_dict.get(model_metabolite)
            if carbon_num is not None:
                emu_name = f"{model_metabolite}__{'1' * carbon_num}"
                if emu_name in emu_name_index_dict:
                    predicted_mid_index = emu_name_index_dict[emu_name]
                    if isinstance(predicted_mid_index, (list, tuple)):
                        predicted_mid_index = predicted_mid_index[0]
                    
                    # 转换实验数据为numpy数组
                    exp_mid_array = np.array(exp_data.mid_list)
                    valid_index_array = np.arange(len(exp_mid_array), dtype=np.int32)
                    
                    loss_operation_list.append([
                        exp_metabolite, predicted_mid_index, exp_mid_array, valid_index_array
                    ])
    
    print(f"  - 损失操作数量: {len(loss_operation_list)}")
    
    # TODO 5: base_prediction_function构建
    print("\n步骤5: 构建基础预测函数...")
    
    # 准备MID数据列表
    complete_predicted_mid_data_list = prepare_mid_data_list(graph_results)
    
    print(f"  - 预测数据列表长度: {len(complete_predicted_mid_data_list)}")
    
    # TODO 6: flux_balance和bounds构建（包括v8的特殊值）
    print("\n步骤6: 构建通量平衡约束和边界条件...")
    
    # 构建通量平衡矩阵
    metabolite_reaction_dict_for_bal = {
        'OAC': ({'v1': 1.0, 'v7': 1.0}, {'v8': 1.0, 'v6': 1.0}),
        'Cit': ({'v2': 1.0}, {'v1': 1.0}),
        'AKG': ({'v3': 1.0, 'v4': 1.0}, {'v2': 1.0}),
        'Suc': ({'v5': 1.0}, {'v4': 1.0}),
        'Fum': ({'v6': 1.0}, {'v5': 1.0, 'v7': 1.0})
    }
    
    flux_balance_matrix, flux_balance_right_side_vector = flux_balance_equation_generator(
        metabolite_reaction_dict=metabolite_reaction_dict_for_bal,
        composite_reaction_dict={},
        flux_name_index_dict=flux_name_index_dict
    )
    
    # 设置通量边界（包括v8的特殊值，基于论文数据v8=50）
    flux_bounds = []
    for flux_name, idx in flux_name_index_dict.items():
        if flux_name == 'v8':
            # v8设置为固定值50（基于论文实验数据）
            flux_bounds.append((49.9, 50.1))
        else:
            flux_bounds.append((0.01, 1000.0))
    
    print(f"  - 流量平衡矩阵形状: {flux_balance_matrix.shape if flux_balance_matrix is not None else 'None'}")
    print(f"  - v8流量设置为特殊值: 50±0.1")
    
    # TODO 7: initial_vector生成（利用flux_balance和bounds）
    print("\n步骤7: 生成初始流量向量...")
    
    # 使用线性规划生成初始向量
    num_vars = len(flux_name_index_dict)
    c = np.random.random(num_vars) - 0.5  # 随机目标函数
    
    # 线性规划求解
    result = linprog(
        c=c,
        A_eq=flux_balance_matrix if flux_balance_matrix is not None and len(flux_balance_matrix) > 0 else None,
        b_eq=flux_balance_right_side_vector if flux_balance_right_side_vector is not None and len(flux_balance_right_side_vector) > 0 else None,
        bounds=flux_bounds,
        method='highs'
    )
    
    if result.success:
        initial_flux_vector = result.x
        print("  - 线性规划成功生成初始流量向量")
    else:
        # 如果线性规划失败，使用论文中的实验值作为初始值
        initial_flux_vector = np.array([100, 100, 50, 50, 50, 125, 75, 50], dtype=float)
        print("  - 使用论文实验值作为初始流量向量")
    
    for flux_name, idx in flux_name_index_dict.items():
        print(f"    {flux_name}: {initial_flux_vector[idx]:.2f}")
    
    # TODO 8: 执行optimization（initial_vector, base_prediction_function, flux_balance, bounds）
    print("\n步骤8: 执行优化操作...")
    
    # 定义目标函数
    def objective_function(x):
        curr_predicted_mid_data_list = [arr.copy() for arr in complete_predicted_mid_data_list]
        
        # 计算预测
        base_prediction_function(x, curr_predicted_mid_data_list, operation_list)
        
        # 计算损失
        loss = solver_objective_func(
            flux_vector=x,
            complete_predicted_mid_data_list=curr_predicted_mid_data_list,
            operation_list=operation_list,
            mix_operation_list=[],  # 无混合操作
            mix_ratio_multiplier=1.0,
            loss_operation_list=loss_operation_list,
            loss_code=1,  # 使用均方差损失
            optimal_cross_entropy=0.0
        )
        return loss
    
    # 计算初始目标函数值
    initial_obj = objective_function(initial_flux_vector)
    print(f"  - 初始目标函数值: {initial_obj:.6f}")
    
    # 执行优化
    opt_result = run_optimization_with_constraints(
        objective_function=objective_function,
        initial_flux_vector=initial_flux_vector,
        flux_bounds=flux_bounds,
        flux_balance_matrix=flux_balance_matrix,
        flux_balance_rhs=flux_balance_right_side_vector,
        verbose=True
    )
    
    if opt_result.success:
        optimized_flux_vector = opt_result.x
        optimized_obj = opt_result.fun
        print(f"  - 优化成功! 最终目标函数值: {optimized_obj:.6f}")
    else:
        optimized_flux_vector = initial_flux_vector
        optimized_obj = initial_obj
        print(f"  - 优化失败，使用初始值")
    
    # TODO 9: 结果展示与输出
    print("\n步骤9: 结果展示与输出...")
    
    print("\n=== 流量分析结果 ===")
    print("反应流量对比 (初始值 → 优化后):")
    for flux_name, idx in flux_name_index_dict.items():
        initial_val = initial_flux_vector[idx]
        optimized_val = optimized_flux_vector[idx]
        change = ((optimized_val - initial_val) / initial_val) * 100 if initial_val != 0 else 0
        print(f"  {flux_name}: {initial_val:.2f} → {optimized_val:.2f} ({change:+.1f}%)")
    
    # 计算最终预测结果
    final_predicted_mid_data_list = [arr.copy() for arr in complete_predicted_mid_data_list]
    base_prediction_function(optimized_flux_vector, final_predicted_mid_data_list, operation_list)
    
    # 与实验数据比较
    print("\n=== 预测vs实验数据比较 ===")
    for exp_metabolite, exp_data in experimental_mid_data_obj_dict.items():
        if exp_metabolite in user_metabolite_to_standard_name_dict:
            model_metabolite = user_metabolite_to_standard_name_dict[exp_metabolite]
            carbon_num = complete_metabolite_dim_dict.get(model_metabolite)
            if carbon_num is not None:
                emu_name = f"{model_metabolite}__{'1' * carbon_num}"
                if emu_name in emu_name_index_dict:
                    predicted_mid_index = emu_name_index_dict[emu_name]
                    if isinstance(predicted_mid_index, (list, tuple)):
                        predicted_mid_index = predicted_mid_index[0]
                    
                    predicted_mid = final_predicted_mid_data_list[predicted_mid_index]
                    experimental_mid = exp_data.mid_list
                    
                    # 计算RMSE
                    min_length = min(len(predicted_mid), len(experimental_mid))
                    rmse = np.sqrt(np.mean((predicted_mid[:min_length] - np.array(experimental_mid[:min_length])) ** 2))
                    
                    print(f"\n{exp_metabolite} (模型: {model_metabolite}):")
                    print(f"  RMSE: {rmse:.4f}")
                    print("  质量分布 (M+0 to M+n):")
                    print("    预测    实验")
                    for i in range(min_length):
                        print(f"  M+{i}: {predicted_mid[i]:.4f}  {experimental_mid[i]:.4f}")
    
    print(f"\n=== 分析总结 ===")
    print(f"初始目标函数值: {initial_obj:.6f}")
    print(f"优化后目标函数值: {optimized_obj:.6f}")
    print(f"改进程度: {((initial_obj - optimized_obj) / initial_obj * 100):.2f}%")
    print(f"优化状态: {'成功' if opt_result.success else '失败'}")
    
    # 返回结果字典
    results = {
        'initial_flux_vector': initial_flux_vector,
        'optimized_flux_vector': optimized_flux_vector,
        'initial_objective': initial_obj,
        'optimized_objective': optimized_obj,
        'optimization_success': opt_result.success,
        'flux_name_index_dict': flux_name_index_dict,
        'predicted_mid_data': final_predicted_mid_data_list,
        'experimental_mid_data': experimental_mid_data_obj_dict,
        'target_emu_name_list': target_emu_name_list,
        'rmse_data': {}
    }
    
    print("\n代谢通量分析流程完成!")
    return results

    #TODO: 结果展示与输出

if __name__ == '__main__':
    results = common_test_pipeline()
