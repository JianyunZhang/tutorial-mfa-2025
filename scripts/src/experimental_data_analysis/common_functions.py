from scripts.src.core.solver.solver_construction_functions.solver_constructor import common_solver_constructor, \
    base_solver_constructor
from scripts.src.core.common.classes import MFAConfig
from ..common.functions import update_parameter_object, special_result_label_converter
from ..common.result_output_functions import solver_output

from scripts.src.common.config import Keywords
from scripts.src.common_parallel_solver.common_parallel_solver import common_solver

from .inventory import RunningMode, DataModelType
from .specific_data_model_combination.common_data_model_loader import common_data_model_function_loader
from .result_processing_functions import CurrentFinalResult, result_label_generator, \
    experimental_mid_and_raw_data_plotting, normal_result_process, hct116_result_process, \
    traditional_method_result_process
from . import config
Direct = config.Direct

# 该文件的主要功能：
# 提供数据加载和处理功能
# 构建和管理求解器对象
# 支持多种运行模式：
# 通量分析
# 原始数据绘图
# 求解器输出
# 结果处理

def mfa_data_loader(data_model_result_label_generator, data_model_parameter_key_list, total_param_list, data_wrap_obj):
    """
        加载代谢流分析(MFA)数据
        Args:
            data_model_result_label_generator: 结果标签生成器函数
            data_model_parameter_key_list: 数据模型参数键列表
            total_param_list: 总参数列表
            data_wrap_obj: 数据封装对象
        Returns:
            mfa_data_dict: MFA数据字典
            result_information_dict: 结果信息字典
    """
    mfa_data_dict = {}
    result_information_dict = {}
    for param_dict in total_param_list:
        # 从数据封装对象获取数据集
        mfa_data = data_wrap_obj.return_dataset(param_dict)
        # 生成结果标签
        result_label = data_model_result_label_generator(*[param_dict[key] for key in data_model_parameter_key_list])
        mfa_data_dict[result_label] = mfa_data
        result_information_dict[result_label] = param_dict
    return mfa_data_dict, result_information_dict


def unoptimized_initial_solution_loader(result_obj, model_label, data_label, config_label, parameter_dict):
    unoptimized_result_label = special_result_label_converter(data_label, Keywords.unoptimized)
    unoptimized_solution_array = result_obj.final_solution_data_dict[unoptimized_result_label]
    return unoptimized_solution_array


def solver_dict_constructor(parameter_label_content_dict):
    """
        构建求解器字典
        Args:
            parameter_label_content_dict: 参数标签内容字典
        Returns:
            target_solver_dict: 目标求解器字典
            same_model_dict: 相同模型标记字典
            same_data_dict: 相同数据标记字典
    """
    target_solver_dict = {}
    same_model_dict = {}
    same_data_dict = {}
    previous_label = None
    for result_label, (_, (mfa_model, mfa_data, mfa_config), _, _) in parameter_label_content_dict.items():
        # 构建求解器对象
        current_solver_obj = common_solver_constructor(
            mfa_model, mfa_data, mfa_config, name=result_label, verbose=False)
        target_solver_dict[result_label] = current_solver_obj
        # 标记相同模型和数据
        if previous_label is None:
            same_model_dict[result_label] = False
            previous_label = result_label
        else:
            same_model_dict[result_label] = True
        if result_label.startswith(previous_label) and result_label.endswith(Keywords.unoptimized):
            same_data_dict[result_label] = True
        else:
            same_data_dict[result_label] = False
    return target_solver_dict, same_model_dict, same_data_dict


def raw_experimental_data_plotting(solver_dict, result_information_dict, final_result_obj):
    target_metabolite_data_obj_dict = {}
    for result_label, solver_obj in solver_dict.items():
        if result_label == Keywords.unoptimized:
            continue
        current_data_obj_dict = {}
        for emu_name in solver_obj.target_experimental_mid_data_dict.keys():
            experimental_name = solver_obj.emu_name_experimental_name_dict[emu_name]
            current_data_obj_dict[experimental_name] = solver_obj.experimental_mid_data_obj_dict[experimental_name]
        target_metabolite_data_obj_dict[result_label] = current_data_obj_dict
    experimental_mid_and_raw_data_plotting(target_metabolite_data_obj_dict, result_information_dict, final_result_obj)


def result_display(solver_dict, final_result_obj, data_model_name):
    if data_model_name in {DataModelType.hct116_cultured_cell_line, DataModelType.hct116_cultured_cell_line_with_glns_m}:
        result_process_func = hct116_result_process
    elif data_model_name.name.endswith(Keywords.traditional_method):
        result_process_func = traditional_method_result_process
    else:
        result_process_func = normal_result_process
    final_result_obj.final_process(result_process_func, solver_dict)


def experimental_data_analysis_common_dispatcher(
        data_model_name, data_model_object, final_result_obj, test_mode, parallel_num, running_mode):
    """
        实验数据分析的通用分发器
        Args:
            data_model_name: 数据模型名称
            data_model_object: 数据模型对象
            final_result_obj: 最终结果对象
            test_mode: 测试模式标志
            parallel_num: 并行进程数
            running_mode: 运行模式
    """
    # 初始化模型和配置
    mfa_model = data_model_object.mfa_model_obj
    data_wrap_obj = data_model_object.data_wrap_obj
    total_param_list = data_model_object.total_param_list

    # 更新参数配置
    common_parameter = config.CommonParameters()
    specific_parameter = data_model_object.SpecificParameter()
    common_parameter = update_parameter_object(common_parameter, specific_parameter)

    # 创建MFA配置
    mfa_config = MFAConfig(
        common_parameter.common_flux_range, common_parameter.specific_flux_range_dict,
        common_parameter.test_dynamic_constant_flux_list, common_parameter.test_preset_constant_flux_value_dict,
        common_parameter.common_mix_ratio_range, common_parameter.mix_ratio_multiplier, config.solver_type,
        common_parameter.solver_config_dict)

    each_case_target_optimization_num, parallel_parameter_dict = config.running_settings(test_mode)
    if parallel_parameter_dict is not None and parallel_num is not None:
        parallel_parameter_dict[Keywords.processes_num] = parallel_num

    mfa_data_dict, result_information_dict = mfa_data_loader(
        data_model_object.project_name_generator, data_model_object.keyword_list,
        total_param_list, data_wrap_obj)
    parameter_label_content_dict = {}
    for result_label, mfa_data_obj in mfa_data_dict.items():
        data_label = result_label
        result_information = result_information_dict[result_label]
        other_information_dict = {}
        new_mfa_config = mfa_config.copy()
        # if result_label.endswith(Keywords.unoptimized):
        #     new_mfa_config.update_miscellaneous_config({
        #         Keywords.unoptimized: True
        #     })
        if Keywords.squared_loss in result_label:
            new_mfa_config = MFAConfig(
                common_parameter.common_flux_range, common_parameter.specific_flux_range_dict,
                common_parameter.test_dynamic_constant_flux_list, common_parameter.test_preset_constant_flux_value_dict,
                common_parameter.common_mix_ratio_range, common_parameter.mix_ratio_multiplier, config.solver_type,
                common_parameter.squared_loss_config_dict)
        if Keywords.miscellaneous in result_information:
            new_mfa_config.update_miscellaneous_config(result_information[Keywords.miscellaneous])
        parameter_label_content_dict[result_label] = (
            (None, data_label, None), (mfa_model, mfa_data_obj, new_mfa_config),
            result_information, other_information_dict
        )

    # 根据运行模式执行不同的分析
    if running_mode == RunningMode.flux_analysis:
        common_solver(
            parameter_label_content_dict, test_mode, final_result_obj, each_case_target_optimization_num,
            config.report_interval, parallel_parameter_dict, load_results=config.load_previous_results,
            predefined_initial_solution_matrix_loader=unoptimized_initial_solution_loader
        )
    else:
        solver_dict, same_model_dict, same_data_dict = solver_dict_constructor(parameter_label_content_dict)
        if running_mode == RunningMode.raw_experimental_data_plotting:
            raw_experimental_data_plotting(solver_dict, result_information_dict, final_result_obj)
        elif running_mode == RunningMode.solver_output:
            solver_output(
                solver_dict, result_information_dict, final_result_obj, same_model_dict, same_data_dict)
        elif running_mode == RunningMode.result_process:
            result_display(solver_dict, final_result_obj, data_model_name)
        else:
            raise ValueError()


def data_analysis_main(running_mode, data_model_name, test_mode, parallel_num):
    """
        数据分析主函数
        Args:
            running_mode: 运行模式
            data_model_name: 数据模型名称
            test_mode: 测试模式标志
            parallel_num: 并行进程数
    """

    if data_model_name is None:
        data_model_name = config.data_model_name
    # 加载数据模型和创建结果对象
    data_model_object = common_data_model_function_loader(data_model_name)
    final_result_obj = CurrentFinalResult(
        Direct.output_direct, Direct.common_data_direct, data_model_name, data_model_object)
    # 调用分发器执行分析
    experimental_data_analysis_common_dispatcher(
        data_model_name, data_model_object, final_result_obj, test_mode, parallel_num, running_mode)

