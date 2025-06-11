from ..common.config import ModelKeyword
from ..common.functions import tissue_specific_name_constructor
from ..common.classes import TransformDict

from .model_class import EMUMIDDimDict, MFAModel, UserDefinedModel
from .common_model_process_functions import model_preprocess, compart_all_metabolites, \
    flux_balance_equation_generator, metabolite_carbon_number_verification
from .emu_analyzer_functions import emu_equation_analyzer, emu_matrix_equation_generator, emu_dependency_analyzer
# These functions are essential for constructing and processing MFA models,
# enabling the analysis of metabolic fluxes in different tissues and conditions.


# 这是model模块里面最重要的
# This function constructs a common MFA model from a user-defined model.
# It preprocesses the user-defined model to generate various dictionaries and sets required for MFA.
# It generates EMU equations, matrix equations, and dependencies.
# It returns an MFAModel object containing all the necessary data for MFA.
# 该调用传递了关键参数 symmetrical_metabolite_set，从代码中可见包含 {'SUC_m', 'FUM_m'} (琥珀酸和延胡索酸)。
def common_model_constructor(user_defined_model: UserDefinedModel, new_format=False):
    # 该函数构建了一个数学框架，用于后续:
    # 接收通量值作为输入
    # 计算预测的MID分布
    # 准备与实验数据比较的基础
    # 1. 预处理模型结构
    (
        metabolite_reaction_dict, product_reaction_dict_for_emu, composite_reaction_dict,
        flux_name_index_dict, input_metabolite_name_set, complete_metabolite_dim_dict,
        target_metabolite_name_list) = model_preprocess(
        user_defined_model.reaction_list, user_defined_model.symmetrical_metabolite_set,
        user_defined_model.added_input_metabolite_set, user_defined_model.emu_excluded_metabolite_set,
        user_defined_model.balance_excluded_metabolite_set, user_defined_model.target_metabolite_list,
        user_defined_model.composite_reaction_list)

    all_target_metabolite_name_set = list(complete_metabolite_dim_dict.keys() - input_metabolite_name_set)
    # all_target_metabolite_name_set = {
    #     key: value for key, value in complete_metabolite_dim_dict.items() if key not in input_metabolite_name_set}
    complete_tissue_metabolite_compartment_dict, metabolite_bare_metabolite_name_dict = compart_all_metabolites(
        complete_metabolite_dim_dict.keys(), user_defined_model.model_compartment_set)
    bare_metabolite_dim_dict = metabolite_carbon_number_verification(
        complete_metabolite_dim_dict, metabolite_bare_metabolite_name_dict)
    flux_balance_matrix, flux_balance_right_side_vector = flux_balance_equation_generator(
        metabolite_reaction_dict, composite_reaction_dict, flux_name_index_dict)

    # 2. 构建EMU方程系统
    emu_mid_equation_dict, input_emu_dict = emu_equation_analyzer(
        product_reaction_dict_for_emu, target_metabolite_name_list, input_metabolite_name_set,
        complete_metabolite_dim_dict)
    output_emu_mid_equation_dict = emu_matrix_equation_generator(emu_mid_equation_dict, input_emu_dict)
    emu_name_dependency_dict, complete_emu_obj_index_dict, new_input_emu_dict = emu_dependency_analyzer(
        product_reaction_dict_for_emu, input_metabolite_name_set, complete_metabolite_dim_dict)
    assert input_emu_dict.keys() == new_input_emu_dict.keys()
    complete_emu_dim_dict = EMUMIDDimDict()
    # 3. 创建最终的MFA模型对象
    mfa_model = MFAModel(
        input_emu_dict, target_metabolite_name_list, output_emu_mid_equation_dict, emu_name_dependency_dict,
        composite_reaction_dict, complete_emu_dim_dict, complete_emu_obj_index_dict, flux_name_index_dict,
        flux_balance_matrix, flux_balance_right_side_vector, bare_metabolite_dim_dict,
        metabolite_bare_metabolite_name_dict, complete_tissue_metabolite_compartment_dict,
        input_metabolite_name_set, all_target_metabolite_name_set,
        user_defined_model.model_metabolite_to_standard_name_dict, user_defined_model)
    # 3. 对称代谢物的意义
    # 对称代谢物（如琥珀酸）的碳原子排列有多种等价方式（例如ABCD或DCBA）
    # 在计算13C标记分布时，必须考虑这种对称性
    # 例如在TCA循环中，琥珀酸→延胡索酸反应中，由于琥珀酸的对称性，延胡索酸可能有不同的碳标记模式
    # 4. 后续处理
    # 处理完对称代谢物后，common_model_constructor() 继续使用处理后的数据进行：
    # EMU方程分析 (emu_equation_analyzer)
    # EMU矩阵方程生成 (emu_matrix_equation_generator)
    # EMU依赖性分析 (emu_dependency_analyzer)
    # 最终构建并返回完整的MFA模型
    # 这些步骤确保了在代谢通量分析中正确处理代谢物对称性问题。
    return mfa_model