from ...common.packages import np
from ...common.config import CoreConstants, ParamName, ModelKeyword
from ...common.functions import remove_numerical_error, natural_dist, isdigit, full_emu_name_constructor, \
    mix_flux_name_constructor, tissue_specific_name_constructor, compartmental_mid_name_constructor, np_log_eps, \
    check_if_subsequence



# 这个是融合data和model数据可能不匹配的问题，这个很重要
# mixing_equation_constructor 函数是解决质谱实验数据与代谢网络模型中代谢物不匹配问题的关键组件。
# 该函数处理以下几种常见的不匹配情况：
def mixing_equation_constructor(
        experimental_mid_data_obj_dict,
        model_target_metabolite_compartment_dict,
        model_metabolite_to_standard_name_dict,
        metabolite_bare_metabolite_name_dict,
        specific_flux_range_dict,
        common_mix_ratio_range,
        mix_ratio_multiplier,
        list_of_case_name=None):
    """
    Experimental metabolite: one metabolite that can be distinguished by mass spec. eg. 3PG/2PG
    Elemental metabolite: one metabolite that for a whole cell. eg. oxaloacetate
    Compartmental metabolite: one metabolite in each compartment. eg. OAC_m

    Compartment mixing: OAC_m + OAC_c -> OAC
    Experimental mixing: 3PG + 2PG -> 3PG/2PG

    If some element metabolite does not exist, just jump it.

    :param mix_ratio_multiplier:
    :param common_mix_ratio_range:
    :param specific_flux_range_dict:
    :param experimental_mid_data_obj_dict:
    :param mixed_compartment_list:
    :param model_target_metabolite_compartment_dict:
    :param model_metabolite_to_standard_name_dict:
    :return:
    """

    # 步骤2: 从实验数据角度构建混合方程
    def data_mixing_equation_dict_generator(_experimental_mid_data_obj_dict):
        _data_mixing_equation_dict = {}
        for experimental_metabolite_name, experimental_mid_data in _experimental_mid_data_obj_dict.items():
            # 检查实验代谢物是否是多种代谢物的组合
            if experimental_mid_data.excluded_from_mfa:
                continue
            nested_metabolite_compartment_list = []
            if experimental_mid_data.combined:
                standard_name_list = experimental_mid_data.combined_standard_name_list
            else:
                standard_name_list = [experimental_mid_data.name]

            # 在不同细胞组分中寻找匹配
            for elemental_standard_name in standard_name_list:
                # 检查该代谢物在模型的不同组分中是否存在
                current_tissue_elemental_name_list = []
                for tissue_name in experimental_mid_data.tissue:
                    current_compartmental_elemental_name_list = []
                    tissue_specific_elemental_standard_name = tissue_specific_name_constructor(
                        elemental_standard_name, tissue_name)
                    for compartment_name in experimental_mid_data.compartment:
                        tissue_specific_elemental_complete_name = compartmental_mid_name_constructor(
                            tissue_specific_elemental_standard_name, compartment_name)
                        if tissue_specific_elemental_complete_name in model_mixing_equation_dict:
                            current_compartmental_elemental_name_list.append((
                                tissue_specific_elemental_complete_name,
                                model_mixing_equation_dict[tissue_specific_elemental_complete_name]))
                            # elemental_to_experimental_metabolite_dict_from_data[elemental_complete_name] = \
                            #     experimental_metabolite_name
                    if len(current_compartmental_elemental_name_list) != 0:
                        current_tissue_elemental_name_list.append((
                            tissue_specific_elemental_standard_name, current_compartmental_elemental_name_list))
                if len(current_tissue_elemental_name_list) != 0:
                    nested_metabolite_compartment_list.append((
                        tissue_specific_name_constructor(elemental_standard_name, experimental_mid_data.tissue),
                        current_tissue_elemental_name_list))
            if len(nested_metabolite_compartment_list) != 0:
                _data_mixing_equation_dict[experimental_metabolite_name] = \
                    nested_metabolite_compartment_list
        return _data_mixing_equation_dict

    def check_if_new_balance_tuple_exist(new_balance_tuple, _mix_ratio_balance_dict):
        subseq_bool = False
        reverse = False
        target_item = None
        for exist_balance_tuple in _mix_ratio_balance_dict.keys():
            subseq_bool, reverse = check_if_subsequence(new_balance_tuple, exist_balance_tuple)
            if subseq_bool:
                target_item = exist_balance_tuple
                break
        return subseq_bool, reverse, target_item

    # 函数最核心的部分是通过嵌套的结构表示代谢物之间的复杂关系：
    # 嵌套依赖分析：通过analyze_nested_dependence函数递归处理复杂的混合关系
    def analyze_nested_dependence(_data_mixing_list, _metabolite_name):
        if _data_mixing_list is None:
            return _metabolite_name
        returned_obj_list = []
        for new_metabolite_name, new_nested_obj in _data_mixing_list:
            returned_obj_list.append(analyze_nested_dependence(new_nested_obj, new_metabolite_name))
        if len(returned_obj_list) == 1:
            return returned_obj_list[0]
        else:
            result_dict = {}
            new_balance_list = []
            for index, returned_metabolite in enumerate(returned_obj_list):
                mix_ratio_name = mix_flux_name_constructor(_metabolite_name, index)
                result_dict[mix_ratio_name] = returned_metabolite
                new_balance_list.append(mix_ratio_name)
                if mix_ratio_name not in mix_ratio_name_index_dict:
                    mix_ratio_name_index_dict[mix_ratio_name] = len(mix_ratio_name_index_dict)
            new_balance_tuple = tuple(sorted(new_balance_list))
            if new_balance_tuple not in mix_ratio_balance_dict:
                new_balance_exist, replace, target_item = check_if_new_balance_tuple_exist(
                    new_balance_tuple, mix_ratio_balance_dict)
                if not new_balance_exist:
                    mix_ratio_balance_dict[new_balance_tuple] = None
                elif new_balance_exist and replace:
                    del mix_ratio_balance_dict[target_item]
                    mix_ratio_balance_dict[new_balance_tuple] = None
            return result_dict

    def nested_mix_equation_dict_generator(_data_mixing_equation_dict):
        _nested_mix_equation_dict = {}
        for experimental_metabolite_name, data_mixing_list in _data_mixing_equation_dict.items():
            if experimental_metabolite_name not in _nested_mix_equation_dict:
                _nested_mix_equation_dict[experimental_metabolite_name] = analyze_nested_dependence(
                    data_mixing_list, experimental_metabolite_name)
        return _nested_mix_equation_dict

    # 步骤1: 从模型角度构建混合字典
    model_mixing_equation_dict = {}
    for tissue_name, each_tissue_metabolite_compartment_dict in model_target_metabolite_compartment_dict.items():
        for compartment_name, model_target_metabolite_set in each_tissue_metabolite_compartment_dict.items():
            # 构建模型中每个代谢物的完整名称(包含组织和细胞器信息)
            for model_metabolite_name in model_target_metabolite_set:
                bared_metabolite_name = metabolite_bare_metabolite_name_dict[model_metabolite_name]
                standard_metabolite_name = model_metabolite_to_standard_name_dict[bared_metabolite_name]
                tissue_specific_elemental_complete_name = compartmental_mid_name_constructor(
                    tissue_specific_name_constructor(standard_metabolite_name, tissue_name),
                    compartment_name)
                if tissue_specific_elemental_complete_name not in model_mixing_equation_dict:
                    model_mixing_equation_dict[tissue_specific_elemental_complete_name] = []
                model_mixing_equation_dict[tissue_specific_elemental_complete_name].append(
                    (model_metabolite_name, None))

    # 混合比例变量创建：为每个混合关系创建混合比例变量（如MIX:3PG_2PG_0）
    # Mapping from experimental_metabolite_name to model_metabolite_name
    # {'experimental_metabolite1': {
    #     'MIX:exp_1': 'elemental_metabolite1',
    #     'MIX:exp_2': {   # 'elemental_metabolite2'
    #         'MIX_ele1': 'compartmental_metabolite1',
    #         'MIX_ele2': 'compartmental_metabolite2'}
    #     }
    #  'experimental_metabolite2': ...
    # }
    mix_ratio_balance_dict = {}
    mix_ratio_name_index_dict = {}
    if list_of_case_name is not None:
        nested_mix_equation_dict = {}
        for case_name in list_of_case_name:
            data_mixing_equation_dict = data_mixing_equation_dict_generator(experimental_mid_data_obj_dict[case_name])
            nested_mix_equation_dict[case_name] = nested_mix_equation_dict_generator(data_mixing_equation_dict)
    else:
        data_mixing_equation_dict = data_mixing_equation_dict_generator(experimental_mid_data_obj_dict)
        nested_mix_equation_dict = nested_mix_equation_dict_generator(data_mixing_equation_dict)
    mix_ratio_balance_list = list(mix_ratio_balance_dict.keys())

    if common_mix_ratio_range is None:
        common_mix_ratio_range = (0, 1)
    mix_range = (common_mix_ratio_range[0] * mix_ratio_multiplier, common_mix_ratio_range[1] * mix_ratio_multiplier)
    updated_specific_flux_range_dict = dict(specific_flux_range_dict)
    for current_mix_ratio_name in mix_ratio_name_index_dict.keys():
        updated_specific_flux_range_dict[current_mix_ratio_name] = mix_range
    # 返回结果
    # nested_mix_equation_dict: 记录实验代谢物如何由模型代谢物组成
    # mix_ratio_name_index_dict: 混合比例变量的索引
    # mix_ratio_balance_list: 混合比例的平衡约束
    # updated_specific_flux_range_dict: 包含混合比例变量范围的更新字典
    return nested_mix_equation_dict, mix_ratio_name_index_dict, mix_ratio_balance_list, updated_specific_flux_range_dict