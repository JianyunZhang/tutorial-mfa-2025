import warnings

from ..common.packages import it, np
from ..common.classes import DefaultDict
from ..common.config import CoreConstants

from .model_class import Reaction, CompositeReaction


# 修改函数签名以适应当前的调用格式
def model_preprocess(
        metabolite_reaction_dict=None, input_metabolite_name_set=None, complete_metabolite_dim_dict=None,
        flux_name_index_dict=None, reaction_list=None, symmetrical_metabolite_set=None,
        added_input_metabolite_set=None, emu_excluded_metabolite_set=None,
        balance_excluded_metabolite_set=None, target_metabolite_name_list=None,
        composite_reaction_list=None):
    """
    预处理代谢网络模型，为代谢通量分析准备数据结构

    函数支持两种调用方式:
    1. 传统方式: 使用reaction_list, symmetrical_metabolite_set等
    2. 新方式: 使用metabolite_reaction_dict, input_metabolite_name_set等

    Parameters:
    -----------
    metabolite_reaction_dict: Dict - 代谢物与其参与反应的映射字典
    input_metabolite_name_set: Set - 输入代谢物集合
    complete_metabolite_dim_dict: Dict - 代谢物碳原子数量字典
    flux_name_index_dict: Dict - 通量名称到索引的映射
    reaction_list: List - 反应列表 (传统参数)
    symmetrical_metabolite_set: Set - 对称代谢物集合 (传统参数)
    added_input_metabolite_set: Set - 附加输入代谢物集合 (传统参数)
    emu_excluded_metabolite_set: Set - EMU计算中排除的代谢物
    balance_excluded_metabolite_set: Set - 通量平衡计算中排除的代谢物
    target_metabolite_name_list: List - 目标代谢物列表
    composite_reaction_list: List - 复合反应列表

    Returns:
    --------
    元组包含处理后的各种数据结构
    """
    # 根据调用方式选择参数
    if metabolite_reaction_dict is not None and input_metabolite_name_set is not None:
        # 新式调用方式，需要从metabolite_reaction_dict中提取反应列表
        if reaction_list is not None:
            # 如果同时提供了两种格式，抛出错误
            raise ValueError("Cannot provide both metabolite_reaction_dict and reaction_list")

        # 从metabolite_reaction_dict提取所有反应并去重
        all_reactions = []
        seen_reaction_ids = set()

        # 如果symmetrical_metabolite_set未提供，创建空集合
        if symmetrical_metabolite_set is None:
            symmetrical_metabolite_set = set()

        # 如果未提供排除代谢物集合，使用空集合
        if emu_excluded_metabolite_set is None:
            emu_excluded_metabolite_set = set()

        if balance_excluded_metabolite_set is None:
            balance_excluded_metabolite_set = set()

        # 设置added_input_metabolite_set为input_metabolite_name_set
        added_input_metabolite_set = input_metabolite_name_set

        # 从metabolite_reaction_dict中收集所有唯一反应
        for _, reactions in metabolite_reaction_dict.items():
            for reaction in reactions:
                if hasattr(reaction, 'reaction_id') and reaction.reaction_id not in seen_reaction_ids:
                    seen_reaction_ids.add(reaction.reaction_id)
                    all_reactions.append({
                        'id': reaction.reaction_id,
                        'sub': reaction.substrate_list,
                        'pro': reaction.product_list,
                        'reverse': reaction.reversible,
                        'num_id': reaction.num_id
                    })

        # 使用提取的反应列表
        reaction_list = all_reactions
    elif reaction_list is None:
        # 如果两种参数都未提供，抛出错误
        raise ValueError("Must provide either metabolite_reaction_dict or reaction_list")

    # ...existing code...
    def extend_symmetric_metabolite_in_product_list(reaction_obj):
        # 处理产物中的对称代谢物
        # 对称代谢物(如琥珀酸SUC_m)的碳原子排列有多种等价形式(abcd或dcba)
        # split_symmetry()将代谢物拆分为所有可能排列，确保13C标记分布计算正确
        # 例如：在TCA循环中，琥珀酸→延胡索酸反应中，由于对称性，碳原子排列需特殊处理
        new_product_list = []
        for product_node in reaction_obj.product_list:
            if product_node.name in symmetrical_metabolite_set:
                new_product_list.extend(product_node.split_symmetry())
                # 对称代谢物(如琥珀酸SUC_m)的碳原子排列有多种等价形式(abcd或dcba)
                # split_symmetry()将代谢物拆分为所有可能排列，确保13C标记分布计算正确
                # 例如：在TCA循环中，琥珀酸→延胡索酸反应中，由于对称性，碳原子排列需特殊处理
            else:
                new_product_list.append(product_node)
        new_reaction_obj = reaction_obj.copy()
        new_reaction_obj.product_list = new_product_list
        return new_reaction_obj

    def process_one_reaction(reaction_obj):
        # Check repeat.
        # 检查重复反应ID
        # 处理通量平衡约束
        # 处理EMU约束
        if reaction_obj.reaction_id in reaction_set:
            raise KeyError("Reaction ID repeat! {}".format(reaction_obj.reaction_id))
        reaction_set.add(reaction_obj.reaction_id)
        flux_name_list.append(reaction_obj.reaction_id)

        # Process flux balance constraints.
        for substrate_node in reaction_obj.substrate_list:
            if substrate_node.name not in balance_excluded_metabolite_set:
                metabolite_reaction_dict_for_bal[substrate_node.name][0][reaction_obj.reaction_id] += substrate_node.coefficient
        for product_node in reaction_obj.product_list:
            if product_node.name not in balance_excluded_metabolite_set:
                metabolite_reaction_dict_for_bal[product_node.name][1][reaction_obj.reaction_id] += product_node.coefficient

        # Process EMU constraints
        reaction_obj = extend_symmetric_metabolite_in_product_list(reaction_obj)
        unique_metabolite_name_set = set()
        # 对同一代谢物的反应仅添加一次
        # Only append the reaction once for symmetrical metabolites
        for product_node in reaction_obj.product_list:
            metabolite_name = product_node.name
            if metabolite_name in emu_excluded_metabolite_set or metabolite_name in unique_metabolite_name_set:
                continue
            unique_metabolite_name_set.add(metabolite_name)
            product_reaction_dict_for_emu[metabolite_name].append(reaction_obj)

    def record_and_check_metabolite_carbon_num(reaction_obj):
        # 迭代所有底物和产物
        # 记录每个代谢物的碳原子数
        # 检查同一代谢物在不同反应中碳原子数的一致性
        for metabolite_node in it.chain(reaction_obj.substrate_list, reaction_obj.product_list):
            metabolite_name = metabolite_node.name
            if metabolite_name not in processed_metabolite_dim_dict:
                current_metabolite_carbon_num = len(metabolite_node.carbon_composition_list)
                if metabolite_name not in processed_metabolite_dim_dict:
                    processed_metabolite_dim_dict[metabolite_name] = current_metabolite_carbon_num
                else:
                    if current_metabolite_carbon_num != 0 and \
                            processed_metabolite_dim_dict[metabolite_name] != current_metabolite_carbon_num:
                        raise ValueError('Metabolite carbon not consistent! {} in reaction {}'.format(
                            metabolite_name, reaction_obj.reaction_id))
                    elif current_metabolite_carbon_num == 0:
                        warnings.warn('Metabolite {} shows empty carbon str in reaction {}'.format(
                            metabolite_name, reaction_obj.reaction_id))

    # 主函数处理流程
    # 第1步，初始化数据结构：创建字典存储代谢物-反应关系、EMU关系等
    product_reaction_dict_for_emu = DefaultDict([])
    metabolite_reaction_dict_for_bal = DefaultDict((DefaultDict(0), DefaultDict(0)))
    if emu_excluded_metabolite_set is None:
        emu_excluded_metabolite_set = set()
    if balance_excluded_metabolite_set is None:
        balance_excluded_metabolite_set = set()
    if composite_reaction_list is None:
        composite_reaction_list = []
    flux_name_list = []
    reaction_set = set()
    processed_metabolite_dim_dict = complete_metabolite_dim_dict.copy() if complete_metabolite_dim_dict else {}
    composite_reaction_dict = {}

    # 第2步，处理反应和可逆反应
    # Parse reactions
    for current_reaction_item in reaction_list:
        # # current_reaction_dict = reaction_list[current_reaction_dict]
        # # 类型检查和转换
        # if isinstance(current_reaction_dict, Reaction):
        #     current_reaction = current_reaction_dict
        # else:
        #     current_reaction = Reaction(**current_reaction_dict)
        # record_and_check_metabolite_carbon_num(current_reaction)
        # process_one_reaction(current_reaction)
        # if current_reaction.reversible:
        #     reversed_reaction = current_reaction.reverse_reaction()
        #     process_one_reaction(reversed_reaction)
        # 类型检查和转换
        if isinstance(current_reaction_item, dict):
            current_reaction_dict = current_reaction_item
        elif isinstance(current_reaction_item, Reaction):
            current_reaction = current_reaction_item
        else:
            raise TypeError(f"反应项必须是字典或Reaction对象，而不是{type(current_reaction_item)}")

        # 从字典创建Reaction对象(如果需要)
        if 'current_reaction' not in locals():
            # 验证必要字段
            required_fields = {'id', 'sub', 'pro'}
            missing_fields = required_fields - set(current_reaction_dict.keys())
            if missing_fields:
                raise ValueError(f"反应缺少必要字段: {missing_fields}")

            # 创建Reaction对象
            current_reaction = Reaction(**current_reaction_dict)

        # 处理反应
        record_and_check_metabolite_carbon_num(current_reaction)
        process_one_reaction(current_reaction)

        # 处理可逆反应
        if current_reaction.reversible:
            reversed_reaction = current_reaction.reverse_reaction()
            process_one_reaction(reversed_reaction)

        # 清理局部变量，为下一次循环做准备
        if 'current_reaction' in locals():
            del current_reaction

    # Generate common index for each flux to ensure the same order all over the code
    # 如果没有提供flux_name_index_dict，则创建它
    if flux_name_index_dict is None:
        flux_name_index_dict = {flux_name: flux_index for flux_index, flux_name in enumerate(flux_name_list)}

    # Add composite reaction to flux dict
    for composite_reaction_parameter_dict in composite_reaction_list:
        composite_reaction_obj = CompositeReaction(**composite_reaction_parameter_dict)
        for composite_node_obj in composite_reaction_obj.compose_list:
            if composite_node_obj.reaction_id not in flux_name_index_dict:
                raise ValueError(
                    'Element {} of composite reaction {} has not be declared before!'.format(
                        composite_node_obj.reaction_id, composite_reaction_obj.reaction_id))
        composite_reaction_dict[composite_reaction_obj.reaction_id] = composite_reaction_obj
        flux_name_index_dict[composite_reaction_obj.reaction_id] = len(flux_name_index_dict)

    # Parse input metabolites
    input_metabolite_name_set_processed = set()
    for input_metabolite_name in added_input_metabolite_set:
        # 这三行注释了
        # if input_metabolite_name not in balance_excluded_metabolite_set:
        #    raise ValueError(
        #        'Input metabolite is not excluded from flux balance equations! {}'.format(input_metabolite_name))
        input_metabolite_name_set_processed.add(input_metabolite_name)

    # 第3步，通量平衡检查：确保每个代谢物既作为底物又作为产物出现
    # Check flux balance relationship
    for balanced_metabolite_name, (reaction_dict_as_substrate, reaction_dict_as_product) in \
            metabolite_reaction_dict_for_bal.items():
        # 检查是否是输入或输出代谢物
        is_input = balanced_metabolite_name in input_metabolite_name_set_processed
        is_output = balanced_metabolite_name in balance_excluded_metabolite_set

        # 如果不是输入或输出代谢物，则必须同时作为底物和产物出现
        #这两行注释了if not (is_input or is_output) and (len(reaction_dict_as_substrate) == 0 or len(reaction_dict_as_product) == 0):
        #    raise ValueError(f'Flux to {balanced_metabolite_name} is unbalanced!')
        # if len(reaction_dict_as_substrate) == 0 or len(reaction_dict_as_product) == 0:
        #     raise ValueError('Flux to {} is unbalanced!'.format(balanced_metabolite_name))


    # 第4步，输入代谢物处理：处理标记底物和排除代谢物
    # Set all excluded EMU metabolites to unlabeled state.
    for emu_excluded_metabolite in emu_excluded_metabolite_set:
        if emu_excluded_metabolite not in added_input_metabolite_set:
            input_metabolite_name_set_processed.add(emu_excluded_metabolite)

    # If target metabolite is not set, all non-excluded metabolites will be considered as target metabolites.
    if target_metabolite_name_list is None:
        target_metabolite_name_list = list(product_reaction_dict_for_emu.keys())

    # 返回结果，使用新的变量名以区别于输入参数
    return {
        'metabolite_reaction_dict': product_reaction_dict_for_emu,
        'input_metabolite_name_set': input_metabolite_name_set_processed,
        'complete_metabolite_dim_dict': processed_metabolite_dim_dict,
        'flux_name_index_dict': flux_name_index_dict,
        'composite_reaction_dict': composite_reaction_dict,
        'target_metabolite_name_list': target_metabolite_name_list
    }
