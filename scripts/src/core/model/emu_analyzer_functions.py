from ..common.packages import Counter, PriorityQueue
from ..common.classes import DefaultDict
from ..common.config import CoreConstants
from ..common.functions import check_if_biomass_flux

from .model_class import EMUElement
from ..common.classes import DictList
from ..common.config import CoreConstants


def decode_emu_name(emu_name):
    metabolite_name, carbon_list_str, *param_list = emu_name.split(CoreConstants.emu_carbon_list_str_sep)
    carbon_list = [int(char) for char in list(carbon_list_str)]
    return EMUElement(metabolite_name, carbon_list)

def find_overlapped_emu(current_emu, reaction_obj):
   # 内部实现
   # 函数的核心逻辑分为两部分：
   # find_overlapped_emu（主函数）
   #
   # 步骤：
   # 初始化最终结果列表 final_overlapped_emu_list。
   # 遍历反应中的产物节点，找到与 current_emu 代谢物名称匹配的产物。
   # 从 current_emu 中提取目标碳原子集合。
   # 调用 dfs_find_all_combination 获取所有底物 EMU 组合。
   # 将每个组合与对应产物的系数配对，添加到最终结果中。
   # 作用：协调搜索过程，并将结果与产物的化学计量信息结合。
    """
        该函数在反应对象 reaction_obj 中找到所有与当前 EMU (current_emu) 重叠的 EMU 组合。
        对于对称代谢物，它们在 reaction_obj 中已被扩展为两个系数为 0.5 的对称项。
        此函数已更新为考虑底物的系数，不再假设所有底物系数均为 1。

        :param current_emu: 当前的 EMU 对象，包含代谢物名称和选中的碳原子列表
        :param reaction_obj: 反应对象，包含底物列表和产物列表
        :return: 一个列表，包含元组 (产物系数, 重叠 EMU 组合)
    """
    """
    This function will find all combinations of overlapped EMU for current_emu in a reaction reaction_obj. For
    those symmetrical metabolites, they have been extended to two 0.5 symmetrical items in the reaction_obj.
    This function has been updated to consider the coefficients of substrates as well.

    :param current_emu:
    :param reaction_obj:
    :return:
    """

   #
   # dfs_find_all_combination（辅助函数）
   #
   # 这是一个深度优先搜索（DFS）函数，用于递归地搜索所有可能的底物 EMU 组合，以覆盖目标 EMU 的碳原子。
   #
   # 参数：
   # rest_carbon_atom_key_set：未映射的碳原子集合。
   # substrate_node_search_list：待搜索的底物节点列表。
   # current_combination：当前已找到的 EMU 组合。
   # result_combination_list：存储所有最终组合的列表。
   # 逻辑：
   # 如果所有碳原子都已映射（rest_carbon_atom_key_set 为空），将当前组合添加到结果中。
   # 如果没有更多底物可搜索（substrate_node_search_list 为空），停止搜索。
   # 否则，遍历底物列表：
   # 检查当前底物是否包含未映射的碳原子。
   # 如果包含，则创建一个新的 EMU（记录底物名称和选中的碳原子），更新未映射的碳原子集合和剩余底物列表，并递归调用自身。
   # 作用：通过递归探索所有可能的底物组合，确保找到所有覆盖目标碳原子的方案。
    def dfs_find_all_combination(
            rest_carbon_atom_key_set, substrate_node_search_list, current_combination, result_combination_list):
        """
                该函数使用深度优先搜索 (DFS) 找到所有能够生成特定产物 EMU 的底物 EMU 组合。
                现在考虑底物的系数，将其纳入EMU组合中。

                :param rest_carbon_atom_key_set: 未映射的碳原子集合，例如 {'a', 'b', 'c'}
                :param substrate_node_search_list: 需要搜索的底物节点列表，例如 [Node('AKG', 'abcde')]
                :param current_combination: 当前已找到的 EMU 组合，用于生成已映射的碳原子
                :param result_combination_list: 存储所有最终组合的列表
                :return: 无返回值，结果存储在 result_combination_list 中
        """
        """
        This function uses iterative DFS to search all EMU combinations that can generate a certain EMU
        of produce metabolite in a reaction_obj. Now considering the coefficients of substrates.

        :param rest_carbon_atom_key_set: Rest of carbon atom that has not been mapped, such as {'a', 'b', 'c'}
        :param substrate_node_search_list: List of reaction node that need to be search, such as [Node('AKG', 'abcde')]
        :param current_combination: Current EMU combination that can generate mapped carbon atoms.
        :param result_combination_list: Final list of all combinations.
        :return: None. The output is result_combination_list.
        """
        # 如果所有碳原子都已映射，则将当前组合添加到结果列表中
        if not rest_carbon_atom_key_set:
            result_combination_list.append(current_combination)
            return
        # 如果没有更多的底物结点可搜索，则停止
        elif not substrate_node_search_list:
            return
        else:
            # 遍历底物列表中的每个节点
            for substrate_index, substrate_node in enumerate(substrate_node_search_list):
                # 提取当前循环的底物名称，如"FUM"
                substrate_name = substrate_node.name
                # 提取当前循环的底物的碳原子组成列表，如['a','b','c','d']
                carbon_composition_list = substrate_node.carbon_composition_list
                # 获取底物系数(这一步骤和反应系数的采集和存储有关,属于系数处理的第1步)
                substrate_coefficient = substrate_node.coefficient
                # 初始化选中碳原子的标记列表，如果是4个碳的话就是[0,0,0,0]
                selected_carbon_list = [0] * len(carbon_composition_list)
                # 未映射碳原子的副本
                new_rest_carbon_set = set(rest_carbon_atom_key_set)
                # 检查当前底物中是否包含未映射的碳原子
                for carbon_atom_key in rest_carbon_atom_key_set:
                    try:
                        # 使用列表 list.index() 函数搜索当前底物结点的碳原子列表
                        # 判断 rest_carbon_atom_key_set (来源于reaction_id) 剩下的碳原子是否位于反应底物碳原子列表
                        # 查找碳原子位于 carbon_composition_list 的位置，如果未找到就抛出异常
                        loc = carbon_composition_list.index(carbon_atom_key)
                    except ValueError:
                        # 如果碳原子不在底物中，跳过
                        continue
                    else:
                        # 如果在当前底物结点找到了 rest_carbon_atom_key_set 标记的碳原子
                        # 将 selected_carbon_list 相应位置的 0 改为 1
                        selected_carbon_list[loc] = 1
                        # 更新 rest_carbon_atom_key_set 作为 new_rest_carbon_set
                        # 具体做法是减掉当前的 carbon_atom_key 字母
                        new_rest_carbon_set -= {carbon_atom_key}
                # 如果选中碳原子的标记列表 selected_carbon_list 不全为0
                if sum(selected_carbon_list) > 0:
                    # 创建新的 EMU 字典, 记录底物名称、选中的碳原子和系数 (这一步包含底物系数)
                    new_emu_dict = {
                        'metabolite_name': substrate_name,
                        'selected_carbon_list': selected_carbon_list,
                        'coefficient': substrate_coefficient,  # 确保保存底物系数
                    }
                    # 更新搜索空间为当前底物之后的底物列表
                    new_search_space_list = substrate_node_search_list[substrate_index + 1:]
                    # 将新 EMU 添加到当前组合中
                    new_combination = current_combination + [new_emu_dict]
                    # 递归搜索剩余的碳原子
                    dfs_find_all_combination(
                        new_rest_carbon_set, new_search_space_list, new_combination, result_combination_list)

    # 初始化存储最终重叠 EMU 组合的列表
    final_overlapped_emu_list = []
    # 遍历反应中的所有产物节点
    for product_node in reaction_obj.product_list:
        # 这个集合推导式受到 current_emu.metabolite_name 的筛选
        # 所以它只会针对匹配的产物节点执行（根据结点名称匹配判断，如MAL）
        if product_node.name == current_emu.metabolite_name:
            # 提取当前选中的Product EMU节点碳原子集合
            # 假设Product_list 有2个节点
            # 每个节点的碳原子列表相同均为 carbon_composition_list = ['a','b','c','d']
            # enumerate(product_node.carbon_composition_list) 会生成:
            # [(0, 'a'), (1, 'b'), (2, 'c'), (3, 'd')]
            # 最终的结果集取决于 current_emu.selected_carbon_list 的值：
            # 如果 selected_carbon_list = [1,1,1,1]，结果将是: {'a','b','c','d'}
            # 如果 selected_carbon_list = [1,0,1,0]，结果将是: {'a','c'}
            # 如果 selected_carbon_list = [0,0,0,0]，结果将是空集: set()
            # 此处 selected_carbon_list = [0,1,0,1]，结果是：{'b','d'}
            query_carbon_atom_key_set = {
                carbon_atom_key for index, carbon_atom_key in enumerate(product_node.carbon_composition_list)
                if current_emu.selected_carbon_list[index]}
            # 初始化存储不同组合的列表
            different_combination_list = []
            # 调用 DFS 函数，搜索所有可能的底物 EMU 组合
            dfs_find_all_combination(
                query_carbon_atom_key_set, reaction_obj.substrate_list, [], different_combination_list)
            # 将每个组合与产物系数一起添加到最终列表中
            for combination in different_combination_list:
                # 检查组合中是否存在底物系数，如果存在则考虑底物系数与产物系数的关系
                # 对于包含多个EMU的组合，需要考虑所有底物的系数
                final_overlapped_emu_list.append((product_node.coefficient, combination))
    ##################################################################
    # # 检查并确保所有底物系数都是有效的
    # for i, (product_coefficient, combination) in enumerate(final_overlapped_emu_list):
    #     for j, emu_dict in enumerate(combination):
    #         # 确保系数存在且有效
    #         if 'coefficient' not in emu_dict or emu_dict['coefficient'] is None or emu_dict['coefficient'] <= 0:
    #             # 将缺失或无效的系数设置为默认值1.0
    #             emu_dict['coefficient'] = 1.0
    #             # 更新组合中的字典
    #             final_overlapped_emu_list[i][1][j] = emu_dict
    ##################################################################
    
    # 返回所有重叠 EMU 组合
    # 总结
    # find_overlapped_emu 是一个基于 DFS 的搜索算法，用于在代谢反应中找到与目标 EMU 重叠的所有底物 EMU 组合。
    # 它通过递归探索底物中的碳原子映射，确保全面覆盖所有可能性，并保留产物和底物的化学计量信息，为代谢通量分析提供关键数据支持。
    return final_overlapped_emu_list

# 这个函数涉及系数处理的第3步，累积所有反应底物的系数
def emu_equation_analyzer(
    metabolite_reaction_dict,   # 代谢物作为底物，参与反应的映射字典
    target_metabolite_name_list,    # 目标代谢物名称的列表
    input_metabolite_name_set,  # 输入代谢物的名称集合
    complete_metabolite_dim_dict    # 代谢物到其碳原子数量的映射字典
):
    """分析并生成 EMU 方程，考虑底物系数
        其他部分代码保持不变
        在处理重叠 EMU 时，考虑系数：
        overlapped_emu_list = find_overlapped_emu(current_emu, reaction_obj)
        for emu_obj, coefficient in overlapped_emu_list:
            # 在生成方程时使用系数
            # 例如可以将系数存储在方程对象的属性中，或者直接应用于方程系数矩阵
            # ...
    """

    # 内部函数 1：add_new_emu_to_queue
    # 功能：
    # 将一个新的 EMU 添加到优先队列中，确保碳原子数量较少的 EMU 优先被处理。
    def add_new_emu_to_queue(current_emu: EMUElement):
        # 记录添加到队列的 EMU 数量，每次调用时递增。
        insert_count['count'] += 1
        # 将 EMU 加入优先队列 emu_need_to_add
        # 优先级由元组 (-1 * current_emu.emu_carbon_num, insert_count['count']) 决定：
        # -1 * current_emu.emu_carbon_num：负值使碳原子数量多的 EMU 优先级更高。
        # insert_count['count']：插入顺序作为次级排序条件，确保同碳原子数量的 EMU 按插入顺序处理。
        emu_need_to_add.put(((-1 * current_emu.emu_carbon_num, insert_count['count']), current_emu))

    # 内部函数：sort_emu_equations_dict_by_carbon_num
    # 功能：
    # 对 EMU 方程字典按碳原子数量进行升序排序。
    def sort_emu_equations_dict_by_carbon_num(raw_emu_equations_dict):
        # 提取原始字典的键（碳原子数量）并排序。
        sorted_key_list = sorted(list(raw_emu_equations_dict.keys()))
        # 新建字典，按排序后的键重新组织原始字典内容。
        sorted_mid_equations_dict = {}
        # 然后按碳数升序对其排序。遍历排序后的键列表，将原始字典中的方程按照碳数排序后添加到新的字典中。
        for sorted_key in sorted_key_list:
            sorted_mid_equations_dict[sorted_key] = raw_emu_equations_dict[sorted_key]
        # 返回值：按碳原子数量升序排序的 EMU 方程字典。
        return sorted_mid_equations_dict

    # Initialize dictionaries and counters for EMU equations and input EMUs
    # emu_mid_equations_dict_carbon_num_list = defaultdict(
    #     lambda: defaultdict(lambda: []))

    # 初始化数据结构：
    # 功能：
    # 按碳原子数量分组存储 EMU 方程, 这些变量用于管理和跟踪EMU的状态和队列。
    emu_mid_equations_dict_carbon_num_list = DefaultDict(DefaultDict([]))
    # 记录 EMU 访问次数
    complete_emu_visiting_num_dict = Counter()
    # 存储输入代谢物的 EMU
    input_emu_dict = {}
    # 跟踪插入队列的 EMU 数量
    insert_count = {'count': 0}
    # 优先队列，按碳原子数量处理 EMU
    emu_need_to_add = PriorityQueue()

    # 处理目标代谢物：
    # 功能：
    # 为每个目标代谢物创建 EMU，并根据其是否为输入代谢物进行分类处理。
    for target_metabolite_name in target_metabolite_name_list:
        # 从 complete_metabolite_dim_dict 获取目标代谢物的碳原子数量。
        carbon_num = complete_metabolite_dim_dict[target_metabolite_name]
        # 创建 EMU 对象，选中所有碳原子（用 [1] * carbon_num 表示）。
        target_emu = EMUElement(target_metabolite_name, [1] * carbon_num)
        # 记录该 EMU 被访问。
        complete_emu_visiting_num_dict[target_emu.emu_name] += 1
        # 如果目标代谢物 EMU 是输入代谢物（in input_metabolite_name_set），添加到input_emu_dict。
        if target_emu.metabolite_name in input_metabolite_name_set:
            input_emu_dict[target_emu.emu_name] = target_emu
        else:
            # 否则，调用 add_new_emu_to_queue 将目标 EMU 加入优先队列
            add_new_emu_to_queue(target_emu)

    # Process the EMUs in the queue.
    # 处理队列中的 EMU
    # 功能：处理优先队列中的EMU元素，处理优先队列中的 EMU，生成相关的 EMU 方程。
    # 队列循环：
    # emu_need_to_add.get() 取出优先级最高的 EMU（碳原子数量最少）。
    while not emu_need_to_add.empty():
        _, analyzing_emu = emu_need_to_add.get()
        # 反应遍历：
        # 从 metabolite_reaction_dict 获取当前 EMU 参与的反应。
        # check_if_biomass_flux：跳过生物质生成反应。
        for reaction in metabolite_reaction_dict[analyzing_emu.metabolite_name]:
            reaction_id = reaction.reaction_id
            # if reaction_id == CoreConstants.biomass_flux_id:
            # 对每个反应，检查是否与生物质生成反应相关（check_if_biomass_flux）。如果是，则跳过。
            if check_if_biomass_flux(reaction_id):
                continue

            overlapped_emu_dict_combination_list = find_overlapped_emu(analyzing_emu, reaction)
            # 重叠 EMU 处理：
            # find_overlapped_emu：查找与当前 EMU 重叠的底物 EMU 组合。
            # 遍历每个重叠组合（coefficient, emu_dict_combination）：
            for product_coefficient, emu_dict_combination in overlapped_emu_dict_combination_list:
                emu_combination_list = []
                # 添加一个标志，用于检查是否已经处理了底物系数
                all_substrate_coefficients_valid = True
                # 用于累积所有底物的系数乘积
                combined_substrate_coefficient = 1.0
                
                for overlapped_emu_dict in emu_dict_combination:
                    metabolite_name = overlapped_emu_dict['metabolite_name']
                    selected_carbon_list = overlapped_emu_dict['selected_carbon_list']
                    # 获取底物系数，如果不存在则使用默认值1.0
                    substrate_coefficient = overlapped_emu_dict.get('coefficient', 1.0)
                    # 验证系数是否有效
                    if substrate_coefficient is None or substrate_coefficient <= 0:
                        all_substrate_coefficients_valid = False
                        substrate_coefficient = 1.0  # 使用默认值
                    # 累积系数乘积
                    combined_substrate_coefficient *= substrate_coefficient
            
                    # 构建EMU名称
                    # 构建 overlapped_emu_name：格式为 代谢物名_碳原子列表。
                    overlapped_emu_name = "{}{}{}".format(
                        metabolite_name, CoreConstants.emu_carbon_list_str_sep, "".join(
                            [str(num) for num in selected_carbon_list]))
            
                    # 创建新的EMU元素并保存底物系数
                    # If the overlapped EMU is in input metabolite, add it to input EMU list.
                    # 条件分支：
                    # 如果是输入代谢物，添加到 input_emu_dict（避免重复添加）。
                    if metabolite_name in input_metabolite_name_set:
                        new_emu_element = EMUElement(metabolite_name, selected_carbon_list)
                        # 明确设置系数属性，确保不为None
                        new_emu_element.coefficient = substrate_coefficient
                        if overlapped_emu_name not in input_emu_dict:
                            input_emu_dict[overlapped_emu_name] = new_emu_element
                    # 如果不是输入代谢物且未被访问，添加到队列并记录访问。
                    # If the overlapped EMU is not input metabolite, check if it is visited in the path of current
                    # layer.
                    else:
                        new_emu_element = EMUElement(metabolite_name, selected_carbon_list)
                        # 明确设置系数属性，确保不为None
                        new_emu_element.coefficient = substrate_coefficient
            
                        if overlapped_emu_name not in complete_emu_visiting_num_dict:
                            complete_emu_visiting_num_dict[overlapped_emu_name] += 1
                            add_new_emu_to_queue(new_emu_element)
            
                    # 将新 EMU 添加到 emu_combination_list。
                    emu_combination_list.append(new_emu_element)
            
                # 这里是反应底物系数确定的第4步, 确定最终系数：产物系数与所有底物系数的组合
                final_coefficient = product_coefficient
                # 如果所有底物系数有效，应用组合系数
                if all_substrate_coefficients_valid and len(emu_combination_list) > 0:
                    # 对于多个底物的情况，考虑所有底物系数的影响
                    final_coefficient *= combined_substrate_coefficient
            
                # 生成 EMU 方程：
                # 创建元组 (reaction_id, emu_combination_list, final_coefficient)
                emu_tuple = (reaction_id, emu_combination_list, final_coefficient)
                
                # 按碳原子数量和 EMU 全名存储到 emu_mid_equations_dict_carbon_num_list
                emu_mid_equations_dict_carbon_num_list[analyzing_emu.emu_carbon_num][analyzing_emu.full_name].append(
                    emu_tuple)

    # 排序并返回结果
    # 功能：
    # 对 EMU 方程字典按碳原子数量排序，并返回结果。
    # Sort the EMU equations dictionary by carbon number and return the results.

    # 调用 sort_emu_equations_dict_by_carbon_num 完成排序。
    sorted_emu_equations_dict_carbon_num_list = sort_emu_equations_dict_by_carbon_num(
        emu_mid_equations_dict_carbon_num_list)
    # 返回排序后的 EMU 方程字典和输入 EMU 字典。
    # 在案例中，input_emu_dict 指的是 A 结点，sorted_emu_equations_dict_carbon_num_list 指的是 EMU 列表
    return sorted_emu_equations_dict_carbon_num_list, input_emu_dict
    # emu_mid_equation_dict 作为返回值接收 sorted_emu_equations_dict_carbon_num_list
    # emu_mid_equation_dict, input_emu_dict 作为实参输入 emu_matrix_equation_generator()

    # 总结：
    # emu_equation_analyzer 函数通过以下步骤实现其功能：
    #
    # 初始化数据结构，管理 EMU 的处理和记录。
    # 将目标代谢物的 EMU 添加到队列或输入字典。
    # 使用优先队列按碳原子数量处理 EMU，生成相关的方程。
    # 按碳原子数量排序并返回结果。
    # 该函数在代谢通量分析中用于追踪碳原子流动，特别适用于大规模代谢网络的建模和分析。

# This function analyzes the dependencies between EMUs in metabolic reactions.
# It generates dictionaries for EMU dependencies, EMU object indices, and input EMUs.
# It returns the EMU dependency dictionary, the complete EMU name object index dictionary, and the input EMU dictionary.
def emu_dependency_analyzer(metabolite_reaction_dict, input_metabolite_name_set, complete_metabolite_dim_dict):
    emu_name_dependency_dict = {}
    complete_emu_name_obj_index_dict = {}
    input_emu_dict = {}
    processed_emu_name_dict = {}
    unprocessed_emu_obj_dict = {}

    for metabolite_name in metabolite_reaction_dict.keys():
        carbon_num = complete_metabolite_dim_dict[metabolite_name]
        current_emu_obj = EMUElement(metabolite_name, [1] * carbon_num)
        current_emu_name = current_emu_obj.full_name
        if current_emu_name not in complete_emu_name_obj_index_dict:
            complete_emu_name_obj_index_dict[current_emu_name] = (
                current_emu_obj, len(complete_emu_name_obj_index_dict))
        unprocessed_emu_obj_dict[current_emu_obj] = None

    while len(unprocessed_emu_obj_dict) > 0:
        current_emu_obj = unprocessed_emu_obj_dict.keys().__iter__().__next__()
        current_emu_full_name = current_emu_obj.full_name
        if current_emu_full_name not in emu_name_dependency_dict:
            emu_name_dependency_dict[current_emu_full_name] = {}
        for reaction in metabolite_reaction_dict[current_emu_obj.metabolite_name]:
            reaction_id = reaction.reaction_id
            # if reaction_id == CoreConstants.biomass_flux_id:
            if check_if_biomass_flux(reaction_id):
                continue
            overlapped_emu_dict_combination_list = find_overlapped_emu(current_emu_obj, reaction)
            """
                Each item in overlapped_emu_dict_combination_list reflects appearance of one target metabolite 
                in current reaction.
                Usually len(overlapped_emu_dict_combination_list) == 1 since one metabolite usually appear once 
                in a reaction. The number > 1 appears when a reaction include metabolite with stoichiometric number
                > 1. In this case, if one EMU depends on same EMU more than once, their coefficients need to be summed.
                If depends on different EMUs, all of them need to be recorded
            """
            for coefficient, emu_dict_combination in overlapped_emu_dict_combination_list:
                dependent_emu_list = []
                for overlapped_emu_dict in emu_dict_combination:
                    metabolite_name = overlapped_emu_dict['metabolite_name']
                    selected_carbon_list = overlapped_emu_dict['selected_carbon_list']
                    # overlapped_emu_name = "{}{}{}".format(
                    #     metabolite_name, CoreConstants.emu_carbon_list_str_sep, "".join(
                    #         [str(num) for num in selected_carbon_list]))
                    overlapped_emu_obj = EMUElement(metabolite_name, selected_carbon_list)
                    overlapped_emu_name = overlapped_emu_obj.full_name
                    if metabolite_name in input_metabolite_name_set:
                        if overlapped_emu_name not in input_emu_dict:
                            input_emu_dict[overlapped_emu_name] = overlapped_emu_obj
                    elif overlapped_emu_name not in processed_emu_name_dict:
                        unprocessed_emu_obj_dict[overlapped_emu_obj] = None
                    if overlapped_emu_name not in complete_emu_name_obj_index_dict:
                        complete_emu_name_obj_index_dict[overlapped_emu_name] = (
                            overlapped_emu_obj, len(complete_emu_name_obj_index_dict))
                    dependent_emu_list.append(overlapped_emu_obj)
                if len(dependent_emu_list) > 1:
                    sorted_dependent_emu_list = sorted(dependent_emu_list, key=lambda x: x.full_name)
                    sorted_dependent_emu_name_list = [
                        dependent_emu.full_name for dependent_emu in sorted_dependent_emu_list]
                    convolution_emu_obj = current_emu_obj.copy_to_convolution(sorted_dependent_emu_list)
                    convolution_emu_full_name = convolution_emu_obj.full_name
                    if convolution_emu_full_name not in complete_emu_name_obj_index_dict:
                        complete_emu_name_obj_index_dict[convolution_emu_full_name] = (
                            convolution_emu_obj, len(complete_emu_name_obj_index_dict))
                    if convolution_emu_full_name not in emu_name_dependency_dict:
                        emu_name_dependency_dict[convolution_emu_full_name] = {
                            convoluted_emu_name: [(CoreConstants.convolution_id, 1)]
                            for convoluted_emu_name in sorted_dependent_emu_name_list}
                    dependent_emu_name = convolution_emu_full_name
                else:
                    the_only_dependent_emu = dependent_emu_list[0]
                    dependent_emu_name = the_only_dependent_emu.full_name
                if current_emu_full_name not in emu_name_dependency_dict:
                    emu_name_dependency_dict[current_emu_full_name] = {}
                if dependent_emu_name not in emu_name_dependency_dict[current_emu_full_name]:
                    emu_name_dependency_dict[current_emu_full_name][dependent_emu_name] = []
                emu_name_dependency_dict[current_emu_full_name][dependent_emu_name].append((reaction_id, coefficient))

        processed_emu_name_dict[current_emu_full_name] = None
        del unprocessed_emu_obj_dict[current_emu_obj]

    return emu_name_dependency_dict, complete_emu_name_obj_index_dict, input_emu_dict

# 总结：
# 这段代码的核心功能是生成代谢流分析中EMU矩阵方程，用于描述代谢物在代谢网络中的流动。
# 通过解码代谢物信息、组合不同的EMU和反应信息，最终生成符合特定条件的矩阵方程。
# 这些矩阵将用于后续的代谢流计算和优化分析。
# 把抽象的EMU关系变为矩阵定义
# 总结：这段代码的核心功能是生成代谢流分析中EMU矩阵方程，用于描述代谢物在代谢网络中的流动。
# 通过解码代谢物信息、组合不同的EMU和反应信息，最终生成符合特定条件的矩阵方程。
# 这些矩阵将用于后续的代谢流计算和优化分析。
# 参数说明：
# emu_mid_equations_dict_carbon_num_list：按碳数分组的代谢物方程字典，每个代谢物有其相关的代谢方程。
# input_emu_dict：包含输入代谢物（EMU）的字典。
# This function generates matrix equations for EMU analysis.
# Parameters:
# emu_mid_equations_dict_carbon_num_list: A dictionary of EMU equations sorted by carbon number.
# input_emu_dict: A dictionary of input EMUs.
# Returns: A dictionary of matrix equations for EMU analysis.
# 代码的主要功能是将代谢通量分析中的 EMU（Elementary Metabolite Unit）方程按碳原子数量转换为矩阵形式，用于后续的线性代数求解。
def emu_matrix_equation_generator(
        emu_mid_equations_dict_carbon_num_list,
        input_emu_dict):
    """
        将EMU方程转换为矩阵形式,用于代谢通量分析。正确考虑底物系数
        输入的EMU方程被转换为线性方程组 A * x = B * y 的矩阵形式，其中：
        - x 表示当前层的代谢物同位素分布
        - y 表示输入或低层次的EMU
        - A 和 B 是系数矩阵

        参数：
            emu_mid_equations_dict_carbon_num_list: 嵌套字典,按碳原子数量分组的EMU方程
            input_emu_dict: 输入EMU信息的字典

        返回：
            emu_matrix_equation_dict: 按碳原子数量存储矩阵方程信息的字典
    """
    # 初始化一个字典，存储所有EMU对象的解码信息
    # 复制输入的input_emu_dict作为初始值
    complete_emu_object_dict = dict(input_emu_dict)

    # 初始化结果字典，用于按碳原子数量存储矩阵方程信息
    emu_matrix_equation_dict = {}

    # 遍历每个碳原子数量及其对应的EMU方程字典
    for carbon_num, emu_mid_equations_dict in emu_mid_equations_dict_carbon_num_list.items():
        # 创建当前碳数层 (carbon_num) 的EMU字典列表，使用DictList()以支持键值访问和索引
        this_layer_emu_dict_list = DictList()

        # 遍历当前碳数下的所有EMU名称
        for emu_name in emu_mid_equations_dict.keys():
            # 如果EMU不在complete_emu_object_dict中，解码并添加
            if emu_name not in complete_emu_object_dict:
                complete_emu_object_dict[emu_name] = decode_emu_name(emu_name)
            # 将EMU对象添加到当前层的字典列表中
            this_layer_emu_dict_list[emu_name] = complete_emu_object_dict[emu_name]

        # 初始化输入和低层次EMU的字典列表
        # A * x = b * y
        input_and_lower_layer_emu_dict_list = DictList()

        # 初始化矩阵A和矩阵B的系数位置字典，使用嵌套默认字典，默认值为0
        matrix_a_flux_location_dict = DefaultDict(DefaultDict(0))
        matrix_b_flux_location_dict = DefaultDict(DefaultDict(0))

        # 遍历当前碳原子数层 emu_mid_equations_dict 中的每个EMU
        # analyzed_emu_name：当前分析的EMU名称
        # mid_equations_list：对应的生成方程列表
        for analyzed_emu_name, mid_equations_list in emu_mid_equations_dict.items():
            # 获取行索引：获取正在分析的 EMU analyzed_emu_index 在本层 EMU 列表 this_layer_emu_dict_list 中的索引
            # 这个索引将成为矩阵 A 和 B 中的行索引，对应 A*x = B*y 方程组中的一个方程
            analyzed_emu_index = this_layer_emu_dict_list.index(analyzed_emu_name)

            # 内层循环：遍历构成 analyzed_emu_name 的每一个具体反应步骤 (输入的前置反应底物)
            # reaction_id：反应过程的 ID (例如v0到v6)
            # emu_combination_list：与该反应 ID 所对应反应的底物 EMU 组合 (EMUElement 对象列表)
            # param_list：反应底物所对应的系数列表 (coefficient)
            for reaction_id, emu_combination_list, *param_list in mid_equations_list:
                # flux_value_tensor = flux_tensor[flux_name_index_dict[reaction_id]] * param_list[0]
                # 获取反应底物系数。如果提供了参数 param_list，则使用第一个参数作为系数，否则默认为1.0
                # 这是反应底物系数处理的第5步, 系数整合到矩阵方程中:
                if len(param_list) == 1:
                    coefficient = param_list[0]
                else:
                    coefficient = 1.0
                # 处理源反应底物 EMU emu_combination_list
                # Case 1: 只有一个底物 EMU
                if len(emu_combination_list) == 1:
                    # 获取唯一的源 EMU 对象
                    emu_object = emu_combination_list[0]
                    # 获取源 EMU 的名称
                    dependent_emu_name = emu_object.emu_name
                    # Subcase 1.1: 如果源 EMU 在当前层中, 这意味着源 EMU 也是 X 向量(矩阵) 的一部分
                    if dependent_emu_name in this_layer_emu_dict_list:
                        # 获取源 EMU 在 x 向量中的索引，这将是矩阵 A 中的列索引
                        current_layer_emu_col_index = this_layer_emu_dict_list.index(dependent_emu_name)
                        # 在 matrix_a_flux_location_dict 中记录: 在矩阵 A 的 (analyzed_emu_index, current_layer_emu_col_index) 位置
                        # 反应 reaction_id 贡献了 + coefficient 的值, 这表明 X 向量内部的依赖关系 (非对角线元素)
                        # A 矩阵的系数都是正数, 除了对角线是负数, 所以这里的系数是相加
                        # 在矩阵A和B中应用反应系数
                        matrix_a_flux_location_dict[
                            (analyzed_emu_index, current_layer_emu_col_index)][reaction_id] += coefficient
                    # Subcase 1.2: 如果源 EMU 在输入和低层次 EMU 列表中, 这意味着源 EMU 是 Y 向量(矩阵) 的一部分
                    elif dependent_emu_name in complete_emu_object_dict:
                        # 这意味着源 EMU 是 Y 向量(矩阵) 的一部分
                        input_and_lower_layer_emu_dict_list[dependent_emu_name] = emu_object
                        # 确保源 EMU 在输入和低层次 EMU 列表中, 获取源 EMU 在 Y 向量中的索引, 这将是矩阵 B 中的列索引
                        current_layer_input_col_index = input_and_lower_layer_emu_dict_list.index(dependent_emu_name)
                        # 在 matrix_b_flux_location_dict 中记录:
                        # 在矩阵 B 的 (analyzed_emu_index, current_layer_input_col_index) 位置
                        # 反应 reaction_id 贡献了 -coefficient 的值, 注意这里的负号, 这表明 Y 向量(矩阵) 内部的依赖关系
                        matrix_b_flux_location_dict[
                            (analyzed_emu_index, current_layer_input_col_index)][reaction_id] -= coefficient
                    # Subcase 1.3: 源 EMU 未知, 如果源 EMU 不在当前层，也不在输入和低层次 EMU 列表中, 表明目前的源 EMU 数据有问题
                    else:
                        raise ValueError()
                # Case 2: 多个底物 EMU/ 源 EMU 作为反应底物组合 (后续用于卷积计算)
                # 这种反应多发生在 TCA 循环中, 例如: OAA + Acetyl-CoA -> Citrate, 两个较小的分子结合形成一个较大的分子，他们的同位素模式将发生卷积
                else:
                    # 检查所有参与卷积的 EMU 是否来源于输入或更低的层(即属于 Y 向量/矩阵)
                    for emu_object in emu_combination_list:
                        dependent_emu_name = emu_object.emu_name
                        if dependent_emu_name not in complete_emu_object_dict:
                            raise ValueError()
                    # 对 EMU 列表进行排序, 以确保组合名称的唯一性
                    sorted_emu_combination_list = sorted(emu_combination_list)
                    # 创建一个可以代表该 EMU 组合的唯一名称
                    dependent_emu_complete_name = CoreConstants.convolution_emu_sep.join(
                        [emu.full_name for emu in sorted_emu_combination_list])
                    # 如果该组合是全新的, 则添加进入代表 Y 向量的列表中
                    if dependent_emu_complete_name not in input_and_lower_layer_emu_dict_list:
                        input_and_lower_layer_emu_dict_list[dependent_emu_complete_name] = sorted_emu_combination_list
                    # 获取这个组合在 Y 向量中的索引，作为矩阵 B 的列索引
                    current_layer_input_col_index = input_and_lower_layer_emu_dict_list.index(
                        dependent_emu_complete_name)
                    # 在矩阵 B (matrix_b_flux_location_dict) 中进行记录
                    # 矩阵 B 的 (analyzed_emu_index, current_layer_input_col_index) 位置, 反应 reaction_id 贡献了 -coefficient
                    # 这和单个输入/底层 EMU 的情况比较相似
                    matrix_b_flux_location_dict[
                        (analyzed_emu_index, current_layer_input_col_index)][reaction_id] -= coefficient
                # 最后处理矩阵A 的对角线元素
                # 非常重要: 对于每一个对 analyzed_emu_name 的生成有贡献的反应 (reaction_id)，无论其来源如何，都会在 matrix_a_flux_location_dict 中记录
                # 在矩阵 A 的对角线位置, 反应 贡献了 -coefficient
                # 对角线元素通常代表了与 analyzed_emu_name 自身相关的项，例如总生成速率或消耗速率（取决于方程的排列方式）
                # 这里的负号表明它可能代表移到等式左侧的生成项
                matrix_a_flux_location_dict[
                    (analyzed_emu_index, analyzed_emu_index)][reaction_id] -= coefficient

        # 构建最终的矩阵方程
        # 最终生成的 matrix_a_flux_location_dict 和 matrix_b_flux_location_dict 存储了构建稀疏矩阵 A 和 B 所需的所有信息
        # 包含行索引、列索引、关联的反应 ID 和 系数 coefficient
        matrix_a_dim = len(emu_mid_equations_dict)
        matrix_b_col = len(input_and_lower_layer_emu_dict_list)
        emu_matrix_equation_dict[carbon_num] = (
            this_layer_emu_dict_list,
            input_and_lower_layer_emu_dict_list,
            matrix_a_flux_location_dict,
            matrix_b_flux_location_dict,
            matrix_a_dim,
            matrix_b_col
        )

    return emu_matrix_equation_dict