import numpy as np

from scripts.src.core.common.config import CoreConstants
from scripts.src.core.model.model_class import Reaction, Node
from common_test_pipeline import run_mfa_pipeline
from common_test_pipeline import load_hct116_experimental_data
from common_test_pipeline import run_mfa_pipeline_with_experimental_data


def setup_tca_network():
    """Setup TCA cycle metabolic network"""
    # 定义反应方程，代谢物的碳原子数量，以及反应对应的代号

    # 先定义代谢反应方程
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

    # 反应3: AKG → Glu (碳原子转换: abcde → abcde)
    v3 = Reaction(id='v3',
                     sub=[Node(name='AKG', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd', 'e'])],
                     pro=[Node(name='Glu', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd', 'e'])],
                     reverse=False)

    # 反应4: AKG → Suc + CO2 (碳原子转换: abcde → bcde + a)
    v4 = Reaction(id='v4',
                     sub=[Node(name='AKG', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd', 'e'])],
                     pro=[Node(name='Suc', coefficient=1, carbon_composition_string=['b', 'c', 'd', 'e']),
                          Node(name='CO2', coefficient=1, carbon_composition_string=['a'])],
                     reverse=False)

    # 反应5: Suc → Fum (考虑碳原子对称性)
    # 对称分子 (如琥珀酸转化为延胡索酸) 通过设置多个产物节点并分配合适的系数来处理,这是反应底物系数处理的第2步:
    # Suc(abcd) -> Fum(abcd) 和 Suc(abcd) -> Fum(dcba) 各占 50%
    v5 = Reaction(id='v5',
                     sub=[Node(name='Suc', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd'])],
                     pro=[Node(name='Fum', coefficient=0.5, carbon_composition_string=['a', 'b', 'c', 'd']),
                          Node(name='Fum', coefficient=0.5, carbon_composition_string=['d', 'c', 'b', 'a'])],
                     reverse=False)

    # 反应6: Fum → OAC (考虑碳原子对称性)
    # Fum(abcd) -> OAC(abcd) 和 Fum(dcba) -> OAC(abcd) 各占 50%
    # 需要定义两个Fum输入节点，代表其对称性
    v6 = Reaction(id='v6',
                     sub=[Node(name='Fum', coefficient=0.5, carbon_composition_string=['a', 'b', 'c', 'd']),
                          Node(name='Fum', coefficient=0.5, carbon_composition_string=['d', 'c', 'b', 'a'])],
                     pro=[Node(name='OAC', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd'])],
                     reverse=False)

    # 反应7: OAC → Fum (考虑碳原子对称性)
    # OAC(abcd) -> Fum(abcd) 和 OAC(abcd) -> Fum(dcba) 各占 50%
    v7 = Reaction(id='v7',
                     sub=[Node(name='OAC', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd'])],
                     pro=[Node(name='Fum', coefficient=0.5, carbon_composition_string=['a', 'b', 'c', 'd']),
                          Node(name='Fum', coefficient=0.5, carbon_composition_string=['d', 'c', 'b', 'a'])],
                     reverse=False)

    # 反应8: Asp → OAC (abcd → abcd)
    v8 = Reaction(id='v8',
                     sub=[Node(name='Asp', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd'])],
                     pro=[Node(name='OAC', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd'])],
                     reverse=False)


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

    # 用于反应分析的代谢物对应反应方程列表
    flux_balance_reaction_dict = {
        # 输入代谢物
        'Asp': [v8],
        'AcCoA': [v1],

        # TCA循环代谢物
        'Cit': [v1, v2],
        'AKG': [v2, v3, v4],
        'Suc': [v4, v5],
        'Fum': [v5, v7, v6],
        'OAC': [v6, v8, v1, v7],

        # 输出代谢物
        'Glu': [v3],
        'CO2': [v2, v4]
    }

    # 定义代谢物碳原子数
    complete_metabolite_dim_dict = {
        # 输入代谢物
        'Asp': 4,       # 天冬氨酸 (4C)
        'AcCoA': 2,     # 乙酰CoA (2C)

        # TCA循环代谢物
        'Cit': 6,       # 柠檬酸 (6C)
        'AKG': 5,       # α-酮戊二酸 (5C)
        'Suc': 4,       # 琥珀酸 (4C)
        'Fum': 4,       # 延胡索酸 (4C)
        'OAC': 4,       # 草酰乙酸 (4C)

        # 输出代谢物
        'Glu': 5,       # 谷氨酸 (5C)
        #'CO2': 1        # 二氧化碳 (1C)
    }

    # Define flux indices
    flux_name_index_dict = {
        'v1': 0,  # OAC + AcCoA → Cit (100)
        'v2': 1,  # Cit → AKG + CO2 (100)
        'v3': 2,  # AKG → Glu (50)
        'v4': 3,  # AKG → Suc + CO2 (50)
        'v5': 4,  # Suc → Fum (50)
        'v6': 5,  # Fum → OAC (125)
        'v7': 6,  # OAC → Fum (75)
        'v8': 7   # Asp → OAC (50)
    }

    return metabolite_reaction_dict, flux_balance_reaction_dict, complete_metabolite_dim_dict, flux_name_index_dict


def test_tca_cycle_mfa():
    """使用实验数据测试TCA循环的代谢通量分析"""
    print("启动TCA循环代谢通量分析与实验数据整合...\n")

    # 1. 加载TCA网络(包含反应方程、代谢物碳原子数、反应流量索引)
    metabolite_reaction_dict, flux_balance_reaction_dict, complete_metabolite_dim_dict, flux_name_index_dict = setup_tca_network()

    # 2. 设置输入代谢物集合
    input_metabolite_name_set = {'AcCoA', 'Asp'}

    # 设置目标代谢物产物名称列表
    # Define target metabolites
    target_metabolite_name_list = ['Glu']

    # 3. 设置反应流量向量（初始系数）
    flux_vector = np.array([100, 100, 50, 50, 50, 125, 75, 50], dtype=float)

    # 4. 加载HCT116细胞系实验数据
    experimental_mid_data_obj_dict, model_metabolite_to_standard_name_dict = load_hct116_experimental_data(
        experiment_name="HCT116_WQ2101",
        condition="ctrl",
        index="average" # 使用平均值而不是单个样本
    )


    # 这里有问题

    # 确保将所有实验数据中的代谢物名称添加到映射字典中
    # for exp_metabolite in experimental_mid_data_obj_dict.keys():
    #     if exp_metabolite not in model_metabolite_to_standard_name_dict:
    #         # 尝试通过简单比较找到匹配
    #         for model_metabolite in complete_metabolite_dim_dict.keys():
    #             if model_metabolite.lower() in exp_metabolite.lower() or exp_metabolite.lower() in model_metabolite.lower():
    #                 model_metabolite_to_standard_name_dict[exp_metabolite] = model_metabolite
    #                 print(f"自动添加映射: {exp_metabolite} → {model_metabolite}")
    #                 break


    # 5. 构建模型代谢物到标准名称的映射
    # 如果model_metabolite_to_standard_name_dict中没有所需的映射，需要手动添加
    user_metabolite_to_standard_name_dict = {
        # 实验数据名称: 模型中的名称
        # 输入代谢物
        'd-aspartate': 'Asp',  # 天冬氨酸
        'acetyl-CoA': 'AcCoA',  # 乙酰CoA (2C)

        # TCA循环代谢物
        'citrate': 'Cit',  # 柠檬酸 (6C)
        '2-oxoglutarate': 'AKG',  # α-酮戊二酸 (5C)
        'succinate': 'Suc',  # 琥珀酸 (4C)
        'fumarate': 'Fum',  # 延胡索酸 (4C)
        'oxaloacetic-acid': 'OAC',  # 草酰乙酸 (4C)

        # 输出代谢物
        'glutamate': 'Glu',  # 谷氨酸 (5C)
        #'carbon-dioxide': 'CO2'  # 二氧化碳 (1C)
    }

    # # 将用户定义的代谢物映射合并到模型映射中
    # print("\n更新代谢物名称映射关系:")
    # added_count = 0
    # existing_count = 0
    #
    # # 遍历用户定义的映射关系
    # for exp_name, model_name in user_metabolite_to_standard_name_dict.items():
    #     if exp_name in model_metabolite_to_standard_name_dict:
    #         # 检查映射是否一致
    #         if model_metabolite_to_standard_name_dict[exp_name] == model_name:
    #             print(f"  已存在映射: {exp_name} → {model_name}")
    #             existing_count += 1
    #         else:
    #             # 更新不同的映射
    #             old_value = model_metabolite_to_standard_name_dict[exp_name]
    #             model_metabolite_to_standard_name_dict[exp_name] = model_name
    #             print(f"  更新映射: {exp_name} → {model_name} (原为: {old_value})")
    #             added_count += 1
    #     else:
    #         # 添加新的映射
    #         model_metabolite_to_standard_name_dict[exp_name] = model_name
    #         print(f"  添加映射: {exp_name} → {model_name}")
    #         added_count += 1
    #
    # print(f"\n映射更新完成: 添加/更新了 {added_count} 个映射, 已有 {existing_count} 个相同映射")
    # print(f"最终映射字典包含 {len(model_metabolite_to_standard_name_dict)} 个条目")
    #
    #
    # # 打印最终映射字典
    # print("\n最终代谢物映射关系:")
    # for exp_name, model_name in model_metabolite_to_standard_name_dict.items():
    #     print(f"  {exp_name} → {model_name}")




    # 6. 设置反应流量约束范围
    # 如果不指定特定反应的反应流量范围，代码会使用默认值(线性规划阶段为0.1-5000.0，SLSQP优化阶段为1.0-1000.0，线性规划是必须比SLSQP要小)
    # 但这可能不符合特定反应的生物学特性，导致模型预测偏离实际情况
    specific_flux_range_dict = {
        'v1': (0.1, 1000),  # OAC + AcCoA → Cit
        'v2': (0.1, 1000),  # Cit → AKG + CO2
        'v3': (0.1, 1000),  # AKG → Glu
        'v4': (0.1, 1000),  # AKG → Suc + CO2
        'v5': (0.1, 1000),  # Suc → Fum
        'v6': (0.1, 1000),  # Fum → OAC
        'v7': (0.1, 1000),  # OAC → Fum
        'v8': (0.1, 1000)   # Asp → OAC
    }

    # 7. 设置输入代谢物标记模式
    # 例如：[1,2-13C2]乙酰辅酶A，表示第1和第2个碳原子被13C标记
    AcCoA_labeled_input_metabolite_dict = {
        "AcCoA": [
            {
                "ratio_list": [0, 1],  # 第2个碳原子标记
                "abundance": 0.25,     # 25%的分子有这种标记
            },
            {
                "ratio_list": [1, 1],  # 两个碳原子都标记
                "abundance": 0.25,     # 25%的分子有这种标记
            },
            {
                "ratio_list": [0, 0],  # 没有标记
                "abundance": 0.5,      # 50%的分子没有标记
            }
        ]
    }

    # 8. 运行集成实验数据的MFA分析
    print("\n运行集成实验数据的MFA分析...")
    results = run_mfa_pipeline_with_experimental_data(metabolite_reaction_dict=metabolite_reaction_dict,
                                                      input_metabolite_name_set=input_metabolite_name_set,
                                                      complete_metabolite_dim_dict=complete_metabolite_dim_dict,
                                                      target_metabolite_name_list=target_metabolite_name_list,
                                                      flux_name_index_dict=flux_name_index_dict,
                                                      flux_vector=flux_vector,

                                                      user_metabolite_to_standard_name_dict=user_metabolite_to_standard_name_dict,
                                                      flux_balance_reaction_dict=flux_balance_reaction_dict,
                                                      experimental_mid_data_obj_dict=experimental_mid_data_obj_dict,    # 来源于加载HCT116细胞系实验数据
                                                      model_metabolite_to_standard_name_dict=user_metabolite_to_standard_name_dict, # 替换为模型自定义代谢物
                                                      specific_flux_range_dict=specific_flux_range_dict,    # 自行设置的反应流量约束范围
                                                      common_mix_ratio_range=(0.05, 0.95),  # 通用混合比例范围
                                                      mix_ratio_multiplier=1.0,    # 混合比例乘数

                                                      input_metabolite_dict=AcCoA_labeled_input_metabolite_dict,    # 输入代谢物字典
                                                      verbose=True) # 是否输出详细信息

    # 9. 分析结果
    print("\n结果分析:")

    # 比较优化前后的反应流量值
    initial_flux = results.get('initial_flux_vector')
    optimized_flux = results.get('optimized_flux_vector')

    print("\n反应流量值比较 (初始线性规划 vs 优化后):")
    for flux_name, idx in flux_name_index_dict.items():
        initial = initial_flux[idx] if initial_flux is not None else flux_vector[idx]
        optimized = optimized_flux[idx] if optimized_flux is not None else initial
        print(f"{flux_name}: {initial:.4f} -> {optimized:.4f}")

    # # 获取混合方程信息 nested_mix_equation_dict
    # nested_mix_equation_dict = results.get('nested_mix_equation_dict', {})
    # print(f"\n代谢物混合方程数量: {len(nested_mix_equation_dict)}")
    # for metabolite, equation in list(nested_mix_equation_dict.items())[:3]:  # 仅显示前3个
    #     print(f"  - {metabolite}: {equation}")
    #
    # # 获取预测的MID数据和目标EMU
    # target_emu_name_list = results.get('target_emu_name_list', [])
    # target_indices = results.get('target_indices', [])
    # predicted_mid_data = results.get('predicted_mid_data', [])
    #
    # # 打印部分目标代谢物的预测MID
    # print("\n预测的代谢物MID:")
    # for i, idx in enumerate(target_indices):
    #     emu_name = target_emu_name_list[i]
    #     print(f"  - {emu_name}: {predicted_mid_data[idx]}")
    #
    # # 打印预测与实验MID比较
    # comparison_data = results.get('comparison_data', [])
    # if comparison_data:
    #     print("\n优化后预测MID vs 质谱实验MID:")
    #     for item in comparison_data:
    #         print(f"\n{item['model_metabolite']} RMSE: {item['RMSE']:.4f}")
    #         print("  预测\t实验")
    #         for i in range(len(item['predicted_mid'])):
    #             print(f"M+{i}: {item['predicted_mid'][i]:.4f}\t{item['experimental_mid'][i]:.4f}")

    # 10. 返回结果
    # 通过这些修改，我们添加了线性规划初始值生成并使用SLSQP优化算法使预测与实验数据间的loss收敛到最小，同时修复了数据加载问题。
    return results

if __name__ == '__main__':
    test_tca_cycle_mfa()
