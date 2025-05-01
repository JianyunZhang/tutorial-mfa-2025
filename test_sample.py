import numpy as np

from scripts.src.core.common.config import CoreConstants
from scripts.src.core.model.model_class import Reaction, Node
from common_test_pipeline import run_mfa_pipeline


def setup_sample_network():
    """Setup Sample metabolic network"""
    metabolite_reaction_dict = {
        # 输入代谢物
        'A': [
        ],
        'B': [Reaction(id = 'V1',
                       sub = [Node(name = 'A', coefficient = 1, carbon_composition_string = ['a', 'b', 'c'])],
                       pro = [Node(name = 'B', coefficient = 1, carbon_composition_string = ['a', 'b', 'c'])],
                       reverse = False),
              Reaction(id = 'V3',
                       sub = [Node(name = 'D', coefficient = 1, carbon_composition_string = ['a', 'b', 'c'])],
                       pro = [Node(name = 'B', coefficient = 1, carbon_composition_string = ['a', 'b', 'c'])],
                       reverse = False),
        ],
        'C': [Reaction(id = 'V4',
                       sub = [Node(name = 'B', coefficient = 1, carbon_composition_string = ['a', 'b', 'c'])],
                       pro = [Node(name = 'C', coefficient = 1, carbon_composition_string = ['b', 'c']),
                            Node(name = 'E', coefficient = 1, carbon_composition_string = ['a'])],
                       reverse = False)
        ],
        'D': [Reaction(id = 'V2',
                       sub = [Node(name = 'B', coefficient = 1, carbon_composition_string = ['a', 'b', 'c'])],
                       pro = [Node(name = 'D', coefficient = 1, carbon_composition_string = ['a', 'b', 'c'])],
                       reverse = False),
              Reaction(id = 'V5',
                       sub = [Node(name = 'B', coefficient = 1, carbon_composition_string = ['a', 'b', 'c']),
                            Node(name = 'C', coefficient = 1, carbon_composition_string = ['d', 'e'])],
                       pro = [Node(name = 'D', coefficient = 1, carbon_composition_string = ['b', 'c', 'd']),
                            Node(name = 'E', coefficient = 1, carbon_composition_string = ['a']),
                            Node(name = 'E', coefficient = 1, carbon_composition_string = ['e'])],
                       reverse = False)
        ],
        'E': [Reaction(id = 'V4',
                       sub = [Node(name = 'B', coefficient = 1, carbon_composition_string = ['a', 'b', 'c'])],
                       pro = [Node(name = 'C', coefficient = 1, carbon_composition_string = ['b', 'c']),
                            Node(name = 'E', coefficient = 1, carbon_composition_string = ['a'])],
                       reverse = False),
              Reaction(id = 'V5',
                       sub = [Node(name = 'B', coefficient = 1, carbon_composition_string = ['a', 'b', 'c']),
                            Node(name ='C', coefficient = 1, carbon_composition_string = ['d', 'e'])],
                       pro = [Node(name = 'D', coefficient = 1, carbon_composition_string = ['b', 'c', 'd']),
                            Node(name = 'E', coefficient = 1, carbon_composition_string = ['a']),
                            Node(name = 'E', coefficient = 1, carbon_composition_string = ['e'])],
                       reverse = False)
        ],
        'F':[Reaction(id = 'V6',
                       sub = [Node(name = 'D', coefficient = 1, carbon_composition_string=['a', 'b', 'c'])],
                       pro = [Node(name = 'F', coefficient = 1, carbon_composition_string=['a', 'b', 'c'])],
                       reverse = False)
        ]
    }

    # 定义代谢物碳原子数
    complete_metabolite_dim_dict = {
        'A': 3,  # A有3个碳原子
        'B': 3,  # B有3个碳原子
        'C': 2,  # C有2个碳原子
        'D': 3,  # D有3个碳原子
        'E': 1,  # E有1个碳原子
        'F': 3   # F有3个碳原子
    }

    # Define flux indices
    flux_name_index_dict = {
        'V1': 0,
        'V2': 1,
        'V3': 2,
        'V4': 3,
        'V5': 4,
        'V6': 5
    }

    return metabolite_reaction_dict, complete_metabolite_dim_dict, flux_name_index_dict


def test_sample_mfa():
    """Sample metabolic flux analysis test"""
    print("Starting Simple network MFA analysis...\n")

    # Setup Sample network
    metabolite_reaction_dict, complete_metabolite_dim_dict, flux_name_index_dict = setup_sample_network()

    # Define input metabolites
    input_metabolite_name_set = {'A'}

    # Define target metabolites
    target_metabolite_names = ['F']

    # Define flux values
    flux_vector = np.array([100, 110, 50, 20, 20, 80])

    # Define custom labeling
    # 自定义碳原子标记，只考虑A第二个碳原子标记的情况，其他输入代谢物均考虑为自然丰度标记
    custom_labeling = {
        # 为A的所有可能EMU组合提供数据
        # 每个数组的元素和为1（表示概率分布）
        # 数组长度必须与EMU中的碳原子数量匹配
        # 键名必须与EMU方程中使用的标识符一致
        # A的单碳EMU：[M+0, M+1]
        # 'A__100': np.array([1, 0]),  # A的第1个碳原子
        # 'A__010': np.array([0, 1]),  # A的第2个碳原子
        # 'A__001': np.array([1, 0]),  # A的第3个碳原子
        'A__010': np.convolve(np.convolve(np.convolve(np.array([1]),np.array([1, 0])),np.array([0, 1])),np.array([1,0])),
        'A__001': np.convolve(np.convolve(np.convolve(np.array([1]),np.array([1, 0])),np.array([1, 0])),np.array([1,0]))

        # A的双碳EMU：[M+0, M+1, M+2]
        # 'A__011': np.array([0, 1, 0]),  # A的第2,3个碳原子
        # 'A__101': np.array([1, 0, 0]),  # A的第1,3个碳原子
        # 'A__110': np.array([0, 1, 0]),  # A的第1,2个碳原子

        # A的三碳EMU：[M+0, M+1, M+2, M+3]
        # 'A__111': np.array([0, 1, 0, 0])  # A的第1,2,3个碳原子的同位素分布
    }

    A_labeled_ratio_list = [0, 1, 0]
    A_labeled_input_metabolite_dict = {
        "A": [
            {
                "ratio_list": A_labeled_ratio_list,
                "abundance": 1,
            }
        ]
    }

    # Run the pipeline
    results = run_mfa_pipeline(
        metabolite_reaction_dict,
        input_metabolite_name_set,
        complete_metabolite_dim_dict,
        target_metabolite_names,
        flux_name_index_dict,
        flux_vector,
        # custom_labeling=custom_labeling,
        input_metabolite_dict=A_labeled_input_metabolite_dict,
        verbose=True
    )

    # Additional validation
    print("\nResult validation:")
    target_indices = results['target_indices']
    target_emu_name_list = results['target_emu_name_list']
    complete_predicted_mid_data_list = results['predicted_mid_data']

    for target_name, target_idx in zip(target_emu_name_list, target_indices):
        if target_idx < len(complete_predicted_mid_data_list):
            mid_sum = sum(complete_predicted_mid_data_list[target_idx])
            print(f"{target_name} MID sum: {mid_sum:.6f} (should be close to 1)")


if __name__ == '__main__':
    test_sample_mfa()