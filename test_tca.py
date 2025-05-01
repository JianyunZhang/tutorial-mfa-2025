import numpy as np

from scripts.src.core.common.config import CoreConstants
from scripts.src.core.model.model_class import Reaction, Node
from common_test_pipeline import run_mfa_pipeline


def setup_tca_network():
    """Setup TCA cycle metabolic network"""
    metabolite_reaction_dict = {
        # 输入代谢物
        'Asp': [],
        'AcCoA': [],

        # TCA循环代谢物
        'Cit': [
            # 反应1: OAC + AcCoA → Cit (碳原子转换: abcd + ef → dcbfea)
            Reaction(id='v1',
                     sub=[Node(name='OAC', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd']),
                          Node(name='AcCoA', coefficient=1, carbon_composition_string=['e', 'f'])],
                     pro=[Node(name='Cit', coefficient=1, carbon_composition_string=['d', 'c', 'b', 'f', 'e', 'a'])],
                     reverse=False)
        ],
        'AKG': [
            # 反应2: Cit → AKG + CO2 (碳原子转换: abcdef → abcde + f)
            Reaction(id='v2',
                     sub=[Node(name='Cit', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd', 'e', 'f'])],
                     pro=[Node(name='AKG', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd', 'e']),
                          Node(name='CO2', coefficient=1, carbon_composition_string=['f'])],
                     reverse=False)
        ],
        'Glu': [
            # 反应3: AKG → Glu (碳原子转换: abcde → abcde)
            Reaction(id='v3',
                     sub=[Node(name='AKG', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd', 'e'])],
                     pro=[Node(name='Glu', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd', 'e'])],
                     reverse=False)
        ],
        'Suc': [
            # 反应4: AKG → Suc + CO2 (碳原子转换: abcde → bcde + a)
            Reaction(id='v4',
                     sub=[Node(name='AKG', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd', 'e'])],
                     pro=[Node(name='Suc', coefficient=1, carbon_composition_string=['b', 'c', 'd', 'e']),
                          Node(name='CO2', coefficient=1, carbon_composition_string=['a'])],
                     reverse=False)
        ],
        'Fum': [
            # 反应5: Suc → Fum (考虑碳原子对称性)
            # Suc(abcd) -> Fum(abcd) 和 Suc(abcd) -> Fum(dcba) 各占 50%
            Reaction(id='v5',
                     sub=[Node(name='Suc', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd'])],
                     pro=[Node(name='Fum', coefficient=0.5, carbon_composition_string=['a', 'b', 'c', 'd']),
                          Node(name='Fum', coefficient=0.5, carbon_composition_string=['d', 'c', 'b', 'a'])],
                     reverse=False),
            # 反应7: OAC → Fum (考虑碳原子对称性)
            # OAC(abcd) -> Fum(abcd) 和 OAC(abcd) -> Fum(dcba) 各占 50%
            Reaction(id='v7',
                     sub=[Node(name='OAC', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd'])],
                     pro=[Node(name='Fum', coefficient=0.5, carbon_composition_string=['a', 'b', 'c', 'd']),
                          Node(name='Fum', coefficient=0.5, carbon_composition_string=['d', 'c', 'b', 'a'])],
                     reverse=False)
        ],
        'OAC': [
            # 反应6: Fum → OAC (考虑碳原子对称性)
            # Fum(abcd) -> OAC(abcd) 和 Fum(dcba) -> OAC(abcd) 各占 50%
            # 需要定义两个Fum输入节点，代表其对称性
            Reaction(id='v6',
                     sub=[Node(name='Fum', coefficient=0.5, carbon_composition_string=['a', 'b', 'c', 'd']),
                          Node(name='Fum', coefficient=0.5, carbon_composition_string=['d', 'c', 'b', 'a'])],
                     pro=[Node(name='OAC', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd'])],
                     reverse=False),
            # 反应8: Asp → OAC (abcd → abcd)
            Reaction(id='v8',
                     sub=[Node(name='Asp', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd'])],
                     pro=[Node(name='OAC', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd'])],
                     reverse=False)
        ],

        # 输出代谢物
        'CO2': [
            # 反应2: Cit → AKG + CO2 (碳原子转换: abcdef → abcde + f)
            Reaction(id='v2',
                     sub=[Node(name='Cit', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd', 'e', 'f'])],
                     pro=[Node(name='AKG', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd', 'e']),
                          Node(name='CO2', coefficient=1, carbon_composition_string=['f'])],
                     reverse=False),
            Reaction(id='v4',
                     sub=[Node(name='AKG', coefficient=1, carbon_composition_string=['a', 'b', 'c', 'd', 'e'])],
                     pro=[Node(name='Suc', coefficient=1, carbon_composition_string=['b', 'c', 'd', 'e']),
                          Node(name='CO2', coefficient=1, carbon_composition_string=['a'])],
                     reverse=False)
        ]  # 二氧化碳 (1C)
    }

    # 定义代谢物碳原子数
    complete_metabolite_dim_dict = {
        'AcCoA': 2,     # 乙酰CoA (2C)
        'OAC': 4,       # 草酰乙酸 (4C)
        'Asp': 4,       # 天冬氨酸 (4C)
        'Cit': 6,       # 柠檬酸 (6C)
        'AKG': 5,       # α-酮戊二酸 (5C)
        'Glu': 5,       # 谷氨酸 (5C)
        'Suc': 4,       # 琥珀酸 (4C)
        'Fum': 4,       # 延胡索酸 (4C)
        'CO2': 1        # 二氧化碳 (1C)
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

    return metabolite_reaction_dict, complete_metabolite_dim_dict, flux_name_index_dict


def test_tca_cycle_mfa():
    """TCA cycle metabolic flux analysis test"""
    print("Starting TCA cycle MFA analysis...\n")

    # Setup TCA network
    metabolite_reaction_dict, complete_metabolite_dim_dict, flux_name_index_dict = setup_tca_network()

    # Define input metabolites
    input_metabolite_name_set = {'AcCoA', 'Asp'}

    # Define target metabolites
    target_metabolite_name_list = ['Glu']

    # Define flux values
    flux_vector = np.array([100, 100, 50, 50, 50, 125, 75, 50])

    # Define custom labeling
    custom_labeling = {
        # AcCoA的双碳EMU: [M+0, M+1, M+2]
        # 'AcCoA__11': np.array([0.5, 0.25, 0.25]),  # 25% 1,2-13C 双碳标记, 剩下 75% 无 13C 标记
        'AcCoA__11': 0.5 * np.convolve(np.convolve(np.array([1]), np.array([1, 0])),np.array([1, 0]))
                     + 0.25 * np.convolve(np.convolve(np.array([1]),np.array([1, 0])),np.array([0, 1]))
                     + 0.25 * np.convolve(np.convolve(np.array([1]),np.array([0, 1])),np.array([0, 1])),
        # AcCoA的单碳EMU: [M+0, M+1]
        # 'AcCoA__10': np.array([0.75, 0.25]),      # 碳1: 100% 无标记
        'AcCoA__10': 0.5 * np.convolve(np.convolve(np.array([1]), np.array([1, 0])),np.array([1, 0]))
                     + 0.25 * np.convolve(np.convolve(np.array([1]), np.array([1, 0])), np.array([1, 0]))
                     + 0.25 * np.convolve(np.convolve(np.array([1]), np.array([0, 1])), np.array([1, 0])),
        # 'AcCoA__01': np.array([0.5, 0.5])    # 25% 2-13C 单碳标记, 剩下 75% 无 13C 标记
        'AcCoA__01': 0.5 * np.convolve(np.convolve(np.array([1]), np.array([1, 0])),np.array([1, 0]))
                     + 0.25 * np.convolve(np.convolve(np.array([1]), np.array([1, 0])), np.array([0, 1]))
                     + 0.25 * np.convolve(np.convolve(np.array([1]), np.array([1, 0])), np.array([0, 1]))
    }


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

    # Run the pipeline with updated parameters
    results = run_mfa_pipeline(
        metabolite_reaction_dict,
        input_metabolite_name_set,
        complete_metabolite_dim_dict,
        target_metabolite_name_list,
        flux_name_index_dict,
        flux_vector,
        #custom_labeling=custom_labeling,
        input_metabolite_dict=AcCoA_labeled_input_metabolite_dict,
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
    test_tca_cycle_mfa()
