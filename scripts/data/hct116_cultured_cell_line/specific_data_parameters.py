from .data_metabolite_to_standard_name_dict import data_metabolite_to_standard_name_dict
from scripts.src.common.config import Direct, DataType, Keywords as CommonKeywords
from scripts.src.core.common.classes import TransformDict
from ..complete_dataset_class import CompleteDataset
from ..common_functions import average_mid_data_dict

# scripts/data/hct116_cultured_cell_line/specific_data_parameters.py
# 定义了一个 SpecificParameters 类，该类继承自 CompleteDataset。
# 这个类的主要功能是处理特定实验的数据参数。以下是该文件的主要功能：
# 包含一些关键字和映射字典，用于实验条件和索引的转换
# 这个类的主要目的是管理和处理特定实验的数据参数，并提供方法来生成项目名称和返回数据集。
class Keyword(object):
    experiments = 'experiments'
    condition = 'condition'
    index = 'index'

    ctrl = 'ctrl'
    carcinoma = 'carcinoma'
    brain = 'brain'

    prefix = 'prefix'

    condition_mapping_dict = TransformDict(**{
        'Ctrl': 'ctrl',
        'DMSO': 'ctrl',
        'NCT25um': '25um'
    })
    index_average_list = [1, 2, 3]


class SpecificParameters(CompleteDataset):
    # 初始化方法 __init__：设置一些实验相关的参数和文件路径。
    def __init__(self):
        super(SpecificParameters, self).__init__()
        self.mixed_compartment_list = ('c', 'm')
        self.current_direct = "{}/hct116_cultured_cell_line".format(Direct.data_direct)
        self.file_path = "{}/13C-Glucose_tracing_Mike.xlsx".format(self.current_direct)
        self.sheet_name_dict = {
            "HCT116_WQ2101_PHGDH": {
                Keyword.experiments: 'HCT116_WQ2101',
                Keyword.prefix: 'HCT116'
            },
            "BT20_WQ2101_PHGDH": {
                Keyword.experiments: 'BT20_WQ2101',
                Keyword.prefix: 'BT20'
            },
            "HCT116_NCT503_PHGDH": {
                Keyword.experiments: 'HCT116_NCT503',
                Keyword.prefix: 'II_13C'
            }
        }
        self.test_experiment_name_prefix = "HCT116_WQ2101_PHGDH"
        self.test_col = 'HCT116_Ctrl_1'
        self._complete_data_parameter_dict_dict = {
            current_sheet_name: {
                'xlsx_file_path': self.file_path,
                'xlsx_sheet_name': current_sheet_name,
                'index_col_name': CommonKeywords.metabolite_name_col,
                'mixed_compartment_list': self.mixed_compartment_list,
                'to_standard_name_dict': data_metabolite_to_standard_name_dict}
            for current_sheet_name in self.sheet_name_dict.keys()}
        self._test_data_parameter_dict_dict = {
            DataType.test: {
                'xlsx_file_path': self.file_path,
                'xlsx_sheet_name': self.test_experiment_name_prefix,
                'index_col_name': CommonKeywords.metabolite_name_col,
                'mixed_compartment_list': self.mixed_compartment_list,
                'to_standard_name_dict': data_metabolite_to_standard_name_dict}}

    # 静态方法 project_name_generator：生成项目名称。
    @staticmethod
    def project_name_generator(experiment_name, condition_name, index_num):
        return '{}__{}__{}'.format(experiment_name, condition_name, index_num)

    # 方法 add_data_sheet：将数据表添加到完整的数据集中。
    def add_data_sheet(self, sheet_name, current_data_dict):
        final_result_dict = self.complete_dataset
        if sheet_name == DataType.test:
            final_result_dict[DataType.test] = current_data_dict
        else:
            current_information_dict = self.sheet_name_dict[sheet_name]
            experiment_name = current_information_dict[Keyword.experiments]
            prefix = current_information_dict[Keyword.prefix]
            if experiment_name not in final_result_dict:
                final_result_dict[experiment_name] = {}
            for data_label, specific_data_dict in current_data_dict.items():
                raw_condition_name, index_str = data_label[len(prefix) + 1:].split('_')
                condition_name = Keyword.condition_mapping_dict[raw_condition_name]
                if condition_name not in final_result_dict[experiment_name]:
                    final_result_dict[experiment_name][condition_name] = {}
                final_result_dict[experiment_name][condition_name][int(index_str)] = specific_data_dict
    # 方法 _complete_return_dataset：根据参数字典返回完整的数据集。
    def _complete_return_dataset(self, param_dict):
        experiment_name = param_dict[Keyword.experiments]
        condition_name = param_dict[Keyword.condition]
        index_num = param_dict[Keyword.index]
        if index_num == CommonKeywords.average:
            final_target_metabolite_data_dict = average_mid_data_dict(
                self.complete_dataset[experiment_name][condition_name], Keyword.index_average_list)
        else:
            final_target_metabolite_data_dict = self.complete_dataset[
                experiment_name][condition_name][index_num]
        final_input_metabolite_data_dict = None
        project_name = self.project_name_generator(experiment_name, condition_name, index_num)
        return project_name, final_target_metabolite_data_dict, final_input_metabolite_data_dict
    # 方法 _test_return_dataset：返回测试数据集
    def _test_return_dataset(self, param_dict=None):
        final_target_metabolite_data_dict = self.complete_dataset[
            DataType.test][self.test_col]
        final_input_metabolite_data_dict = None
        project_name = DataType.test
        return project_name, final_target_metabolite_data_dict, final_input_metabolite_data_dict
