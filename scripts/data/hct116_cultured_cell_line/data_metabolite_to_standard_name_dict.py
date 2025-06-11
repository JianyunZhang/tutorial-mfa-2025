from scripts.src.core.common.classes import TransformDict
#  defines a dictionary that maps various metabolite names to their standard names.
#  This dictionary is implemented using the TransformDict class from the scripts/src/core/common/classes module.
#  The primary purpose of this file is to provide a standardized naming convention for metabolites, which can be used throughout the project to ensure consistency in data processing and analysis.
#  The data_metabolite_to_standard_name_dict dictionary is used to map various metabolite names to their standard names.
#  This dictionary is implemented using the TransformDict class from the scripts/src/core/common/classes module.
#  The primary purpose of this file is to provide a standardized naming convention for metabolites, which can be used throughout the project to ensure consistency in data processing and analysis.
#  Importing the TransformDict class:
#  The TransformDict class is imported from the scripts/src/core/common/classes module.
#  Defining the data_metabolite_to_standard_name_dict dictionary:
#  该字典包含键-值对，其中键是代谢物的各种名称，值是相应的标准名称。
#  This dictionary contains key-value pairs where the keys are various names of metabolites, and the values are their corresponding standard names.
#  The dictionary includes mappings for several metabolites, such as lactate, malate, glucose, and various amino acids, ensuring that different naming conventions are unified under a single standard name.
#  This is useful for data consistency and avoiding confusion when dealing with metabolite data in different parts of the project.
data_metabolite_to_standard_name_dict = TransformDict({
    '(r,s)-lactate': 'lactate',
    '(s)-malate': 'malate',
    '2-oxoglutarate': 'a-ketoglutarate',

    'd-glucose': 'glucose',
    'g6p': 'glucose 6-phosphate',
    'f6p': 'fructose 6-phosphate',
    'fructose-1,6-bisphosphate': 'fructose 1,6-bisphosphate',
    '3pg': '3-phosphoglycerate',
    '2pg': '2-phosphoglycerate',

    'd-ribose-5-phosphate': 'ribose 5-phosphate',
    'd-erythrose 4-phosphate': 'erythrose 4-phosphate',

    'l-lysine': 'lysine',
    'l-tryptophan': 'tryptophan',
    'l-isoleucine': 'isoleucine',
    'l-leucine': 'leucine',
    'l-glutamate': 'glutamate',
    'd-aspartate': 'aspartate',
})
