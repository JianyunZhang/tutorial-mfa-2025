from scripts.src.common.built_in_packages import ValueEnum

# scripts/main.py 是一个参数设置模块，主要用于设置计算相关的命令行参数
# 该模块作为计算功能的参数配置中心，将不同的计算功能模块组织在 computation 命令下，便于统一管理和调用

# ComputationFunction 枚举类
# ValueEnum 是一个自定义枚举类，用于更方便的字符串处理
class ComputationFunction(ValueEnum):
    simulation = 'simulation'       # 用于生成模拟的代谢物同位素分布(MID)数据
    sensitivity = 'sensitivity'     # 用于进行代谢流分析模型的敏感性分析
    experiments = 'experiments'     # 用于处理实验数据分析
    standard_name = 'standard_name' # 用于输出标准化的代谢物名称

# arg_setting 函数: 这个函数设置计算相关的子命令解析器
def arg_setting(subparsers):
    # 当没有指定子命令时显示帮助信息的回调函数
    def print_computation_help(args):
        computation_parser.print_help()

    # 创建主计算命令解析器
    computation_parser = subparsers.add_parser('computation', help='Run computation functions')

    # 创建计算子命令解析器
    computation_subparsers = computation_parser.add_subparsers(
        title='Commands',
        description='Different content for computation',
        help='Decide to run different analysis functions')

    # 导入并设置四个主要功能模块的参数

    # 导入并配置模拟数据生成模块的参数
    # 模拟数据生成 (simulation)
    from scripts.src.simulated_data.simulated_data_generator_main import arg_setting as simulation_arg_setting
    simulation_arg_setting(computation_subparsers, ComputationFunction.simulation)

    # 导入并配置实验数据分析模块的参数
    # 实验数据分析 (experiments)
    from scripts.src.experimental_data_analysis.experimental_data_analysis_main import arg_setting as \
        experiments_arg_setting
    experiments_arg_setting(computation_subparsers, ComputationFunction.experiments)

    # 导入并配置模型敏感性分析模块的参数
    # 模型数据敏感性分析 (sensitivity)
    from scripts.src.model_data_sensitivity.model_data_sensitivity_main import arg_setting as \
        model_data_sensitivity_arg_setting
    model_data_sensitivity_arg_setting(computation_subparsers, ComputationFunction.sensitivity)

    # 导入并配置标准名称输出模块的参数
    # 标准名称输出 (standard_name)
    from scripts.src.standard_name_output.standard_name_output_main import arg_setting as \
        standard_name_output_arg_setting
    standard_name_output_arg_setting(computation_subparsers, ComputationFunction.standard_name)

    # 设置默认行为：显示帮助信息
    computation_parser.set_defaults(func=print_computation_help)
