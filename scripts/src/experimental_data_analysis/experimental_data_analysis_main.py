from .inventory import DataModelType, RunningMode, data_model_comment
from ..common.built_in_packages import argparse, mp

# 这个模块是代谢流分析(MFA)项目中实验数据分析部分的入口点
# 配置实验数据分析的命令行参数系统
# 支持以下参数设置：
# 测试模式 (-t/--test_mode)
# 并行进程数 (-p/--parallel_num)
# 运行模式 (running_mode)
# 数据模型 (data_model)

def arg_setting(computation_subparsers, experiments_enum):
    """
        设置实验数据分析的命令行参数
        Args:
            computation_subparsers: 计算模块的子命令解析器
            experiments_enum: 实验分析的枚举值
    """

    # 定义实验分析命令的处理函数
    def experiments(args):
        main(experimental_data_analysis_parser, args)

    # 创建实验数据分析的命令行解析器
    experimental_data_analysis_parser = computation_subparsers.add_parser(
        experiments_enum.value, help='Run MFA for several experimental data analyses',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog='Definition of data_model_name:\n\n{}'.format(
            '\n'.join([
                f'{enum_item.name:<50}{data_model_comment[enum_item]}'
                for enum_item in DataModelType
            ])
        )
    )

    # 添加测试模式参数
    experimental_data_analysis_parser.add_argument(
        '-t', '--test_mode', action='store_true', default=False,
        help='Whether the code is executed in test mode, which means less sample number and shorter time.'
    )

    # 添加并行进程数参数
    experimental_data_analysis_parser.add_argument(
        '-p', '--parallel_num', type=int, default=None,
        help='Number of parallel processes. If not provided, it will be selected according to CPU cores.'
    )

    # 设置运行模式参数
    running_mode_display = '{}'.format(',  '.join([running_mode.value for running_mode in RunningMode]))
    experimental_data_analysis_parser.add_argument(
        'running_mode', nargs='?', type=RunningMode, choices=list(RunningMode),
        help='Running mode of experimental data analysis', default=None, metavar=running_mode_display)

    # 设置数据模型参数
    experimental_data_analysis_parser.add_argument(
        'data_model', nargs='?', type=DataModelType, choices=list(DataModelType),
        help='The data-model combination that need to calculate. Detailed list is attached below',
        default=None, metavar='data_model_name')

    # 设置默认处理函数
    experimental_data_analysis_parser.set_defaults(func=experiments)


def main(experimental_data_analysis_parser, args):
    """
        主函数：处理实验数据分析的入口
        Args:
            experimental_data_analysis_parser: 实验数据分析的命令行解析器
            args: 解析的命令行参数
    """
    running_mode = args.running_mode
    if running_mode is None:
        # 如果没有指定运行模式，显示帮助信息
        experimental_data_analysis_parser.print_help()
    else:
        # 导入并执行数据分析主函数
        from .common_functions import data_analysis_main
        mp.set_start_method('spawn')    # 设置多进程启动方法
        data_analysis_main(running_mode, args.data_model, args.test_mode, args.parallel_num)
