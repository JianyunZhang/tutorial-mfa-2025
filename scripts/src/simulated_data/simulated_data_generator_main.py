
# 该函数的主要功能：
# 接收主解析器和模拟枚举值作为参数
# 创建模拟数据生成的命令行解析器
# 配置多个命令行参数选项:
# -f/--new_flux: 是否生成新的通量优化
# -n/--with_noise: 是否添加噪声
# --with_glns_m: 是否添加 GLNS_m 通量
# -i/--index: 区分不同模拟数据的索引
# -b/--batch_num: 生成模拟解决方案的批次数量
# 这些参数可以通过命令行方式灵活配置模拟数据的生成过程。
# 例如 python main.py computation simulation -f -n -i 1 -b 5

def arg_setting(computation_subparsers, simulation_enum):
    # 定义内部函数，用于处理模拟命令的实际执行
    def simulation(args):
        # 调用主函数处理模拟数据生成
        main(simulation_parser, args)

    # 创建模拟数据生成的子命令解析器
    simulation_parser = computation_subparsers.add_parser(
        simulation_enum.value,  # 使用枚举值作为命令名
        help='Generate simulated MID data'  # 命令的帮助说明
    )

    # 添加新通量优化参数
    simulation_parser.add_argument(
        '-f', '--new_flux', # 参数的短格式和长格式
        action='store_true', # 标志参数，存在则为True
        default=False, # 默认值为False
        help='Generate new flux optimized from PHDGH mass spectrometry data' # 参数说明
    )

    # 添加噪声参数
    simulation_parser.add_argument(
        '-n', '--with_noise',
        action='store_true',
        default=False,
        help='Add noise to generated MID data to mimic real case.'
    )

    # 添加GLNS_m通量参数
    simulation_parser.add_argument(
        '--with_glns_m', action='store_true', default=False,
        help='Add GLNS_m flux to final model.'
    )

    # 添加索引参数
    simulation_parser.add_argument(
        '-i', '--index',
        type=int,   # 参数类型为整数
        default=None,   # 默认值为None
        help='Index to distinguish different simulated data.'
    )

    # 添加批次数量参数
    simulation_parser.add_argument(
        '-b', '--batch_num',
        type=int,
        default=1,  # 默认生成1个批次
        help='Number to generate batched simulated solutions.'
    )

    # 设置默认的处理函数
    simulation_parser.set_defaults(func=simulation)


def main(simulation_parser, args):
    from .common_functions import simulated_mid_data_generator
    simulated_mid_data_generator(args.new_flux, args.batch_num, args.index, args.with_noise, args.with_glns_m)
