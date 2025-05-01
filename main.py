# 导入 argparse 用于处理命令行参数
import argparse
# 从 scripts.main 导入 arg_setting 函数用于设置计算相关的参数
from scripts.main import arg_setting as computation_arg_setting
# from figures.main import arg_setting as figure_arg_setting

# 这是一个使用 argparse 模块实现命令行参数解析的主程序文件
def main():
    # 创建 ArgumentParser 对象，设置程序名和描述
    parser = argparse.ArgumentParser(
        prog = 'MFA_development',
        description = 'Code for development of MFA by Shiyu Liu.')
    # 创建子命令解析器，允许程序支持多个子命令
    # 增加属性：在命令行中，该py文件希望用户能给他一个参数，最终将之转化为：args.square
    subparsers = parser.add_subparsers(
        title = 'Commands',
        description = 'Different command for this code',
        help = 'Decide to run analysis or to generate figures based on analysis')
    # 调用 computation_arg_setting 设置计算相关的参数
    computation_arg_setting(subparsers)
    # 注释掉的 figure_arg_setting 可能用于图形生成相关的参数设置
    # figure_arg_setting(subparsers)

    # 参数处理逻辑
    # 解析命令行参数到 args 对象
    # 属性给与args实例：add_argument 返回到 args 子类实例
    args = parser.parse_args()
    # 使用 try-except 结构处理命令行参数：
    try:
        current_func = args.func
    # 如果没有提供有效的子命令（args.func 不存在）
    except AttributeError:
        # 显示帮助信息
        parser.print_help()
    # 如果提供了有效的子命令
    else:
        # 执行对应的功能函数
        args.func(args)

# 通过 if __name__ == '__main__' 确保代码只在直接运行时执行
if __name__ == '__main__':
    main()
