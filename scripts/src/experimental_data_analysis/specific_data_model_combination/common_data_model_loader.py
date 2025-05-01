from ..inventory import DataModelType

# 该文件的主要功能
# 提供统一的数据模型加载接口
# 支持多种实验数据模型类型：
# 肾癌体内灌注模型（多个变体）
# 肺肿瘤体内灌注模型
# 结肠癌细胞系模型（多个变体）
# HCT116细胞系模型（多个变体）
# 根据传入的模型枚举类型动态导入对应的模型模块
# 这个加载器的设计使得系统可以灵活地处理不同类型的实验数据模型，同时保持了代码的模块化和可维护性。

def common_data_model_function_loader(model_name):
    """
        根据模型名称动态加载对应的数据模型模块
        Args:
            model_name: DataModelType枚举类型，指定要加载的模型名称
        Returns:
            data_model_object: 加载的数据模型模块对象
        Raises:
            ValueError: 当提供的model_name不在支持的类型中时抛出
    """
    # 肾癌体内灌注相关模型
    if model_name == DataModelType.renal_carcinoma_invivo_infusion:
        from . import renal_carcinoma_invivo_infusion as data_model_object
    elif model_name == DataModelType.renal_carcinoma_invivo_infusion_squared_loss:
        from . import renal_carcinoma_invivo_infusion_squared_loss as data_model_object
    elif model_name == DataModelType.renal_carcinoma_invivo_infusion_traditional_method:
        from . import renal_carcinoma_invivo_infusion_traditional_method as data_model_object
    elif model_name == DataModelType.renal_carcinoma_invivo_infusion_with_glns_m:
        from . import renal_carcinoma_invivo_infusion_with_glns_m as data_model_object
    elif model_name == DataModelType.renal_carcinoma_invivo_infusion_with_glns_m_traditional_method:
        from . import renal_carcinoma_invivo_infusion_with_glns_m_traditional_method as data_model_object

    # 肺肿瘤体内灌注模型
    elif model_name == DataModelType.lung_tumor_invivo_infusion:
        from . import lung_tumor_invivo_infusion as data_model_object

    # 结肠癌细胞系相关模型
    elif model_name == DataModelType.colon_cancer_cell_line:
        from . import colon_cancer_cell_line as data_model_object
    elif model_name == DataModelType.colon_cancer_cell_line_squared_loss:
        from . import colon_cancer_cell_line_squared_loss as data_model_object
    elif model_name == DataModelType.colon_cancer_cell_line_traditional_method:
        from . import colon_cancer_cell_line_traditional_method as data_model_object
    elif model_name == DataModelType.colon_cancer_cell_line_with_glns_m:
        from . import colon_cancer_cell_line_with_glns_m as data_model_object
    elif model_name == DataModelType.colon_cancer_cell_line_with_glns_m_traditional_method:
        from . import colon_cancer_cell_line_with_glns_m_traditional_method as data_model_object

    # HCT116细胞系相关模型
    elif model_name == DataModelType.hct116_cultured_cell_line:
        from . import hct116_cultured_cell_line as data_model_object
    elif model_name == DataModelType.hct116_cultured_cell_line_with_glns_m:
        from . import hct116_cultured_cell_line_with_glns_m as data_model_object
    else:
        raise ValueError()

    return data_model_object
