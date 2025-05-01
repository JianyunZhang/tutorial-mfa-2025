from ..common.built_in_packages import ValueEnum

# 数据模型类型枚举，定义了不同类型的实验数据模型
class DataModelType(ValueEnum):
    # 结肠癌细胞系相关模型
    colon_cancer_cell_line = 'colon_cancer_cell_line'
    colon_cancer_cell_line_squared_loss = 'colon_cancer_cell_line_squared_loss'
    colon_cancer_cell_line_traditional_method = 'colon_cancer_cell_line_traditional_method'
    colon_cancer_cell_line_with_glns_m = 'colon_cancer_cell_line_with_glns_m'
    colon_cancer_cell_line_with_glns_m_traditional_method = 'colon_cancer_cell_line_with_glns_m_traditional_method'
    # 肾癌体内灌注相关模型
    renal_carcinoma_invivo_infusion = 'renal_carcinoma_invivo_infusion'
    renal_carcinoma_invivo_infusion_squared_loss = 'renal_carcinoma_invivo_infusion_squared_loss'
    renal_carcinoma_invivo_infusion_traditional_method = 'renal_carcinoma_invivo_infusion_traditional_method'
    renal_carcinoma_invivo_infusion_with_glns_m = 'renal_carcinoma_invivo_infusion_with_glns_m'
    renal_carcinoma_invivo_infusion_with_glns_m_traditional_method = 'renal_carcinoma_invivo_infusion_with_glns_m_traditional_method'
    # 其他类型模型
    lung_tumor_invivo_infusion = 'lung_tumor_invivo_infusion'
    hct116_cultured_cell_line = 'hct116_cultured_cell_line'
    hct116_cultured_cell_line_with_glns_m = 'hct116_cultured_cell_line_with_glns_m'

# 数据模型的详细说明字典
data_model_comment = {
    DataModelType.colon_cancer_cell_line:
        'Data from cultured 8 colon cancer cell line and incubated with 13C in media',
    DataModelType.colon_cancer_cell_line_squared_loss:
        'Data from cultured 8 colon cancer cell line and incubated with 13C in media, optimizing with mean squared loss',
    DataModelType.colon_cancer_cell_line_traditional_method:
        'Data from cultured 8 colon cancer cell line and incubated with 13C in media, '
        'processing with traditional method',
    DataModelType.colon_cancer_cell_line_with_glns_m:
        'Data from cultured 8 colon cancer cell line and incubated with 13C in media. Model includes GLNS_m',
    DataModelType.colon_cancer_cell_line_with_glns_m_traditional_method:
        'Data from cultured 8 colon cancer cell line and incubated with 13C in media, '
        'processing with traditional method. Model includes GLNS_m',
    DataModelType.hct116_cultured_cell_line:
        'Data from cultured cancer cell line (HCT116) and incubated with 13C in media',
    DataModelType.hct116_cultured_cell_line_with_glns_m:
        'Data from cultured cancer cell line (HCT116) and incubated with 13C in media, Model includes GLNS_m',
    DataModelType.renal_carcinoma_invivo_infusion:
        'Data from patient with renal carcinoma and infused with 13C in vivo',
    DataModelType.renal_carcinoma_invivo_infusion_squared_loss:
        'Data from patient with renal carcinoma and infused with 13C in vivo, optimizing with mean squared loss',
    DataModelType.renal_carcinoma_invivo_infusion_traditional_method:
        'Data from patient with renal carcinoma and infused with 13C in vivo, processing with traditional method',
    DataModelType.renal_carcinoma_invivo_infusion_with_glns_m:
        'Data from patient with renal carcinoma and infused with 13C in vivo. Model includes GLNS_m',
    DataModelType.renal_carcinoma_invivo_infusion_with_glns_m_traditional_method:
        'Data from patient with renal carcinoma and infused with 13C in vivo, processing with traditional method. '
        'Model includes GLNS_m',
    DataModelType.lung_tumor_invivo_infusion:
        'Data from patient with lung tumor and infused with 13C in vivo',
}

# 运行模式枚举，定义了不同的分析运行模式
class RunningMode(ValueEnum):
    flux_analysis = 'flux_analysis'
    result_process = 'result_process'
    raw_experimental_data_plotting = 'raw_experimental_data_plotting'
    solver_output = 'solver_output'

