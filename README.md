# MFA Tutorial 2025 - 代谢通量分析教程

<div>
  <h1>🧬 代谢通量分析(MFA)系统</h1>
  <p>基于同位素标记实验数据的代谢通量分析工具</p>
  
  [![Python Version](https://img.shields.io/badge/python-3.9+-blue.svg)](https://python.org)
  [![Status](https://img.shields.io/badge/status-active-brightgreen.svg)]()
</div>

## 📋 目录

- [项目简介](#项目简介)
- [核心特性](#核心特性)
- [安装要求](#安装要求)
- [快速开始](#快速开始)
- [使用示例](#使用示例)
- [项目架构](#项目架构)
- [API文档](#api文档)
- [贡献指南](#贡献指南)
- [许可证](#许可证)

## 🧪 项目简介

**MFA Tutorial 2025** 是一个专业的代谢通量分析(Metabolic Flux Analysis, MFA)系统，采用基本代谢单元(Elementary Metabolite Unit, EMU)方法，通过同位素标记实验数据来量化细胞内代谢反应的通量分布。

### 🎯 主要用途

- **代谢网络建模**：构建和分析复杂的代谢网络模型
- **通量估计**：基于实验数据精确估计代谢反应速率
- **同位素追踪**：追踪¹³C标记在代谢网络中的传播
- **实验数据集成**：处理质谱仪测量的MID(Mass Isotopomer Distribution)数据
- **结果可视化**：生成详细的分析报告和图表

## ✨ 核心特性

### 🔬 先进的算法实现
- **EMU方法**：高效的同位素标记模式追踪
- **SLSQP优化**：序列最小二乘规划算法优化通量值
- **矩阵方程求解**：快速求解大规模线性方程组
- **Numba加速**：关键计算部分使用JIT编译加速

### 📊 数据处理能力
- **多格式数据支持**：支持Excel、CSV等多种数据格式
- **实验数据验证**：自动检查数据完整性和一致性
- **噪声处理**：内置噪声模拟和数据清洗功能
- **批量分析**：支持多组实验条件的批量处理

### 🧮 建模功能
- **代谢网络构建**：直观的代谢物和反应定义
- **对称性处理**：自动处理分子对称性(如琥珀酸、延胡索酸)
- **约束管理**：灵活的通量边界和质量平衡约束
- **模型验证**：自动检查模型的数学一致性

### 📈 分析工具
- **敏感性分析**：计算通量控制系数
- **置信区间**：估计参数不确定性
- **模型比较**：支持多个模型的统计比较
- **结果可视化**：生成专业的图表和报告

## 🛠 安装要求

### 系统要求
- Python 3.8+
- Windows 10/11, macOS 10.14+, 或 Linux

### 依赖包
```bash
# 核心计算库
numpy>=1.21.0
scipy>=1.7.0
pandas>=1.3.0

# 优化和加速
numba>=0.56.0
scikit-learn>=1.0.0

# 数据处理
openpyxl>=3.0.0
xlsxwriter>=3.0.0

# 可视化(可选)
matplotlib>=3.5.0
seaborn>=0.11.0
plotly>=5.0.0
```

### 安装步骤

1. **克隆仓库**
```bash
git clone https://github.com/yourusername/tutorial-mfa-2025.git
cd tutorial-mfa-2025
```

2. **创建虚拟环境**
```bash
python -m venv mfa_env
source mfa_env/bin/activate  # Linux/macOS
# 或
mfa_env\Scripts\activate     # Windows
```

3. **安装依赖**
```bash
pip install -r requirements.txt
```

## 🚀 快速开始

### 基础示例：TCA循环分析

```python
import numpy as np
from scripts.src.core.model.model_class import Reaction, Node
from common_test_pipeline import run_mfa_pipeline_with_experimental_data

# 1. 定义代谢网络
def setup_tca_network():
    # 定义反应：OAC + AcCoA → Cit
    v1 = Reaction(
        id='v1',
        sub=[
            Node(name='OAC', coefficient=1, carbon_composition_string=['a','b','c','d']),
            Node(name='AcCoA', coefficient=1, carbon_composition_string=['e','f'])
        ],
        pro=[
            Node(name='Cit', coefficient=1, carbon_composition_string=['d','c','b','f','e','a'])
        ],
        reverse=False
    )
    # ... 定义更多反应
    
    return metabolite_reaction_dict, flux_balance_reaction_dict, \
           complete_metabolite_dim_dict, flux_name_index_dict

# 2. 运行MFA分析
metabolite_reaction_dict, flux_balance_reaction_dict, \
complete_metabolite_dim_dict, flux_name_index_dict = setup_tca_network()

# 设置初始通量值
flux_vector = np.array([100, 100, 50, 50, 50, 125, 75, 50], dtype=float)

# 设置输入代谢物标记
input_metabolite_dict = {
    "AcCoA": [
        {"ratio_list": [0, 1], "abundance": 0.25},  # [2-13C]AcCoA
        {"ratio_list": [1, 1], "abundance": 0.25},  # [1,2-13C2]AcCoA  
        {"ratio_list": [0, 0], "abundance": 0.5}    # 无标记
    ]
}

# 运行完整分析
results = run_mfa_pipeline_with_experimental_data(
    metabolite_reaction_dict=metabolite_reaction_dict,
    input_metabolite_name_set={'AcCoA', 'Asp'},
    complete_metabolite_dim_dict=complete_metabolite_dim_dict,
    target_metabolite_name_list=['Glu'],
    flux_name_index_dict=flux_name_index_dict,
    flux_vector=flux_vector,
    input_metabolite_dict=input_metabolite_dict,
    # ... 其他参数
)

# 查看结果
print("优化后的通量值:")
for flux_name, idx in flux_name_index_dict.items():
    print(f"{flux_name}: {results['optimized_flux_vector'][idx]:.4f}")
```

### 实验数据集成示例

```python
from common_test_pipeline import load_hct116_experimental_data

# 加载HCT116细胞系实验数据
experimental_data, metabolite_mapping = load_hct116_experimental_data(
    experiment_name="HCT116_WQ2101",
    condition="ctrl", 
    index="average"
)

# 定义代谢物映射
user_metabolite_mapping = {
    'glutamate': 'Glu',
    'citrate': 'Cit', 
    'succinate': 'Suc',
    # ... 更多映射
}

# 运行集成分析
results = run_mfa_pipeline_with_experimental_data(
    # ... 模型参数
    experimental_mid_data_obj_dict=experimental_data,
    model_metabolite_to_standard_name_dict=user_metabolite_mapping,
    verbose=True
)
```

## 📁 项目架构

```
tutorial-mfa-2025/
├── 📂 scripts/                    # 核心代码
│   ├── 📂 src/
│   │   ├── 📂 core/               # 核心功能
│   │   │   ├── 📂 common/         # 通用配置和工具
│   │   │   ├── 📂 model/          # 代谢模型相关
│   │   │   ├── 📂 data/           # 数据处理
│   │   │   └── 📂 solver/         # 优化求解器
│   │   ├── 📂 experimental_data_analysis/  # 实验数据分析
│   │   └── 📂 simulated_data/     # 模拟数据生成
│   ├── 📂 data/                   # 实验数据
│   │   └── 📂 hct116_cultured_cell_line/
│   └── 📂 model/                  # 模型定义
├── 📄 common_functions.py         # 通用函数
├── 📄 common_test_pipeline.py     # 测试流程
├── 📄 test_tca.py                # TCA循环示例
├── 📄 test_sample.py             # 样本测试
└── 📄 architecture_explanation.md # 架构说明
```

### 核心组件说明

| 组件 | 功能描述 |
|------|----------|
| **Model** | 代谢网络建模、EMU分析、矩阵方程生成 |
| **Solver** | SLSQP优化、EMU图构建、预测计算 |
| **Data** | 实验数据加载、MID数据处理、格式转换 |
| **Analysis** | 结果分析、敏感性分析、可视化 |

## 📚 API文档

### 主要函数

#### `run_mfa_pipeline()`
运行基础MFA分析流程

**参数:**
- `metabolite_reaction_dict`: 代谢物反应字典
- `input_metabolite_name_set`: 输入代谢物集合  
- `complete_metabolite_dim_dict`: 代谢物碳原子数字典
- `target_metabolite_name_list`: 目标代谢物列表
- `flux_name_index_dict`: 通量名称索引字典
- `flux_vector`: 初始通量向量
- `input_metabolite_dict`: 输入代谢物标记字典

**返回:**
- `dict`: 包含EMU方程、预测MID数据等结果的字典

#### `run_mfa_pipeline_with_experimental_data()`
集成实验数据的MFA分析流程

**额外参数:**
- `experimental_mid_data_obj_dict`: 实验MID数据字典
- `model_metabolite_to_standard_name_dict`: 代谢物映射字典
- `specific_flux_range_dict`: 通量范围约束字典

**返回:**
- `dict`: 包含优化结果、比较数据、控制系数等的完整结果字典

### 模型定义类

#### `Reaction`
代谢反应类

```python
Reaction(
    id='reaction_id',
    sub=[Node(...)],      # 底物列表
    pro=[Node(...)],      # 产物列表  
    reverse=False         # 是否可逆
)
```

#### `Node` 
代谢物节点类

```python
Node(
    name='metabolite_name',
    coefficient=1.0,
    carbon_composition_string=['a','b','c']  # 碳原子标记
)
```

## 🔬 使用示例

### 1. 简单代谢路径分析

```python
# 定义简单的三步代谢路径: A → B → C → D
from scripts.src.core.model.model_class import Reaction, Node

reactions = {
    'A': [],  # 输入
    'B': [Reaction(id='v1', sub=[Node('A', 1, ['a','b'])], 
                   pro=[Node('B', 1, ['a','b'])])],
    'C': [Reaction(id='v2', sub=[Node('B', 1, ['a','b'])], 
                   pro=[Node('C', 1, ['a','b'])])],
    'D': [Reaction(id='v3', sub=[Node('C', 1, ['a','b'])], 
                   pro=[Node('D', 1, ['a','b'])])]
}

# 运行分析
results = run_mfa_pipeline(
    metabolite_reaction_dict=reactions,
    input_metabolite_name_set={'A'},
    complete_metabolite_dim_dict={'A':2, 'B':2, 'C':2, 'D':2},
    target_metabolite_name_list=['D'],
    flux_name_index_dict={'v1':0, 'v2':1, 'v3':2},
    flux_vector=np.array([100, 100, 100])
)
```

### 2. 分支代谢网络

```python
# 处理代谢分支: A → B → C
#                    ↘ D
branch_reactions = {
    'A': [],
    'B': [Reaction(id='v1', sub=[Node('A', 1, ['a','b'])], 
                   pro=[Node('B', 1, ['a','b'])])],
    'C': [Reaction(id='v2', sub=[Node('B', 1, ['a','b'])], 
                   pro=[Node('C', 1, ['a','b'])])],
    'D': [Reaction(id='v3', sub=[Node('B', 1, ['a','b'])], 
                   pro=[Node('D', 1, ['a','b'])])]
}

# 设置通量比例: 70%流向C, 30%流向D  
flux_vector = np.array([100, 70, 30])
```

### 3. 处理对称分子

```python
# 琥珀酸等对称分子的处理
symmetric_reaction = Reaction(
    id='v_symmetric',
    sub=[Node('Suc', 1, ['a','b','c','d'])],
    pro=[
        Node('Fum', 0.5, ['a','b','c','d']),  # 50%保持原序
        Node('Fum', 0.5, ['d','c','b','a'])   # 50%翻转序列
    ]
)
```

## 📊 结果解读

### 通量估计结果
```python
# 查看优化后的通量值
optimized_flux = results['optimized_flux_vector']
print("反应通量估计:")
for reaction, flux in zip(flux_names, optimized_flux):
    print(f"{reaction}: {flux:.2f} ± {uncertainty:.2f}")
```

### MID预测vs实验比较
```python
# 分析预测准确性
comparison_data = results['comparison_data']
for data in comparison_data:
    print(f"\n{data['experimental_metabolite']}:")
    print(f"RMSE: {data['RMSE']:.4f}")
    print("M+0  M+1  M+2  M+3  M+4")
    print("预测:", " ".join([f"{x:.3f}" for x in data['predicted_mid']]))
    print("实验:", " ".join([f"{x:.3f}" for x in data['experimental_mid']]))
```

### 敏感性分析
```python
# 查看通量控制系数
flux_control = results['flux_control_coefficients']
print("\n通量控制系数(对1%扰动的响应):")
for flux_name, coefficients in flux_control.items():
    print(f"\n{flux_name}:")
    for metabolite, coefficient in coefficients.items():
        print(f"  {metabolite}: {coefficient:.4f}")
```

## 🤝 贡献指南

我们欢迎社区贡献！请遵循以下步骤：

1. **Fork** 此仓库
2. 创建特性分支 (`git checkout -b feature/AmazingFeature`)
3. 提交更改 (`git commit -m 'Add some AmazingFeature'`)
4. 推送到分支 (`git push origin feature/AmazingFeature`)
5. 打开 **Pull Request**

### 代码规范
- 遵循PEP 8 Python代码风格
- 添加适当的文档字符串
- 编写单元测试
- 更新相关文档

### 报告问题
使用GitHub Issues报告bug或请求新功能，请包含：
- 详细的问题描述
- 重现步骤
- 系统环境信息
- 相关代码片段

## 📄 许可证

本项目采用许可证 - 查看 [LICENSE](LICENSE) 文件了解详情。

## 🙏 致谢

- 感谢所有贡献者的努力
- 基于科学界在代谢流分析领域的前沿研究
- 特别感谢EMU方法的开发者们

## 📞 联系方式

- **项目主页**: https://github.com/yourusername/tutorial-mfa-2025
- **文档**: https://tutorial-mfa-2025.readthedocs.io
- **问题反馈**: https://github.com/yourusername/tutorial-mfa-2025/issues

---

<div align="center">
  <p>如果这个项目对您有帮助，请给我们一个 ⭐!</p>
</div>

