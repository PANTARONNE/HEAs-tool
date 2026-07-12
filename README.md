# HEA 数据集与计算工作流

本项目用于生成高熵合金（HEA）FCC(111) 表面结构、管理结构与吸附数据集，并驱动 VASP/OpenMX 计算及 Hamiltonian 提取。

## 快速入口

- [Hamiltonian 计算工作流](docs/CALCULATION_WORKFLOW_HAMILTONIAN.md)：结构生成、slab 弛豫、表面索引、OpenMX 计算及 Hamiltonian 入库。
- [吸附能计算工作流](docs/CALCULATION_WORKFLOW_ADSORPTION.md)：FCC 位点识别、吸附结构生成、VASP 弛豫及吸附能入库。
- [数据集命令手册](docs/DATASET_WORKFLOW.md)：`hea_dataset.py` 各子命令的完整用法。
- [最小输入工作流](docs/MINIMAL_INPUT_WORKFLOW.md)：以最少人工输入运行各步骤。

文档中的命令默认从项目根目录执行。

## 目录说明

```text
HEAs/
├── hea_dataset.py                  # 数据集管理主入口
├── build_hea_surface.py            # HEA FCC(111) SQS 表面生成
├── add_fcc_adsorbate.py            # FCC 位点检测与吸附物放置
├── cif_to_openmx.py                # CIF 转 OpenMX 输入
├── extract_openmx_hamiltonian.py   # OpenMX Hamiltonian 提取
├── plot_hamiltonian_heatmap.py     # Hamiltonian 热图绘制
├── docs/                           # 工作流与使用文档
├── scripts/                        # 自动化、Slurm 与 VASP 辅助脚本
├── DFT_DATA19/                     # OpenMX 2019 PAO/VPS 数据库
├── dataset/                        # 主数据集
├── equalRatio/                     # 等比例样本数据集
└── randomRatio/                    # 随机比例样本数据集
```

`dataset/`、`equalRatio/` 和 `randomRatio/` 是相互独立的 SQLite 数据集根目录，不应只移动其中的索引或样本子目录。

## 常用命令

```bash
# 查看数据集管理命令
python hea_dataset.py --help

# 运行结构弛豫与 Hamiltonian 工作流
bash scripts/run-hamiltonian-workflow.sh --help

# 运行吸附能工作流
bash scripts/run-adsorption-workflow.sh --help
```

Python 环境定义见 `environment.yml`。VASP/OpenMX 与 Slurm 环境需要按集群实际配置准备。
