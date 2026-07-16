# HEA 吸附能自动化计算工作流

本工作流在 [Hamiltonian 计算工作流](CALCULATION_WORKFLOW_HAMILTONIAN.md) 已完成的基础上（slab 已弛豫且总能量已入库），针对**指定吸附物种**计算其在各 FCC 位点的吸附能。通过 `scripts/run-adsorption-workflow.sh` 一键驱动，在 Slurm 登录节点运行，自动阻塞等待每一阶段完成后再推进。以下命令均从项目根目录执行：

1. 表面位点识别（FCC hollow 位点枚举）
2. 吸附结构生成（每位点一个初始吸附构型）
3. 对每个初始吸附结构进行 VASP 结构优化 → 收敛判定 → 自动续算
4. 计算吸附能并入库

**必需输入：** `--surface-id`（已完成工作流 1 的组成）与 `-a`（吸附物种）。

---

## 依赖与前置条件

| 依赖 | 说明 |
|---|---|
| Python（conda 环境 `ase`） | ASE、numpy |
| vaspkit | VASP 输入文件生成 |
| VASP（`vasp_gam`） | slab+吸附物弛豫 |
| Slurm | 作业管理，需要 `sbatch`/`squeue` |

**数据集前置状态（由工作流 1 产生）：**

- `dataset/<surface_id>/structures/01_relaxed_slab.cif` 已存在（slab 已弛豫入库）。
- `surfaces.total_energy_eV` 已记录（工作流 1 的 `record-relaxed --energy`）。若缺失，可用 `--slab-energy` 手动传入。
- `dataset/<surface_id>/metadata/top_atoms.jsonl` 已存在（工作流 1 的 `index-surface`），使位点能关联到顶层原子 ID。

项目根目录下相关脚本：

```text
scripts/
  run-adsorption-workflow.sh   # 吸附能全流程自动化驱动脚本
  record-single-site-adsorption.sh  # 手动补记单个位点的吸附能
  vasp-inputs.sh               # 从 CIF 生成 VASP 输入文件（复用）
  vasp-gam.slurm               # VASP 提交脚本，含自检自投逻辑（复用）
  check-convergence.sh         # 收敛判定（纯读取，无副作用，复用）
hea_dataset.py                 # 数据集管理工具
add_fcc_adsorbate.py           # FCC 位点检测与吸附物放置
```

---

## 一键运行（推荐）

> **运行前必读：** 打开 `scripts/run-adsorption-workflow.sh`，在脚本顶部的
> `ADSORBATE_REF_ENERGY` 数组中填入本次所用物种的参考能量（eV），否则脚本会报错退出。

```bash
# 计算 N 在所有 FCC 位点的吸附能（比例已在工作流 1 确定）
bash scripts/run-adsorption-workflow.sh \
  --surface-id Co_13-Cr_13-Fe_13-Mn_12-Ni_13 \
  -a N

# 只计算指定位点（例如第 5 个）
bash scripts/run-adsorption-workflow.sh \
  --surface-id Co_13-Cr_13-Fe_13-Mn_12-Ni_13 \
  -a NH2 --site 5

# 计算多个指定位点
bash scripts/run-adsorption-workflow.sh \
  --surface-id Co_13-Cr_13-Fe_13-Mn_12-Ni_13 \
  -a NH --site 1 3 5

# 命令行临时提供参考能量与 slab 能量（不改脚本）
bash scripts/run-adsorption-workflow.sh \
  --surface-id Co_13-Cr_13-Fe_13-Mn_12-Ni_13 \
  -a N --ref-energy -16.63 --slab-energy -412.35

# 集群环境后台运行，日志重定向
nohup bash scripts/run-adsorption-workflow.sh \
  --surface-id Co_13-Cr_13-Fe_13-Mn_12-Ni_13 \
  -a N \
  --root ./dataset \
  --workspace ./workspace \
  > ./logs/Co-Cr-Fe-Mn-Ni_N.log 2>&1 &
```

**完整参数列表：**

| 参数 | 默认值 | 说明 |
|---|---|---|
| `--surface-id ID` | 必填 | 已完成工作流 1 的组成 ID |
| `-a`/`--adsorbate SP` | 必填 | 吸附物种：`N` \| `NH` \| `NH2` \| `NH3` |
| `--site SEL...` | 全部 | 位点编号（如 `--site 1 3 5`），省略或 `all` 表示全部位点 |
| `--root DIR` | `dataset` | 数据集根目录 |
| `--workspace DIR` | `workspace` | 计算工作区 |
| `--slab-energy E` | 从数据集读取 | 手动覆盖 E(slab)，单位 eV |
| `--ref-energy E` | 从脚本内数组读取 | 手动覆盖本次吸附物参考能量，单位 eV |
| `--height A` | `1.25` | 吸附物 N 原子距 hollow 平面高度（Å） |
| `--nh A` | `1.02` | NHx 的 N–H 键长（Å） |
| `--poll SECONDS` | `60` | Slurm 作业状态轮询间隔 |
| `--max-gen N` | `5` | 每个位点 VASP 弛豫最大续算代数 |
| `--skip-detect` | — | 复用已有 `fcc_sites.jsonl`，跳过第 1 步 |
| `--force` | — | 对已有终止标记的位点强制重投 |

---

## 吸附能定义

脚本按下式计算吸附能：

```
E_ads = E(slab + 吸附物)  -  E(slab)  -  E_ref(吸附物)
```

- **E(slab + 吸附物)：** 本工作流第 3 步 VASP 弛豫得到的最终总能量（取 OUTCAR 的
  `energy without entropy`，与工作流 1 一致）。
- **E(slab)：** 工作流 1 登记的干净 slab 总能量，默认从 `index.sqlite` 的
  `surfaces.total_energy_eV` 读取。
- **E_ref(吸附物)：** 吸附物参考（如气相）总能量，**由使用者手动填入**脚本顶部的
  `ADSORBATE_REF_ENERGY` 数组：

  ```bash
  declare -A ADSORBATE_REF_ENERGY=(
      [N]=""      # eV  例如 0.5 * E(N2)
      [NH]=""     # eV
      [NH2]=""    # eV
      [NH3]=""    # eV
  )
  ```

  只需填写本次用到的物种；缺失且未提供 `--ref-energy` 时脚本报错退出。

---

## 步骤详解

### 步骤 1：表面位点识别

`run-adsorption-workflow.sh` 调用：

```bash
python hea_dataset.py detect-sites \
  --root dataset --surface-id <surface_id>
```

**逻辑：** 从 `01_relaxed_slab.cif` 沿表面法向识别层结构，用第三层原子的面内位置标记 FCC hollow 位点，再将位点吸附到最近三个顶层原子构成的 hollow 中心（`relaxed` 模式，对弛豫后的表面更鲁棒）。4×4 表面得到 **16 个位点**。生成：

```text
dataset/<surface_id>/metadata/
  fcc_sites.jsonl    # 每个位点的 site_id、site_index、行列、面内/平面坐标、关联顶层原子 ID
  site_grid.npy      # site_id 的二维数组
```

指定 `--skip-detect` 且文件已存在时跳过本步。位点检测是确定性的，重跑得到相同的 `site_index`。

---

### 步骤 2：吸附结构生成

```bash
python hea_dataset.py create-adsorbate-records \
  --root dataset --surface-id <surface_id> \
  --adsorbate <SP> --height 1.25 --nh 1.02
```

**逻辑：** 读取 `01_relaxed_slab.cif` 与 `fcc_sites.jsonl`，为**每个位点**在 hollow 上方放置吸附物（N 原子朝向表面，NHx 的 H 朝外），写出初始吸附结构并登记记录：

```text
dataset/<surface_id>/adsorbates/<SP>/site_XXXX/
  00_initial_adsorbate.cif    # 初始吸附结构（slab 原子在前，吸附物 N/H 追加在后）
  adsorption_energy.json      # 记录占位：energy_status = empty
```

原子顺序为 `[金属(字母序)…, N, H…]`，与 vaspkit 按元素分组的顺序一致，从而保证弛豫后 `record-energy` 的原子数/元素顺序校验通过。

---

### 步骤 3：slab+吸附物 VASP 弛豫

对每个选定位点，`run-adsorption-workflow.sh` 在
`workspace/<surface_id>/adsorbate/<SP>/site_XXXX/` 中：

#### 3a. 准备计算目录

1. 将 `00_initial_adsorbate.cif` 复制为 `<surface_id>_<SP>_site_XXXX.cif`（作为 Slurm 作业名，无空格）。
2. 运行 `vasp-inputs.sh`，生成 `POSCAR`（`z ∈ [0, 0.46]` 的底层原子被固定，吸附物在顶部不受约束）、`INCAR`、`POTCAR`、`KPOINTS`（Gamma-only）及 `vasp-gam.slurm`。
3. 复制 `check-convergence.sh`，初始化 `.gen_count = 1`，提交第一代作业。

**所有选定位点的作业并行提交**，随后脚本统一轮询各位点目录的终止标记，全部离队后再推进。

#### 3b. 自检自投链（`vasp-gam.slurm`）

与工作流 1 完全相同：每次 `vasp_gam` 结束后备份本代输出到 `gen_01/`（`gen_02/`…），运行 `check-convergence.sh` 判定：

```
退出码 0  → 写 CONVERGED，链结束
退出码 10 → 代数 < MAX_GEN? → cp CONTCAR POSCAR → 代数+1 → sbatch 自己
                             否? → 写 NOT_CONVERGED_MAX_GEN，链结束
退出码 20 → 写 CRASHED，链结束（超时/崩溃）
退出码 30 → 写 SCF_ILL，链结束（SCF 病态）
```

`--force` 可对已有终止标记的位点强制重投；缺省则跳过已有 `CONVERGED`（等）的位点，便于中断续跑。

---

### 步骤 4：计算吸附能并入库

对每个 `CONVERGED` 的位点：

1. 用 ASE 将 `CONTCAR` 转为 `01_relaxed_adsorbate.cif`（保留元素顺序）。
2. 从 `OUTCAR` 提取最终总能量（`energy without entropy`）。
3. 计算 `E_ads = E(slab+ads) − E(slab) − E_ref`。
4. 调用 `record-energy` 入库：

```bash
python hea_dataset.py record-energy \
  --root dataset --surface-id <surface_id> \
  --adsorbate <SP> --site-id <N> \
  --relaxed-cif .../01_relaxed_adsorbate.cif \
  --energy <E_ads> --status computed \
  --notes "auto: E_slab=... eV, E_ref=... eV"
```

`record-energy` 会校验：relaxed 结构原子数/元素顺序不变、子结构与 `01_relaxed_slab.cif` 一致、吸附物位移不超过 2.0 Å（超出则报错，不入库）。入库后更新：

```text
dataset/<surface_id>/adsorbates/<SP>/site_XXXX/
  01_relaxed_adsorbate.cif          # 新增（弛豫结果）
  adsorption_energy.json            # adsorption_energy_eV、energy_status=computed、relaxation_validation
```

同时更新 `adsorbate_configs` 与 `fcc_sites` 表中该位点的吸附能与状态。

**未收敛 / 无 CONTCAR / 无能量的位点自动跳过，不写入垃圾数据。** 运行结束打印 `recorded`（入库数）与 `skipped`（跳过数）。

---

## 产生的完整目录结构

```text
dataset/
  index.sqlite
  <surface_id>/
    structures/
      00_initial_sqs.cif
      01_relaxed_slab.cif
    metadata/
      top_atoms.jsonl
      atom_grid.npy
      fcc_sites.jsonl               # 步骤 1 产生（或已存在）
      site_grid.npy
    adsorbates/
      <SP>/
        site_0001/
          00_initial_adsorbate.cif  # 步骤 2
          01_relaxed_adsorbate.cif  # 步骤 4
          adsorption_energy.json    # 步骤 2 创建 → 步骤 4 填充
        site_0002/
        ...

workspace/
  <surface_id>/
    adsorbate/
      <SP>/
        site_0001/
          INCAR  POSCAR  POTCAR  KPOINTS
          vasp-gam.slurm
          check-convergence.sh
          .gen_count
          CONTCAR  OUTCAR  OSZICAR  ...
          01_relaxed_adsorbate.cif  # CONTCAR 转换结果
          CONVERGED                 # 收敛标记（或 CRASHED 等）
          gen_01/                   # 第 1 代归档
          gen_02/                   # 第 2 代归档（如有续算）
        site_0002/
        ...
```

---

## 手动补记单个位点

若自动工作流将某个位点标记为 `CRASHED`，但手动恢复计算后已经在该位点目录得到
有效的 `CONTCAR` 和 `OUTCAR`，可直接计算并记录该点，无需修改或移除终止标记：

```bash
bash scripts/record-single-site-adsorption.sh <surface_id> <site编号>
```

脚本会从标准的 `workspace`/`dataset` 目录自动识别吸附物种，从数据集读取
`E(slab)`，从 `OUTCAR` 读取最后一个 `energy without entropy`，计算吸附能并调用
`record-energy` 入库。如果同一位点存在多个吸附物种，请用 `-a N`（或 `NH`、
`NH2`、`NH3`）明确指定。手动计算目录不在标准位置时可传入 `--calc-dir DIR`；
其余选项可运行 `bash scripts/record-single-site-adsorption.sh --help` 查看。

---

## 故障排查

**报错 `no reference energy for <SP>`：**

在 `scripts/run-adsorption-workflow.sh` 顶部的 `ADSORBATE_REF_ENERGY` 中填入该物种的参考能量，或运行时加 `--ref-energy <值>`。

**报错 `E(slab) is not recorded`：**

该组成在工作流 1 中未记录 slab 总能量。重新运行工作流 1 的 `record-relaxed --energy`，或本次运行加 `--slab-energy <值>`。

**报错 `relaxed slab not found`：**

尚未完成 Hamiltonian 计算工作流。先运行 `scripts/run-hamiltonian-workflow.sh` 得到 `01_relaxed_slab.cif`。

**某位点弛豫未收敛（`NOT_CONVERGED_MAX_GEN`）：**

```bash
grep 'RMS' workspace/<surface_id>/adsorbate/<SP>/site_XXXX/OUTCAR | tail -20
```

若力仍在下降，增大 `--max-gen` 后加 `--force` 重投该位点。

**某位点被跳过入库：**

查看该位点目录内的终止标记（`CRASHED`/`SCF_ILL` 等）与 `gen_01/err_*.log`；崩溃常见原因为超时、内存不足或 POTCAR 缺失。

**`record-energy` 报吸附物位移过大：**

吸附物在弛豫中脱附或迁移到其他位点。检查 `01_relaxed_adsorbate.cif`；如属正常物理过程可用 `--max-adsorbate-displacement` 放宽后手动调用 `record-energy` 入库。

**`record-energy` 报原子顺序不一致：**

vaspkit 生成 POSCAR 时的元素顺序与初始 CIF 不符（少见）。核对 `POSCAR` 第 6 行元素顺序与 `00_initial_adsorbate.cif` 是否一致。
