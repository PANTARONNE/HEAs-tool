# HEA DFT 数据库自动化计算工作流

本工作流通过 `scripts/run-hamiltonian-workflow.sh` 一键驱动以下全部步骤，在 Slurm 登录节点运行，自动阻塞等待每一阶段完成后再推进。以下命令均从项目根目录执行：

1. 生成并注册初始 HEA 表面结构（SQS）
2. VASP slab 弛豫 → 收敛判定 → 自动续算 → 弛豫结构入库
3. 表面原子标注（网格化）
4. OpenMX Hamiltonian 计算 → 哈密顿矩阵入库

---

## 依赖与前置条件

| 依赖 | 说明 |
|---|---|
| Python（conda 环境 `ase`） | ASE、icet（SQS）、numpy |
| vaspkit | VASP 输入文件生成 |
| VASP（`vasp_gam`） | slab 弛豫 |
| OpenMX | Hamiltonian SCF |
| Slurm | 作业管理，需要 `sbatch`/`squeue` |

项目根目录下的脚本文件：

```text
scripts/
  vasp-inputs.sh          # 从 CIF 生成 VASP 输入文件（已有）
  vasp-cpu.slurm          # VASP CPU 提交脚本，含自检自投逻辑（默认）
  vasp-dcu.slurm          # VASP DCU/GPU 提交脚本，含自检自投逻辑
  check-convergence.sh    # 收敛判定（纯读取，无副作用）
  openmx.slurm            # OpenMX 提交模板（已有）
    run-hamiltonian-workflow.sh  # 全流程自动化驱动脚本
hea_dataset.py            # 数据集管理工具
cif_to_openmx.py          # CIF -> OpenMX .dat 转换
```

---

## 一键运行（推荐）

在项目根目录的登录节点上执行：

```bash
# 指定元素，比例随机生成
bash scripts/run-hamiltonian-workflow.sh \
  -e Fe Co Ni Cr Mn \
  --openmx-data /work/home/<user>/openmx3.9/DFT_DATA19

# 指定元素和比例（等比例）
bash scripts/run-hamiltonian-workflow.sh \
  -e Fe Co Ni Cr Mn -r 1 1 1 1 1 \
  --openmx-data /work/home/<user>/openmx3.9/DFT_DATA19

# 组成已存在时，跳过第 1 步从弛豫开始续跑
bash scripts/run-hamiltonian-workflow.sh \
  --surface-id Co13Cr13Fe13Mn12Ni13 \
  --openmx-data /work/home/<user>/openmx3.9/DFT_DATA19

# 只运行到第 3 步（不做 OpenMX）
bash scripts/run-hamiltonian-workflow.sh \
  -e Fe Co Ni Cr Mn --skip-hamilton \
  --openmx-data /work/home/<user>/openmx3.9/DFT_DATA19

# 集群环境无 Terminal 输出，后台计算
nohup bash HEA-tools/scripts/run-hamiltonian-workflow.sh \
  -e Fe Co Ni Cu Zn \
  -r 1 1 1 1 1 \
  --openmx-data ~/DFT_DATA19 \
  --root ./dataset \
  --workspace ./workspace \
  > ./logs/FeCoNiCuZn.log 2>&1 &
```

**完整参数列表：**

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-e`/`--elements EL...` | 必填（与 `--surface-id` 二选一） | 元素种类 |
| `-r`/`--ratios R...` | 随机 | 各元素比例，省略则随机抽样 |
| `--surface-id ID` | — | 跳过第 1 步，使用已有组成 |
| `--root DIR` | `dataset` | 数据集根目录 |
| `--workspace DIR` | `workspace` | 计算工作区 |
| `--openmx-data DIR` | `DFT_DATA19` | 集群上 OpenMX DFT_DATA 目录 |
| `--poll SECONDS` | `60` | Slurm 作业状态轮询间隔 |
| `--max-gen N` | `5` | VASP 弛豫最大续算代数 |
| `--skip-hamilton` | — | 第 3 步后停止，不做 OpenMX |

---

## 步骤详解

### 步骤 0：初始化（首次使用）

```bash
python hea_dataset.py init --root dataset
mkdir -p workspace
```

`init` 创建 SQLite 索引和 `dataset_manifest.json`；`create-sample` 也会自动初始化，所以实际上只在需要单独预建数据集目录时使用。

---

### 步骤 1：生成并注册初始 HEA 结构

`run-hamiltonian-workflow.sh` 调用：

```bash
python hea_dataset.py create-sample --root dataset \
  -e Fe Co Ni Cr Mn [-r 1 1 1 1 1]
```

**逻辑：**

1. 按比例（随机或指定）→ Vegard 定律估算晶格常数 → 在 4×4×4 FCC(111) slab（64 个位点）上用最大余数法分配整数原子数，得到候选组成。
2. 检查该 `surface_id` 是否已在数据库或磁盘上存在：
   - 比例随机时撞车 → 重新抽样重试（上限 500 次）。
   - 比例固定时撞车 → 报错退出。
3. 组成唯一后，随机替代原子 → SQS（icet）优化元素分布 → 按元素字母序排列原子 → 写 `00_initial_sqs.cif`。

注册后的目录结构：

```text
dataset/
  index.sqlite
  dataset_manifest.json
  Co13Cr13Fe13Mn12Ni13/
    manifest.json              # status: created
    structures/
      00_initial_sqs.cif       # 初始 SQS 结构
    metadata/                  # 待后续步骤填充
    adsorbates/
    openmx_slab/
```

`surface_id`（如 `Co_13-Cr_13-Fe_13-Mn_12-Ni_13`）由 `run-hamiltonian-workflow.sh` 自动从输出中解析，无需手动传入后续步骤。

---

### 步骤 2：VASP slab 弛豫

#### 2a. 准备计算目录

`run-hamiltonian-workflow.sh` 在 `workspace/<surface_id>/slab-relax/` 中：

1. 将 `00_initial_sqs.cif` 复制为 `<surface_id>.cif`。
2. 运行 `vasp-inputs.sh <surface_id>.cif`，生成：
   - `POSCAR`（含固定底层原子的 F/F/F 标记，`z ∈ [0, 0.46]` 的原子被固定）
   - `INCAR`（来自 `scripts/INCAR-opt` 模板：ENCUT=500, EDIFF=1e-5, EDIFFG=-0.02, NSW=200, IBRION=2）
   - `POTCAR`、`KPOINTS`（Gamma-only）
   - `vasp-cpu.slurm`（默认，job-name 已设为 `surface_id`；可用 `--slurm vasp-dcu.slurm` 切换）
3. 将 `check-convergence.sh` 复制到该目录。
4. 初始化 `.gen_count = 1`，提交第一代作业。

#### 2b. 自检自投链（`vasp-cpu.slurm` / `vasp-dcu.slurm`）

每次 VASP 运行结束后，slurm 脚本自动执行以下逻辑：

```
mpirun 完成
    ↓
把本代输出备份到 gen_01/ (gen_02/ ...)
    ↓
bash check-convergence.sh
    ↓
退出码 0  → 写 CONVERGED，链结束
退出码 10 → 代数 < MAX_GEN? → cp CONTCAR POSCAR → 代数+1 → sbatch 自己
                             否? → 写 NOT_CONVERGED_MAX_GEN，链结束
退出码 20 → 写 CRASHED，链结束（超时/崩溃，禁止续算）
退出码 30 → 写 SCF_ILL，链结束（SCF 病态，禁止续算）
```

**`check-convergence.sh` 判定逻辑（4 关依次通过）：**

| 关卡 | 信号 | 退出码 |
|---|---|---|
| VASP 是否正常收尾 | `OUTCAR` 末尾有 `General timing and accounting` | 无则 → 20 |
| 离子弛豫是否收敛 | `OUTCAR` 有 `reached required accuracy` | 有则 → 0 |
| SCF 是否病态 | `OSZICAR` 中**连续 ≥ 3 个**离子步的电子迭代数撞 NELM | 是则 → 30 |
| 其余情况 | 正常收尾，SCF 健康，但未收敛 | → 10 |

SCF 病态阈值可通过环境变量 `SCF_STALL_LIMIT` 覆盖（默认 3）。

每代的输出文件归档在 `gen_01/`、`gen_02/` 等子目录，便于追溯历史。

#### 2c. 收敛后入库

`run-hamiltonian-workflow.sh` 检测到 `CONVERGED` 标记后：

1. 用 ASE 将 `CONTCAR` 转为 CIF 格式（保留 POSCAR 的元素顺序）。
2. 调用 `record-relaxed`，校验组成一致、原子数和元素顺序不变后入库：

```bash
python hea_dataset.py record-relaxed \
  --root dataset \
  --surface-id <surface_id> \
  --relaxed-cif workspace/<surface_id>/slab-relax/01_relaxed_slab.cif
```

注册后 manifest 中 `status` 更新为 `slab_relaxed`，`01_relaxed_slab.cif` 落盘：

```text
dataset/<surface_id>/structures/
  00_initial_sqs.cif
  01_relaxed_slab.cif    # 新增
```

---

### 步骤 3：表面原子标注

```bash
python hea_dataset.py index-surface \
  --root dataset \
  --surface-id <surface_id>
```

从 `00_initial_sqs.cif`（定义原子 ID）和 `01_relaxed_slab.cif`（提供弛豫坐标）中识别最外层金属原子，并将其排列为二维网格，供后续 Hamiltonian 提取时规范化索引使用。生成：

```text
dataset/<surface_id>/metadata/
  top_atoms.jsonl     # 每个顶层原子的 atom_id、元素、初始/弛豫坐标、网格位置
  atom_grid.npy       # atom_id 的二维数组
```

---

### 步骤 4：Hamiltonian 计算与入库

`run-hamiltonian-workflow.sh` 在 `workspace/<surface_id>/hamilton/` 中执行：

**4a. 生成 OpenMX 输入：**

```bash
python cif_to_openmx.py \
  dataset/<surface_id>/structures/01_relaxed_slab.cif \
  --data-path /work/.../DFT_DATA19 \
  -o workspace/<surface_id>/hamilton/input.dat
```

固定使用 `input.dat` 文件名，使 `System.Name = input`，因此 OpenMX 输出为 `input.scfout`。

**4b. 提交 OpenMX 作业：**

`run-hamiltonian-workflow.sh` 从 `scripts/openmx.slurm` 模板自动替换：
- 任务名 → `<surface_id>-hamilton`
- 输入文件 → `openmx input.dat`

提交并等待作业离队（通过 `squeue` 轮询）。

**4c. 提取 Hamiltonian 并入库：**

```bash
python hea_dataset.py extract-hamiltonian \
  --root dataset \
  --surface-id <surface_id> \
  --dat workspace/<surface_id>/hamilton/input.dat \
  --scfout workspace/<surface_id>/hamilton/input.scfout
```

输出：

```text
dataset/<surface_id>/openmx_slab/
  hamiltonian_d_surface.npz          # 顶层原子 d 轨道 Hamiltonian
  hamiltonian_d_surface.npz.basis.jsonl
```

---

## 产生的完整目录结构

```text
dataset/
  index.sqlite
  dataset_manifest.json
  <surface_id>/
    manifest.json
    structures/
      00_initial_sqs.cif
      01_relaxed_slab.cif
    metadata/
      top_atoms.jsonl
      atom_grid.npy
    openmx_slab/
      hamiltonian_d_surface.npz
      hamiltonian_d_surface.npz.basis.jsonl

workspace/
  <surface_id>/
    slab-relax/
      INCAR  POSCAR  POTCAR  KPOINTS
      vasp-cpu.slurm
      check-convergence.sh
      .gen_count
      CONTCAR  OUTCAR  OSZICAR  ...
      01_relaxed_slab.cif       # CONTCAR 转换结果
      CONVERGED                 # 收敛标记（或 CRASHED 等）
      gen_01/                   # 第 1 代归档
        OUTCAR  OSZICAR  CONTCAR  out_<jid>.log  ...
      gen_02/                   # 第 2 代归档（如有续算）
    hamilton/
      input.dat
      openmx.slurm
      input.scfout
      input.std                 # OpenMX 标准输出
```

---

## 故障排查

**弛豫未收敛（`NOT_CONVERGED_MAX_GEN`）：**

查看最后一代的输出，确认力是否在下降：

```bash
grep 'RMS' workspace/<surface_id>/slab-relax/OUTCAR | tail -20
```

若趋势良好，可增大 `--max-gen` 继续；若力停滞，检查结构或 INCAR 参数。

**VASP 崩溃（`CRASHED`）：**

```bash
cat workspace/<surface_id>/slab-relax/gen_01/err_<jid>.log
tail -50 workspace/<surface_id>/slab-relax/gen_01/OUTCAR
```

常见原因：超时（增大 `--time`）、内存不足（调整 `NCORE`/`NPAR`）、POTCAR 缺失。

**SCF 病态（`SCF_ILL`）：**

```bash
grep 'DAV\|RMM' workspace/<surface_id>/slab-relax/OSZICAR | tail -100
```

连续多步电子迭代到 NELM 而不收敛，通常需要调整 `ALGO`、`SIGMA` 或初始磁矩。

**OpenMX 无 scfout：**

```bash
cat workspace/<surface_id>/hamilton/input.std
```

检查 DATA.PATH 是否指向集群上实际存在的目录，以及 OpenMX 模块是否正确加载。
