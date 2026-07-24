# 条目完整性检查（check）

`python hea_dataset.py check` 用于检查指定 `surface_id` 的条目是否完整，覆盖从初始结构到吸附能的全部工作流产物。检查同时校验 **SQLite 索引** 与 **磁盘文件**，两者不一致时会明确标出，便于发现孤立的索引行或未入库的文件。

全部完整时进程退出码为 `0`，否则为 `1`，可直接用于自动化流程的条件判断。

---

## 用法

```bash
# 仅检查核心产物
python hea_dataset.py check --root dataset --surface-id <surface_id>

# 同时检查多种吸附物种在所有位点的覆盖情况
python hea_dataset.py check --root dataset --surface-id <surface_id> \
  -a N NH NH2 NH3

# 机器可读输出，便于脚本判定
python hea_dataset.py check --root dataset --surface-id <surface_id> \
  -a N --json
```

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `--surface-id ID` | 必填 | 待检查的组成 |
| `--root DIR` | `dataset` | 数据集根目录 |
| `-a`/`--adsorbates SPECIES...` | — | 附加检查的吸附物种，可指定多个 |
| `--max-list N` | `10` | 每类缺失最多列出的位点数（`0`=全部） |
| `--json` | — | 输出 JSON 而非人类可读表格 |

---

## 核心检查项

缺任意一项即视为不完整。每项都要求索引行与磁盘文件同时存在：

| 检查项 | 判定依据 |
|---|---|
| 初始结构已注册 | `surfaces.initial_cif` + `structures/00_initial_sqs.cif` |
| 弛豫结构已注册 | `surfaces.relaxed_cif` + `structures/01_relaxed_slab.cif` |
| 结构能量已入库 | `surfaces.total_energy_eV` 非空 |
| 表面原子已标注 | `top_atoms` 行 + `metadata/top_atoms.jsonl` + `atom_grid.npy` |
| 吸附位点已标注 | `fcc_sites` 行 + `metadata/fcc_sites.jsonl` + `site_grid.npy` |
| 哈密顿矩阵已入库 | `hamiltonian_exports` 行 + `openmx_slab/hamiltonian_d_surface.npz` |

每项的状态标记：

- `[OK]` — 索引与磁盘一致，产物存在。
- `[MISSING]` — 索引与磁盘均无该产物。
- `[PARTIAL]` — 索引与磁盘不一致（`present in index only` 表示只有索引行没有文件；`present in disk only` 表示只有文件未入库）。`[PARTIAL]` 也计为不完整。

> 结构能量是纯索引字段（无对应文件），因此只判断 `total_energy_eV` 是否非空。

---

## 吸附物种覆盖（可选）

`-a`/`--adsorbates` 指定的每个物种，会遍历该表面的**所有** `fcc_sites`，逐位点核对三项：

1. **已注册配置** — `adsorbate_configs` 中存在对应行。
2. **弛豫构型在盘** — `adsorbates/<物种>/<site_id>/01_relaxed_adsorbate.cif` 存在。
3. **吸附能已入库** — 该行 `adsorption_energy_eV` 非空。

支持一次传入多个物种（如 `-a N NH NH2 NH3`），逐个物种分别汇总。每个物种输出 `configs`/`structures`/`energies` 的 `已完成/总数`，并列出缺失的 `site_id`（数量受 `--max-list` 限制）。只有当所有位点三项齐备时，该物种才标记为 `[COMPLETE]`。

若表面尚未检测吸附位点（`fcc_sites` 为空），则物种覆盖直接判为不完整，并提示先运行 `detect-sites`。

---

## 输出示例

人类可读表格：

```text
============================================================
Surface : Co13Cr13Fe13Mn12Ni13
Status  : sites_detected
============================================================
[OK]      Initial structure registered
[OK]      Relaxed structure registered
[OK]      Slab total energy in index  (-123.456000 eV)
[OK]      Surface atoms indexed  (64 top atoms)
[OK]      Adsorption sites detected  (48 FCC sites)
[MISSING] Hamiltonian matrix stored  (not recorded)
------------------------------------------------------------
Adsorbate N   : configs 48/48  structures 47/48  energies 47/48   [INCOMPLETE]
           missing structure: site_0012
           missing energy   : site_0012
------------------------------------------------------------
Adsorbate NH  : configs 0/48  structures 0/48  energies 0/48   [INCOMPLETE]
           missing config   : site_0001, site_0002, ... (+46 more)
           ...
============================================================
Result  : INCOMPLETE (1 core item(s) missing/inconsistent; 2 adsorbate(s) incomplete)
```

`--json` 输出结构（供脚本消费）：

```json
{
  "surface_id": "...",
  "status": "sites_detected",
  "core_checks": [{"item": "...", "status": "ok|missing|partial", "detail": "..."}],
  "core_complete": false,
  "adsorbates": [{"adsorbate": "N", "total_sites": 48, "registered": 48,
                  "structures": 47, "energies": 47,
                  "missing_config": [], "missing_structure": ["site_0012"],
                  "missing_energy": ["site_0012"], "complete": false}],
  "adsorbates_complete": false,
  "complete": false
}
```

---

## 退出码

| 退出码 | 含义 |
|---|---|
| `0` | 所有被检查的项（核心 + 指定物种，如有）均完整 |
| `1` | 存在缺失或索引/磁盘不一致的项 |

在工作流脚本中可据此判断是否推进后续步骤，例如：

```bash
if python hea_dataset.py check --root dataset --surface-id "$sid" -a N; then
  echo "条目完整，可进入下一阶段"
fi
```
