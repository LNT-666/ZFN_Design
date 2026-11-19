# ZFN WebSuite (`zfn_design`)

**输入一段 DNA 序列或整个基因组，筛出能做 ZFN 的靶点、给出锌指蛋白序列、打分、做全基因组脱靶评估，最后给出可下单的融合蛋白与人源化 CDS。**

纯标准库 Python 3.10+，可离线运行；图形界面与命令行共用同一套 `pipeline`。

- 在线演示：<https://zfn.mobilome.cn/>（本组部署，可能随时停机）
- 文档：`docs/设计方案.md`、`docs/术语表.md`、`docs/数据来源与许可.md`；程序内自带中英文使用说明书（页面右上角「使用说明书」）

---

## ⚠️ 先读这三句

1. **分数不是活性预测。** `Design_score` 只表示"按公开设计规则看起来更像好位点"，与实测亲和力、切割效率没有定量关系。
2. **完整指头库不随仓库分发。** 本仓库只带 3 条 1AAY 种子库；完整库（Persikov 2015 B1H，1,027 条、覆盖 64/64 三联体）需**自行向原作者申请**后放入 `zfn_design/data/finger_library/`，见该目录下的 `README.md`。缺库时程序能跑，但结果不可用于真实设计。
3. **任何设计在合成前都必须湿实验验证**（T7E1 或扩增子深度测序）。程序只负责把候选缩到能做的数量。

---

## 四个阶段

| 阶段 | 做什么 | 入口 |
|---|---|---|
| S1 DesignScore | 枚举候选位点（9 + 5–7 + 9 bp，双链）→ 生成指头阵列 → 六项加权打分 → 排序 | GUI「Design Score」/ `zfn-design run` |
| S2 Target filter | 可选。按 RepeatMasker 注释 + 全基因组拷贝数 ≥ 3 两条判据剔除落在重复/低复杂度区的候选（**先过滤，后评分**） | 同一按钮，勾选 Target filter |
| S3 Off-target score | 对 Design_score 前 N 个设计在全基因组上找近似位点，算 burden → off-target score（与 Design_score 并列，不相加） | GUI「Off-target Score」/ `--offtarget-top` |
| S4 Fusion design | 选设计 → 拼 NLS-ZF-linker-效应域（VP64 / KRAB / FokI-ELD:KKK）→ 人源化 CDS（Kozak、ATG、TGA、避酶切位点）→ 订单表 | GUI「Fusion design」/ `zfn-design fusion` |

设计要点：
- 打分权重：evidence 0.428571、subsite 0.257143、gc_balance 0.1、finger_diversity 0.1、complexity 0.057143、self_complementarity 0.057143；spacer 仅作诊断（权重 0）。
- 全基因组扫描用 Aho-Corasick（`pyahocorasick`，可选；缺失自动退回 `bytes.find`），GRCh38 一遍约 94 s。
- CDS：`GCCACC` Kozak + `ATG` … `TGA`，GC 40–60%，无 ≥6 同聚物，避开 EcoRI/BamHI/HindIII/XbaI/SalI/PstI/NdeI/XhoI/NotI/BsaI 位点。

---

## 安装

```bash
# 1) 直接跑（零依赖，只用标准库）
python -m pip install -e .        # 或直接 python -m zfn_design ...

# 2) 推荐：带 numpy + pyahocorasick（全基因组扫描快很多）
bash packaging/install.sh

# 3) conda
conda env create -f packaging/environment.yml && conda activate zfn

# 4) Docker
docker build -f packaging/Dockerfile -t zfn-websuite .
docker run -p 8797:8797 -v "$PWD/out:/out" zfn-websuite
```

安装后自检：

```bash
zfn-design info          # 版本、权重、融合架构、元件证据状态
zfn-design library info  # 指头库覆盖情况（仓库版只有 3 条种子库）
python packaging/smoke_test.py   # 端到端冒烟：S1 → S4 → 校验 CDS
```

## 用法

```bash
# 图形界面（推荐）
zfn-design gui            # 默认 http://127.0.0.1:8797/

# 命令行：单靶点 + 可选重复过滤
zfn-design run --fasta examples/demo_target.fa --out out/ \
               --genome GRCh38.fa --repeats GRCh38.fa.out

# 命令行：对前 3 个设计做全基因组脱靶评估
zfn-design run --fasta gene.fa --out out/ --genome GRCh38.fa --offtarget-top 3

# 命令行：对已有作业做 Fusion（不重跑设计）
zfn-design fusion out/20261118-101010_600bp_ab12cd --top 3 \
               --arch single_activator paired_nuclease
```

GUI 三步：贴序列（或用示例）→ 点 **Design Score**（S1，可选勾 Target filter 走 S2）→ 想评估脱靶就点 **Off-target Score**（S3），想出蛋白/CDS 就点 **Fusion design**（S4）。每次作业落在 `out/gui/<时间>_<长度>_<运行号>/`，页面给下载按钮与整包 zip。

## 输出文件

| 文件 | 内容 |
|---|---|
| `report.html` / `report.json` | 人类可读报告 / 机器可读完整结果（含 meta.stage、参数、汇总） |
| `sites.tsv` | 候选位点（坐标、链、左右半位点、可设计性） |
| `designs.tsv` | 设计（指头阵列、蛋白、Design_score 分项） |
| `offtarget_genome.tsv` | S3：逐步设计的 burden、近位点分布、脱靶分 |
| `proteins.tsv` / `cds.tsv` / `cds.fasta` | S4：每条链的融合蛋白、CDS、Kozak/起始/终止 |
| `order_form.csv` | 下单表（链、架构、长度、GC、酶切位点自检、状态） |
| `zfn_bundle.zip` | 以上全部打包 |

## 指头库怎么补

见 `zfn_design/data/finger_library/README.md`：把申请到的表命名为 `finger_library_b1h.csv` 放进该目录即可，程序会自动优先使用；`zfn-design info` 应显示 `b1h_persikov2015`、约 1,027 条、64/64 三联体。

## 测试

```bash
python -m pytest -q
```

仓库版没有完整指头库时，依赖 64/64 覆盖的用例会自动 skip（`ZFN_LIBRARY_FULL=1` 可强制要求全库）。

## 许可与引用

- 代码：MIT，见 `LICENSE`。
- 指头权重数据：Persikov et al., *Nucleic Acids Res.* 2015;43:1965–1984 —— **不再分发**，请自行申请。
- 引用本软件：见 `CITATION.cff`。

## 目录

```text
zfn_design/        核心包（pipeline / targets / designability / score / genome / genome_landing /
                   genome_offtarget / fusion / fusion_stage / codonopt / report / gui / cli / manual）
zfn_design/data/   指头库（seed 3 条）、骨架模板、权重、融合元件、说明书
tests/             单元与端到端测试
packaging/         Dockerfile、environment.yml、install.sh、smoke_test.py
examples/          演示序列
tools/             发布与文档构建脚本
docs/              设计方案、术语表、数据来源与许可
```

## 联系

Issues 或邮件联系维护者（见 `CITATION.cff`）。
