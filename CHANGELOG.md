# 更新日志

版本号规则：`主.次.修订`。次版本号增加 = 新增一项对外能力（新的模式、新的输入、新的产物）。

## 0.8.1 — 2026-11-19

**英文界面全量覆盖；版本口径统一**

- Web 界面英文态此前残留中文（S4 面板、Target filter 提示、动态结果串、下拉选项等），
  改为「整句字典 + 子串表」两层翻译，覆盖 141 条子串与全部动态文案；
  四态（首页 / 设计后 / 融合后 / 脱靶后）经无头浏览器逐节点扫描，含中文行数为 0。
- 语言按钮在英文态显示 ZH（原为「中文」，不符合全英文要求）。
- 修复 MutationObserver 在语言切换期间漏翻新渲染节点的问题（下拉、引擎条、设计卡片）。
- 包版本、pyproject.toml、CITATION.cff、README、CHANGELOG 统一到 0.8.1；
  线上服务的 /api/info 同步回报 ersion=0.8.1、
_weights=7。
## 0.8.0 — 2026-11-19

**四阶段流程：DesignScore → Target filter → Off-target → Fusion（首次公开）**

- 界面拆成两个独立按钮：**Design Score（S1+S2）** 与 **Off-target Score（S3）**，各带进度条；
  Target filter（全基因组重复过滤）收进「可选」折叠区，默认关闭。
  后端新增 `/api/offtarget`，S3 只读上一次作业的 `report.json`，不重跑设计。
- **S2 全基因组重复过滤**（新模块 `genome_landing.py`）：候选按 *全基因组落点* 判，
  两条判据——落点压 RepeatMasker 注释区间，或 `site_seq` 全基因组拷贝数 ≥ 3 ——整条剔除。
  Aho-Corasick 单趟扫（GRCh38 约 94 s；装 pyahocorasick 更快，缺失自动退回 `bytes.find`）。
  **先过滤、后评分**：被剔除的候选只进报告，不参与打分。
- **S3 全基因组脱靶分**（`genome_offtarget.py` + `offtarget_stage.py`）：对 Design_score 前 N 个设计
  扫全基因组近似位点，按 0/1/2/3/4 错配加权（1.00/0.25/0.06/0.015/0.004）算 burden，
  `off-target score = max(0, 100 − 25 × burden)`；**与 Design_score 并列、不相加**。
  设计数触到命中上限时报告与表格都标注 burden 为下限。S3 结果页下面直接带 Fusion 区块与整包 zip。
- **S4 Fusion design**（`fusion.py` + `fusion_stage.py` + `codonopt.py`）：三架构
  （ZF-VP64 / ZF-KRAB / ZFN-ELD:KKK），N 端 SV40 NLS、10 aa 柔性 linker（GGGGSGGGGS）、
  成对核酸酶左右链分别用 FokI **ELD / KKK** 强制异二聚；蛋白以 Met 起，
  CDS 为 `GCCACC` Kozak + `ATG` … `TGA`，人源密码子优化（GC 40–60%、同聚物 ≤5、
  避开 EcoRI/BamHI/HindIII/XbaI/SalI/PstI/NdeI/XhoI/NotI/BsaI）。产物：
  `proteins.tsv` / `cds.tsv` / `cds.fasta` / `order_form.csv`，并重打 `zfn_bundle.zip`。
- **中英双语界面 + 内置使用说明书**：页面左上角切语言（记忆偏好）；右上角「使用说明书」打开
  `/manual?lang=zh|en`，页面顶部可下载 Word / Markdown。
- **报告带阶段标题**：`meta.stage` = S1/S2/S3/S4，报告小节按阶段动态编号（Target filter 与
  Off-target 只在真的跑了才出现）。
- **发布形态**：完整 B1H 指头库不再随包分发（申请制数据）；程序保留
  `zfn_design/data/finger_library/`，用户申请到数据后放进去即可，缺库时退回 3 条种子库并在
  `zfn-design info` 里明示。新增 `LICENSE`（MIT）、`CITATION.cff`、根级 `.gitignore`。
- **可复现包**：`packaging/Dockerfile`、`packaging/environment.yml`、`packaging/install.sh`
  （venv 失败自动回退 `pip --target`）、`packaging/smoke_test.py`（S1→S4 端到端冒烟）。
- **测试**：166 项（含 S2/S3/S4 与门禁用例）；没有完整指头库时 5 项依赖库的集成用例自动 skip。
  同时修掉旧套件里 8 项口径过期的断言（组件数 8→7、`n_sites` 语义、脱靶不再进设计分等）。

## 0.7.0 — 2026-11-15

**任意指数阵列（Custom：2–12 指）**

- 阵列长度不再只有"三指 / 六指"：`--fingers N` 支持 2–12，GUI 下拉框加「自定义」。
- 骨架按"三指一模块首尾相接"生成：N 端、两种模块内连接区、C 端全部逐字符取自 1AAY，
  模块之间用 `TGEKP` + 去掉 Met 的 N 端连接（与六指模板同一接法）。
  长度恒为 **28n + 6**：n=3 → 90 aa 正好还原 1AAY，n=6 → 184 aa 与手写六指模板逐字符相同。
- `designability.context_for_slot()` 改成按"模块内位置"映射（第 2 指→F2、第 3 指→F3、
  第 1 指→本位无数据），因此 4/5/7/8… 指都能自动指出哪些指位在 B1H 库里没有本位数据。
- 命令行对 `--fingers` 做 2–12 校验并给出可读错误；报告与 GUI 横幅按指数动态生成限制说明。
- 新增 `tests/test_custom_fingers.py`（15 项），测试总数 125 → 148。

## 0.6.0 — 2026-11-14

**全基因组脱靶评分**

- 新增 `--ot-designs N`：对分数最高的 N 个设计，在整条基因组上搜近似位点
  （正负两链、每侧错配上限可调），按错配加权求和得到基因组尺度的 `burden`，
  替换 `offtarget` 分项后重新打分，替换前的分数记在 `score.notes` 里。
- 新增 `zfn_design/offtarget_genome.py`：装了 numpy 走 18 位置向量化比较，
  没装就退回纯 Python 循环（慢约 40 倍）。
- 权重表扩展到 4 个错配：`1.00 / 0.25 / 0.06 / 0.015 / 0.004`
  （旧表只到 2，3/4 错配会被算成 0，等于白扫）。
- 新产物 `offtargets.tsv`（`design_id / chrom / start0 / strand / spacer / mm_* / weight`）。
- GUI「全基因组模式」面板加两个数字框：对几个设计做脱靶、每侧错配。
- 测试 117 → 125。

## 0.5.0 — 2026-11-14

**全基因组扫描 + 重复序列过滤**

- `scan` 支持 `--genome FASTA`：流式分块扫描（一次只驻留一条染色体），
  每块向右多读「窗口长度 − 1」个碱基保证跨块窗口完整。
- `--repeats`（RepeatMasker `.out` 或 `.bed`）：落在重复注释里的窗口
  **在构造候选之前**就跳过；`--regions BED` 先把范围缩到启动子/外显子。
- 同一段序列第二次出现标 `multi_copy`，不参与设计——兜住注释漏掉的拷贝。
- 新增 `zfn_design/genome.py`；`pipeline.run_genome()` 边扫边写 `sites.tsv`。
- GUI 左栏新增「全基因组模式」面板（填本机路径，不上传）。
- 测试 95 → 117。

## 0.4.0 — 2026-11-14

**单体型与组合型真正分开**

- 修掉一个真 bug：`--architecture` 之前只改了报告里的一个字符串，
  扫描/设计/打分/脱靶全是写死的成对逻辑，两种构型跑出完全一样的数。
- 现在：单体型（`single_activator` / `single_repressor`）只枚举一个半位点，
  没有 spacer、没有搭档、不产生双链断裂；打分少一个 spacer 分项，其余权重按比例重归一化；
  脱靶按「单个半位点像不像」判（`offtarget.search_single`）。
- 融合校验：单体型设计套用二聚化架构（`paired_nuclease`）直接抛 `FusionError`。
- `tests/test_single_mode.py` 8 项断言锁住这条区别。
- 测试 87 → 95。

## 0.3.0 — 2026-11-14

**六指阵列 + spacer 自定义**

- `--fingers 6`：半位点从 9 bp 变 18 bp，蛋白 184 aa（两个 1AAY 三指模块 + TGEKP 连接）。
- `--spacers` 支持 `6` / `5,6` / `5-7` / `5,7-8` 写法；GUI 里是「默认 5–7 / 自定义」两个单选。
- 明写限制：三指只读 9 bp，在人类基因组中不唯一；六指阵列的模块间界面**没有任何证据**，
  只是"结构上说得过去"的草稿。
- 测试 69 → 87。

## 0.2.0 — 2026-11-14

**接入真实 B1H 数据**

- 用户提供 Princeton `F2_union.zip` / `F3_union.zip` 后，随包库从 3 条扩到 **1,027 条，覆盖 64/64**。
- 导入器 `library import-b1h` 支持官方 6 残基连写格式；位置 4（不接触 DNA）统一填 `L`。
- 新增「本位上下文 / 上下文缺口」概念：库里缺某个位置的数据时退用同一三联体在别处的指头，
  并在报告与 `warnings` 里显式标为未验证。
- 该 csv **不可再分发**（申请制数据），发布包默认不含。

## 0.1.0 — 2026-11-14

**首个可用版本**

- 七阶段流水线：`ingest → targetscan → designability → build → score → report → fusion`。
- 零第三方依赖；`rules` 打分器 8 个分项，权重外置为 JSON、和恒为 1。
- 单文件 HTML 报告（零外部资源）、JSON / TSV 输出。
- 本地图形界面 `python -m zfn_design gui`（标准库 `http.server`）。
- 种子指头库 3 条（PDB 1AAY，可被任何人一条命令复核）。
