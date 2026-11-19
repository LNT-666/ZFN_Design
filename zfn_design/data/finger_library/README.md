# 指头库放在这里 / Put the finger library here

本目录用于放置 **B1H 指头特异性全表**。该表不是本项目的代码，也**不随仓库分发**：它来自
Persikov et al. (2015) 发表的数据集，再分发受原始条款约束。请自行向数据作者申请（或从
论文补充材料按原始许可获取）后，把文件放进本目录。

This folder holds the **B1H zinc-finger specificity table**. The table is third-party data,
not part of this repository, and is **not redistributed here** because its original terms do
not permit us to relicense it. Obtain it from the original authors (or from the publication's
supplementary material under its own licence) and place it in this folder.

## 需要什么文件 / What to drop in

- 文件名必须叫 `finger_library_b1h.csv`（区分大小写）。
- 表头必须是这 10 列，顺序一致：

  ```
  helix_id,triplet,helix7,evidence_level,source_type,source_id,assay,affinity_nM,split_position,notes
  ```

- 一行一个指头变体：`triplet` 是它识别的三联体（如 `GCG`），`helix7` 是识别螺旋第 4–10 位
  七肽（如 `RSDELTR`），`evidence_level` 取 `measured` / `literature` / `derived`。

- File must be named `finger_library_b1h.csv` (case sensitive) with exactly this header:

  ```
  helix_id,triplet,helix7,evidence_level,source_type,source_id,assay,affinity_nM,split_position,notes
  ```

## 放好之后 / After you drop it in

```bash
zfn-design info          # library_name 应变成 b1h_persikov2015，records 约 1027，triplets 64/64
```

网页版同理：重启服务后首页引擎行会显示 `指头库 1027 条，覆盖 64/64 个三联体`。

Restart the server and the header line on the web UI should read 1027 records / 64 of 64
triplets.

## 没放会怎样 / Without it

程序不会崩。`zfn_design/library.py` 会退回只含 3 条 Zif268 指头（PDB 1AAY）的种子库：流程能
跑通、报告能生成，但候选位点覆盖极低，结果**不能用于真实设计**。

The program does not crash. `zfn_design/library.py` falls back to a three-finger seed library
derived from PDB 1AAY: the pipeline runs end to end, but coverage is negligible and the output
is **not usable for real design**.

## 为什么不能放仓库里 / Why it is not in the repository

这份表的上游许可没有授权第三方再分发；把 1027 条记录直接提交到公开仓库会侵犯原作者的数据
权利，也会让论文的 Availability 声明站不住。所以仓库只带种子库和这份说明。

The upstream licence does not grant redistribution rights to third parties. Committing the
1027-record table would violate the data authors' rights and would also invalidate the
Availability statement in our manuscript. The repository therefore ships only the seed library
and this note.

## 引用 / Citation

Persikov AV, Wetzel JL, Rowland EF, Oakes BL, Xu DJ, Singh M, Noyes MB. A systematic survey of
the Cys2His2 zinc finger DNA-binding landscape. *Nucleic Acids Res.* 2015;43:1965–1984.
