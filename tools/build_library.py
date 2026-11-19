#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""把 Princeton B1H 数据集（Persikov 2015）编译成随包指头库。

用法::

    python tools/build_library.py --f2 <F2 目录> --f3 <F3 目录>
    python tools/build_library.py --f2 dump/F2_union --f3 dump/F3_union --top-n 8

输出 ``zfn_design/data/finger_library_b1h.csv``：PDB 种子（3 条，可再分发）
+ B1H 导入（每个三联体 × 每个位置取频率最高的 ``--top-n`` 条）。

数据来源与许可见 ``docs/数据来源与许可.md``：B1H 数据免费用于非商业
学术用途、需向作者申请，**不要把它再分发到公开仓库**。
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from zfn_design.library import FingerLibrary  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="编译随包 B1H 指头库")
    ap.add_argument("--f2", required=True, help="Princeton F2 union 数据集目录")
    ap.add_argument("--f3", required=True, help="Princeton F3 union 数据集目录")
    ap.add_argument("--top-n", type=int, default=8, help="每个三联体每个位置保留多少条")
    ap.add_argument("--pos4", default="L", help="第 4 位填充残基（数据未测）")
    ap.add_argument("--out", default=str(ROOT / "zfn_design" / "data" / "finger_library_b1h.csv"))
    args = ap.parse_args(argv)

    seed = FingerLibrary.load_csv(ROOT / "zfn_design" / "data" / "finger_library_seed.csv",
                                  name="seed_1AAY")
    f2 = FingerLibrary.import_b1h_dir(args.f2, split_position="F2", top_n=args.top_n,
                                      pos4=args.pos4, source_id="Persikov2015_B1H_F2")
    f3 = FingerLibrary.import_b1h_dir(args.f3, split_position="F3", top_n=args.top_n,
                                      pos4=args.pos4, source_id="Persikov2015_B1H_F3")

    merged = seed.merged(f2).merged(f3)
    merged.name = "b1h_persikov2015"
    dest = merged.to_csv(args.out)

    stats = merged.coverage_stats()
    print(f"写入 {dest}")
    for k in ("n_records", "n_curated_records", "triplets_covered", "triplets_total"):
        print(f"  {k:20s} {stats[k]}")
    # 位置覆盖度：F1 只有 PDB 种子里的 1 条（GCG）
    by_pos: dict[str, int] = {}
    for r in merged:
        by_pos[r.split_position] = by_pos.get(r.split_position, 0) + 1
    print(f"  {'by_split_position':20s} {by_pos}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
