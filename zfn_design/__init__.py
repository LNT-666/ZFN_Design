"""ZFN_Design —— 从一段 DNA 出发，找 ZFN 位点、装锌指阵列、打分、出报告。

流水线（与设计方案文档的 7 个阶段一一对应）::

    ingest → targetscan → designability → build → score → report → fusion

设计哲学（三句话）
------------------
1. **不假装能预测活性。** 打分是把公开设计规则写成可审计的公式，
   不是亲和力模型；报告里每个分数都能被人工复核。
2. **证据等级不混用。** 每条指头自带 ``evidence_level``，
   实测数据与计算推断永远分开标注，缺口就报缺口，不凑数。
3. **每个阶段都能单独跑、单独测。** 阶段之间只通过数据模型通信。
"""

from __future__ import annotations

__version__ = "0.8.0"
__all__ = ["__version__"]
