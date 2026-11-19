"""Headless smoke test for a fresh ZFN WebSuite install.

Runs the full S1 pipeline on a synthetic sequence, checks the artefacts, then
builds fusion proteins/CDS from the result. Exits non-zero on the first
failure so it can be used as a build / CI gate.

    python packaging/smoke_test.py
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def main() -> int:
    from zfn_design import fusion_stage, pipeline
    from zfn_design.models import SequenceRecord

    seq = (
        "ATGGCCTTGACAGGTCTAGCAAGGTTACCGGATCCGTTAGCACCTGATCGTAGCTAGGCTA"
        "CCGTAGCTAGCTAGGCTTACCGATCGTAGCTAGCTAACGTTAGGCATCGTAGCTAGGCATC"
        "GTAGCTAGCTAGGCATCGTAGCTAGCTAGGCATCGTAGCTAGCTAGGCATCGTAGCTAGGC"
    )
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "run"
        payload = pipeline.run(
            [SequenceRecord(name="smoke", seq=seq)],
            spacers=(6,),
            fingers=3,
            n_designs=3,
            max_sites=20,
            fusion_keys=[],
            command="smoke_test",
        )
        paths = pipeline.write_outputs(payload, out)
        need = ("html", "json", "designs_tsv", "bundle")
        missing = [k for k in need if not paths.get(k) or not Path(paths[k]).exists()]
        if missing:
            print(f"FAIL: missing artefacts {missing}")
            return 1
        meta = json.loads(Path(paths["json"]).read_text(encoding="utf-8"))["meta"]
        print("stage:", meta.get("stage"), "kind:", meta.get("report_kind"))
        print("designs:", len(payload.get("designs") or []))
        if not payload.get("designs"):
            print("FAIL: no designs produced")
            return 1

        res = fusion_stage.run_fusion(out, top=1, arch_keys=["single_activator"])
        summary = res["summary"]
        print("fusion:", summary.get("n_entries"), "entries,", summary.get("n_ok"), "ok")
        if not summary.get("n_ok"):
            print("FAIL: fusion produced no usable CDS")
            return 1
        entry = res["entries"][0]
        chain = entry["chains"][0]
        cds = chain["cds"]
        if not cds.startswith("GCCACCATG") or cds[-3:] not in ("TGA", "TAA", "TAG"):
            print(f"FAIL: CDS framing wrong: {cds[:12]}...{cds[-6:]}")
            return 1
        print("cds:", len(cds), "bp, gc", round(chain["cds_gc"], 2))
    print("SMOKE OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
