#!/usr/bin/env python3
"""打一个可分发的 release 包。

用法（在项目根目录跑）:

    python tools/make_release.py                 # 版本号取自 zfn_design/__init__.py
    python tools/make_release.py --no-verify     # 跳过 pytest
    python tools/make_release.py --with-b1h      # 把 B1H 指头库也打进去（仅内部使用）

产出:

    release/ZFN_Design_v<ver>/            解压即用的目录树
    release/ZFN_Design_v<ver>.zip         上面那个目录的压缩包
    release/ZFN_Design_v<ver>.zip.sha256  压缩包的校验和

设计原则:
- **只用标准库**（和主程序一致），任何装了 Python 3.10+ 的机器都能跑；
- 包里带上源码、测试、文档、启动脚本和 wheel —— 用户既能双击用，也能自己验证；
- 默认**不打** `finger_library_b1h.csv`：那份数据来自 Princeton 的申请制数据集，
  不能再分发（见 docs/数据来源与许可.md）。要打包给自己的服务器用才加 --with-b1h；
- 每次打包都写 MANIFEST.txt（逐文件 sha256）和 RELEASE_NOTES.md（本版说明），
  所以"这个 zip 是从哪来的、里面有什么"永远可查。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
import sys
import zipfile
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RELEASE_DIR = ROOT / "release"

# ---- 哪些进包、哪些不进 -------------------------------------------------

# 顶层就整个跳过的目录（用户结果、构建中间产物、缓存）
SKIP_DIRS = {
    "__pycache__", ".pytest_cache", ".git", ".idea", ".vscode",
    "build", "dist", "release", "unused_tmp",
    ".venv", "venv", "node_modules",
}
# 顶层以这些前缀开头的目录也不进包（out/ out3/ out6/ out6/ out/gui ...）
SKIP_PREFIXES = ("out",)

SKIP_SUFFIXES = {".pyc", ".pyo", ".log", ".tmp"}

# 这些即使存在也要显式排除（数据许可原因）
NEVER_SHIP_B1H = "finger_library_b1h.csv"
NEVER_SHIP = {NEVER_SHIP_B1H}
NEVER_SHIP_DIRS: set[str] = set()

# 必须有，缺了就直接报错
REQUIRED = [
    "README.md",
    "pyproject.toml",
    "zfn_design/__init__.py",
    "zfn_design/cli.py",
    "tools/build_docx.py",
    "启动工作台.bat",
]

TEXT_SUFFIXES = {".py", ".md", ".toml", ".json", ".csv", ".html", ".bat", ".txt", ".tsv", ".cfg", ".ini"}


def read_version() -> str:
    """从 `zfn_design/__init__.py` 读版本号，并和 pyproject.toml 对齐。"""
    init = (ROOT / "zfn_design" / "__init__.py").read_text(encoding="utf-8")
    m = re.search(r'__version__\s*=\s*"([^"]+)"', init)
    if not m:
        raise SystemExit("读不到 zfn_design/__init__.py 里的 __version__")
    ver = m.group(1)
    toml = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    m2 = re.search(r'^version\s*=\s*"([^"]+)"', toml, re.M)
    if not m2:
        raise SystemExit("读不到 pyproject.toml 里的 version")
    if m2.group(1) != ver:
        raise SystemExit(
            f"版本号不一致：zfn_design/__init__.py={ver}，pyproject.toml={m2.group(1)}"
        )
    return ver


def _skip_dir(name: str) -> bool:
    if name in SKIP_DIRS or name.endswith(".egg-info"):
        return True
    low = name.lower()
    return any(low == p or low.startswith(p) for p in SKIP_PREFIXES)


def collect_files(with_b1h: bool = False) -> list[Path]:
    """返回要打进包的文件（相对 ROOT 的 Path 列表，已排序）。"""
    picked: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(ROOT):
        here = Path(dirpath)
        rel_dir = here.relative_to(ROOT)
        if rel_dir == Path("."):
            dirnames[:] = [d for d in dirnames if not _skip_dir(d)]
        else:
            dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
        if any(part in NEVER_SHIP_DIRS for part in rel_dir.parts):
            dirnames[:] = []
            continue
        for fn in filenames:
            p = here / fn
            rel = p.relative_to(ROOT)
            if p.suffix.lower() in SKIP_SUFFIXES:
                continue
            if not with_b1h and fn in NEVER_SHIP:
                continue
            if p.stat().st_size == 0:
                continue
            picked.append(rel)
    return sorted(picked, key=lambda p: str(p).replace("\\", "/"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def human(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024.0
    return f"{n} B"


def strip_bom(path: Path) -> bool:
    """顺手清掉 UTF-8 BOM —— 带 BOM 的 pyproject.toml 会让 pytest 直接报错。"""
    if path.suffix.lower() not in TEXT_SUFFIXES:
        return False
    b = path.read_bytes()
    if b.startswith(b"\xef\xbb\xbf"):
        path.write_bytes(b[3:])
        return True
    return False


def _stage_for_wheel(with_b1h: bool) -> Path:
    """把打包需要的东西复制到一个临时目录，按需删掉申请制数据。

    **必须**这样子做：直接在仓库根目录构建 wheel，`package-data` 的 `data/*.csv`
    会把 ``finger_library_b1h.csv`` 一起打进 wheel —— 即使发布包把它们当普通文件
    排除了，wheel 仍会把它悄悄带出去（真实踩过这个坑）。
    """
    stage = Path(tempfile.mkdtemp(prefix="zfn_wheel_"))
    for name in ("pyproject.toml", "README.md"):
        src = ROOT / name
        if src.exists():
            shutil.copy2(src, stage / name)
    shutil.copytree(
        ROOT / "zfn_design", stage / "zfn_design",
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
    )
    if not with_b1h:
        (stage / "zfn_design" / "data" / NEVER_SHIP_B1H).unlink(missing_ok=True)
    return stage


def _wheel_members(wheel: Path) -> list[str]:
    with zipfile.ZipFile(wheel) as zf:
        return zf.namelist()


def build_wheel(ver: str, with_b1h: bool = False) -> tuple[Path | None, str]:
    """构建 wheel；返回 (路径, 说明)。没有 pip/网络时只警告，不失败。

    每次都重新构建：复用旧 wheel 会让"公开包里其实带着 B1H"这种错误悄无声息地留下来。
    """
    dist = ROOT / "dist"
    dist.mkdir(exist_ok=True)
    target = dist / f"zfn_design-{ver}-py3-none-any.whl"
    stage = _stage_for_wheel(with_b1h)
    outdir = Path(tempfile.mkdtemp(prefix="zfn_wheeldist_"))
    cmd = [sys.executable, "-m", "pip", "wheel", ".", "--no-deps", "-w", str(outdir), "--quiet"]
    try:
        proc = subprocess.run(cmd, cwd=stage, capture_output=True, text=True, timeout=600)
    except Exception as exc:  # noqa: BLE001 - 打包失败不该炸掉整个流程
        return None, f"wheel 构建异常：{exc}"
    finally:
        shutil.rmtree(stage, ignore_errors=True)
    built = sorted(outdir.glob(f"zfn_design-{ver}-*.whl"))
    if proc.returncode != 0 or not built:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-1:] or [""]
        return None, f"wheel 构建失败（{tail[0][:120]}）"
    wheel = built[0]
    has_b1h = any(m.endswith(NEVER_SHIP_B1H) for m in _wheel_members(wheel))
    if has_b1h != with_b1h:
        return None, ("wheel 内容与预期不符："
                      f"{'缺' if with_b1h else '混入了'} {NEVER_SHIP_B1H}")
    shutil.copy2(wheel, target)
    shutil.rmtree(outdir, ignore_errors=True)
    return target, "wheel 已生成（%s）" % ("含 B1H" if with_b1h else "不含 B1H")


def run_tests(timeout: int = 900) -> tuple[bool, str]:
    # 注意：pyproject 的 addopts 里已经有 -q，再叠一个 -q 会变成 -qq，
    # pytest 在 -qq 下**不打印** "N passed" 汇总行。所以这里用 -o addopts= 覆盖掉。
    cmd = [sys.executable, "-m", "pytest", "tests", "-o", "addopts=", "-q",
           "-p", "no:cacheprovider"]
    try:
        proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, timeout=timeout)
    except Exception as exc:  # noqa: BLE001
        return False, f"pytest 跑不起来：{exc}"
    out = (proc.stdout or "") + (proc.stderr or "")
    m = re.search(r"(\d+)\s+passed", out)
    n = m.group(1) if m else "?"
    if proc.returncode == 0:
        if n == "?":
            collected = re.search(r"collected\s+(\d+)\s+items?", out)
            n = collected.group(1) if collected else "?"
        return True, f"{n} passed"
    fails = re.findall(r"(\d+)\s+failed", out)
    return False, f"{n} passed / {fails[0] if fails else '?'} failed（见终端输出）"


def release_notes(ver: str, tests: tuple[bool, str], wheel: tuple[Path | None, str],
                  with_b1h: bool, n_files: int, n_doc_bytes: int, tag: str = "") -> str:
    today = datetime.now().strftime("%Y-%m-%d")
    title = f"ZFN_Design v{ver}" + (f"（{tag}）" if tag else "")
    banner = ""
    if with_b1h:
        banner = (
            "\n> ⚠️ **本包包含 Princeton B1H 指头库**（`zfn_design/data/finger_library_b1h.csv`，1,027 条）。\n"
            "> 这份数据来自 Persikov 等人 2015 年的申请制数据集，作者授权的是**非商业学术用途**。\n"
            "> 请在课题组内部使用，**不要上传到公开仓库、不要转发给组外**。对外分发请用不带该文件的版本。\n"
        )
    test_line = f"pytest **{tests[1]}**（在源码树、含 B1H 库的环境下）" if tests[0] else \
        f"⚠️ 测试未通过或未运行：{tests[1]}"
    if not with_b1h:
        test_line += "\n  　发布包内跑 `pytest tests` 会是 **145 passed / 3 skipped**——" \
                     "三条依赖 B1H 库的用例自动跳过（`-rs` 可看清单），这是预期行为，不是坏了。"
    wheel_line = f"`dist/{wheel[0].name}`" if wheel[0] else f"⚠️ 未生成 wheel（{wheel[1]}）"
    lib_line = (
        "**包含** `finger_library_b1h.csv`（1,027 条，Princeton B1H）——**仅限内部使用，不得转发**"
        if with_b1h else
        "只含 3 条 PDB 1AAY 种子指头（覆盖 2/64 三联体）。完整库需按 `docs/数据来源与许可.md` §2 自行申请导入"
    )
    b1h_license_line = (
        "**本包已附带**（`zfn_design/data/finger_library_b1h.csv`），"
        "仅限课题组内部使用；不要上传公开仓库、不要转发给组外。"
        if with_b1h else
        "随包**不附带**，需要者请按 `docs/数据来源与许可.md` §2 申请后，"
        "用 `zfn-design library import-b1h` 导入。"
    )
    return f"""# {title}

生成时间：{today}　|　文件数：{n_files}　|　文档合计：{human(n_doc_bytes)}
{banner}
一段 DNA → 候选 ZFN 位点 → 锌指阵列设计 → 可审计打分 → 单文件报告。
零第三方依赖（`numpy` 只是可选加速）。

## 怎么打开

1. **双击 `启动工作台.bat`**（Windows）：自动找 Python、起本地服务、开浏览器；
2. **命令行**：`python -m zfn_design gui`，或 `python -m zfn_design run --fasta examples/demo_target.fa --out out`；
3. **先看文档**：`docs/ZFN设计方案.docx`（38 页，含设计原理、打分明细、边界）；
4. **只想讲流程**：`docs/pipeline.html`（单文件，双击即看）。

## 本版状态

- {test_line}
- wheel：{wheel_line}
- 指头库：{lib_line}

## 必须知道的边界

- **分数不是活性预测**。它只表示"按公开设计规则看起来更像靠谱位点"，与实测亲和力/切割效率没有定量关系。
- **三指阵列只读 9 bp**，单独作为单体在人类基因组中并不唯一；要单拷贝特异性就做成对 ZFN，或加长阵列（`--fingers 4..12`，六指读 18 bp），并重新验证。
- **全基因组脱靶分是序列相似度**，不是切割概率；也不替代 BLAST / ZFN-Site。
- 所有输出都必须实验验证：T7E1 或扩增子深度测序 + in silico 脱靶复核 + 至少两次独立重复。

## 数据与许可

代码 MIT；文档 CC BY 4.0。
PDB 1AAY 结构数据属公有领域（CC0）。
Princeton B1H 数据集（Persikov 2015）为**申请制、非商业**用途，{b1h_license_line}
RepeatMasker 注释与基因组序列由使用者自备，仓库与发布包都不包含任何基因组文件。

## 校验

```bash
# Windows PowerShell
Get-FileHash .\\ZFN_Design_v{ver}.zip -Algorithm SHA256
# macOS / Linux
shasum -a 256 ZFN_Design_v{ver}.zip
```

对照 `ZFN_Design_v{ver}.zip.sha256`；包内 `MANIFEST.txt` 有逐文件校验和。
"""


def write_manifest(stage: Path, zippath: Path | None = None) -> int:
    """给 stage 目录里除 MANIFEST.txt 以外的每个文件写校验和。返回文件数。"""
    entries: list[tuple[Path, int, str]] = []
    for p in sorted(stage.rglob("*")):
        if not p.is_file() or p.name == "MANIFEST.txt":
            continue
        entries.append((p, p.stat().st_size, sha256_file(p)))
    total = sum(e[1] for e in entries)
    lines = [
        "# MANIFEST  ZFN_Design",
        f"# 生成时间: {datetime.now().isoformat(timespec='seconds')}",
        f"# 文件数: {len(entries)}（本文件不计）",
        "",
        f"{'sha256':64s}  {'bytes':>10s}  path",
    ]
    for p, size, digest in entries:
        lines.append(f"{digest:64s}  {size:10d}  {p.relative_to(stage).as_posix()}")
    lines.append("")
    lines.append(f"# 合计 {total} bytes（{human(total)}）")
    if zippath and zippath.exists():
        lines.append(f"# zip {zippath.name}: {sha256_file(zippath)}（{human(zippath.stat().st_size)}）")
    (stage / "MANIFEST.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return len(entries)


def make_zip(src_dir: Path, zip_path: Path) -> None:
    base = src_dir.name
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        for p in sorted(src_dir.rglob("*")):
            if p.is_file():
                zf.write(p, arcname=f"{base}/{p.relative_to(src_dir).as_posix()}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="make_release", description="把 ZFN_Design 打成一个可分发的 release 包",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--version", help="覆盖版本号（默认取 __init__.py）")
    ap.add_argument("--with-b1h", action="store_true",
                    help="把 finger_library_b1h.csv 也打进去（申请制数据，仅限内部使用）")
    ap.add_argument("--no-verify", action="store_true", help="跳过 pytest")
    ap.add_argument("--no-wheel", action="store_true", help="不构建/不附带 wheel")
    ap.add_argument("--tag", default=None,
                    help="文件名后缀，用来区分变体（默认：--with-b1h 时用 internal，否则不带）")
    ap.add_argument("--keep", action="store_true", help="保留同名目录（默认先删掉重建）")
    args = ap.parse_args(argv)

    for rel in REQUIRED:
        if not (ROOT / rel).exists():
            raise SystemExit(f"缺少必需文件：{rel}")

    ver = args.version or read_version()
    if not re.fullmatch(r"\d+\.\d+\.\d+", ver):
        raise SystemExit(f"版本号格式不对：{ver}")

    # 带 B1H 库的包名字里带 _internal，别人一眼就知道这不是能随便转发的版本
    tag = args.tag if args.tag is not None else ("internal" if args.with_b1h else "")
    tag = tag.strip("_- ")
    suffix = f"_{tag}" if tag else ""

    print(f"ZFN_Design release  v{ver}")
    print(f"  项目根目录: {ROOT}")

    # 1) 顺手清 BOM（打包前的最后一道防线；只看源码，不动 out/ 里的用户产物）
    bad = []
    for p in ROOT.rglob("*"):
        if not p.is_file() or p.suffix.lower() not in TEXT_SUFFIXES:
            continue
        if any(part in SKIP_DIRS or _skip_dir(part) for part in p.relative_to(ROOT).parts[:-1]):
            continue
        if strip_bom(p):
            bad.append(p)
    if bad:
        print(f"  清理 BOM: {len(bad)} 个文件")

    # 2) 测试
    tests = (True, "跳过（--no-verify）") if args.no_verify else run_tests()
    if args.no_verify:
        print("  测试: 跳过")
    else:
        print(f"  测试: {tests[1]}  {'✓' if tests[0] else '✗'}")
        if not tests[0]:
            raise SystemExit("测试没通过，先修好再打包（或加 --no-verify 强行打）")

    # 3) wheel
    wheel = (None, "跳过（--no-wheel）") if args.no_wheel else build_wheel(ver, args.with_b1h)
    print(f"  wheel: {wheel[1]}" + (f" → {wheel[0].name}" if wheel[0] else ""))

    # 4) 组装目录
    stage = RELEASE_DIR / f"ZFN_Design_v{ver}{suffix}"
    if stage.exists():
        if not args.keep:
            shutil.rmtree(stage)
        else:
            shutil.rmtree(stage)
    stage.mkdir(parents=True, exist_ok=True)

    files = collect_files(with_b1h=args.with_b1h)
    doc_bytes = sum((ROOT / f).stat().st_size for f in files if f.parts[0] == "docs")
    copied = 0
    for rel in files:
        dst = stage / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / rel, dst)
        copied += 1

    # 5) wheel 与说明
    if wheel[0]:
        (stage / "dist").mkdir(exist_ok=True)
        shutil.copy2(wheel[0], stage / "dist" / wheel[0].name)
    n_total = copied + 1 + (1 if wheel[0] else 0)   # + RELEASE_NOTES + wheel
    (stage / "RELEASE_NOTES.md").write_text(
        release_notes(ver, tests, wheel, args.with_b1h, n_total, doc_bytes, tag),
        encoding="utf-8")

    # 6) 压缩。注意：zip 无法包含自己的哈希，所以 MANIFEST 只记包内文件，
    #    整包的 sha256 写在外面的 .sha256 文件里。
    write_manifest(stage, None)
    zip_path = RELEASE_DIR / f"ZFN_Design_v{ver}{suffix}.zip"
    if zip_path.exists():
        zip_path.unlink()
    make_zip(stage, zip_path)
    digest = sha256_file(zip_path)
    (RELEASE_DIR / f"{zip_path.name}.sha256").write_text(
        f"{digest}  {zip_path.name}\n", encoding="utf-8", newline="\n")

    # 8) 报告
    print()
    print(f"  目录: {stage}")
    print(f"  压缩包: {zip_path}  ({human(zip_path.stat().st_size)})")
    print(f"  sha256: {digest}")
    print(f"  文件数: {copied + (2 if wheel[0] else 1)}")
    print()
    print("  提示: 解压后双击 `启动工作台.bat` 即可；也可以用包内 wheel 装进任何虚拟环境：")
    if wheel[0]:
        print(f"        pip install dist/{wheel[0].name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
