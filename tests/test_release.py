"""打包脚本（tools/make_release.py）的测试。

只测"该进包的进了、不该进包的没进、清单对得上"，不测 zip 能不能解压
（那是标准库的事）。
"""

from __future__ import annotations

import hashlib
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def _load():
    spec = importlib.util.spec_from_file_location("make_release", ROOT / "tools" / "make_release.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod


@pytest.fixture(scope="module")
def mr():
    return _load()


def test_version_is_readable_and_consistent(mr):
    ver = mr.read_version()
    assert ver.count(".") == 2
    assert all(part.isdigit() for part in ver.split("."))


def test_required_files_exist(mr):
    for rel in mr.REQUIRED:
        assert (ROOT / rel).exists(), f"缺少必需文件 {rel}"


def test_collect_files_excludes_junk(mr):
    files = [p.as_posix() for p in mr.collect_files()]
    assert any(f == "README.md" for f in files)
    assert any(f.startswith("zfn_design/") and f.endswith(".py") for f in files)
    assert any(f.startswith("tests/") for f in files)
    assert any(f.startswith("docs/") for f in files)
    # 不该进来
    assert not any("__pycache__" in f for f in files)
    assert not any(f.endswith(".pyc") for f in files)
    assert not any(f.startswith("out/") or f.startswith("out3/") for f in files)
    assert not any("egg-info" in f for f in files)
    assert not any(f.startswith("release/") for f in files)


def test_b1h_library_is_not_shipped_by_default(mr):
    """申请制数据，默认不进发布包（--with-b1h 才进）。"""
    default = [p.as_posix() for p in mr.collect_files(with_b1h=False)]
    assert not any("finger_library_b1h.csv" in f for f in default)
    # 种子库必须在
    assert any("finger_library_seed.csv" in f for f in default)


def test_manifest_lists_every_file_with_checksum(mr, tmp_path):
    stage = tmp_path / "pkg"
    (stage / "sub").mkdir(parents=True)
    (stage / "a.txt").write_text("hello", encoding="utf-8")
    (stage / "sub" / "b.bin").write_bytes(b"\x00\x01\x02")
    n = mr.write_manifest(stage, None)
    assert n == 2
    text = (stage / "MANIFEST.txt").read_text(encoding="utf-8")
    assert "文件数: 2" in text
    expect = hashlib.sha256(b"hello").hexdigest()
    assert expect in text
    assert "a.txt" in text and "sub/b.bin" in text
    # MANIFEST 不记录自己
    assert "MANIFEST.txt" not in text


def test_sha_matches_manifest_value(mr, tmp_path):
    f = tmp_path / "x.md"
    f.write_text("校验我", encoding="utf-8")
    assert mr.sha256_file(f) == hashlib.sha256(f.read_bytes()).hexdigest()


def test_bom_is_stripped(mr, tmp_path):
    f = tmp_path / "bom.py"
    f.write_bytes(b"\xef\xbb\xbfx = 1\n")
    assert mr.strip_bom(f) is True
    assert f.read_bytes() == b"x = 1\n"
    assert mr.strip_bom(f) is False  # 第二次没有 BOM 了


def test_skip_rules(mr):
    assert mr._skip_dir("out")
    assert mr._skip_dir("out3")
    assert mr._skip_dir("build")
    assert mr._skip_dir("zfn_design.egg-info")
    assert not mr._skip_dir("zfn_design")
    assert not mr._skip_dir("docs")
