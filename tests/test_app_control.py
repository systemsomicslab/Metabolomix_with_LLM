"""Smart App Control が未署名アセンブリを阻む状況を、Console 起動**前**に判定する。

実機の Smart App Control 状態に依存させない——判定は「ポリシー状態」と
「PE のバイト列」の 2 入力に対する純関数として書き、両方を注入して試験する。
背景は docs/HISTRY.md 2026-09-17(7)、kb
failures/msdial-console-smart-app-control-blocks-project-save.md。
"""
from __future__ import annotations

import struct

import pytest

from metabolomix.core import app_control


def _pe(*, signed: bool, pe32_plus: bool = True) -> bytes:
    """Certificate Table の有無だけが違う最小の PE を組む。"""
    pe_off = 0x80
    b = bytearray(0x400)
    b[0:2] = b"MZ"
    struct.pack_into("<I", b, 0x3C, pe_off)
    b[pe_off:pe_off + 4] = b"PE\x00\x00"
    struct.pack_into("<H", b, pe_off + 4, 0x8664 if pe32_plus else 0x014C)
    struct.pack_into("<H", b, pe_off + 0x18, 0x20B if pe32_plus else 0x10B)
    dd = pe_off + 0x18 + (0x70 if pe32_plus else 0x60)
    # データディレクトリ index 4 = Certificate Table (RVA, size)
    rva, size = (0x1000, 0x2000) if signed else (0, 0)
    struct.pack_into("<II", b, dd + 4 * 8, rva, size)
    return bytes(b)


def _write(tmp_path, name, data):
    p = tmp_path / name
    p.write_bytes(data)
    return p


# ---------- 署名判定 ----------

def test_a_pe_without_a_certificate_table_is_unsigned(tmp_path):
    assert app_control.has_authenticode_signature(_write(tmp_path, "u.dll", _pe(signed=False))) is False


def test_a_pe_with_a_certificate_table_is_signed(tmp_path):
    assert app_control.has_authenticode_signature(_write(tmp_path, "s.dll", _pe(signed=True))) is True


def test_a_pe32_binary_is_read_with_the_right_directory_offset(tmp_path):
    """PE32 と PE32+ でデータディレクトリの位置が違う。取り違えると誤判定する。"""
    assert app_control.has_authenticode_signature(
        _write(tmp_path, "s32.dll", _pe(signed=True, pe32_plus=False))) is True
    assert app_control.has_authenticode_signature(
        _write(tmp_path, "u32.dll", _pe(signed=False, pe32_plus=False))) is False


def test_an_unreadable_file_is_not_claimed_to_be_signed(tmp_path):
    """読めないものを「署名あり」と言わない（判定不能は未署名側へ倒さない）。"""
    assert app_control.has_authenticode_signature(tmp_path / "absent.dll") is None


# ---------- ブロック判定 ----------

def test_enforcing_policy_with_an_unsigned_assembly_is_blocked(tmp_path):
    dll = _write(tmp_path, "MsdialLcImMsApi.dll", _pe(signed=False))
    assert app_control.project_save_blocked(dll, policy_state=1) is True


def test_enforcing_policy_with_a_signed_assembly_is_allowed(tmp_path):
    """公式配布版は署名済み。ポリシーの有無だけで塞ぐと正当な構成を誤って止める。"""
    dll = _write(tmp_path, "MsdialLcImMsApi.dll", _pe(signed=True))
    assert app_control.project_save_blocked(dll, policy_state=1) is False


@pytest.mark.parametrize("policy_state", [0, 2, None])
def test_a_policy_that_is_not_enforcing_never_blocks(tmp_path, policy_state):
    """Off(0) / 評価(2) / 取得不能(None) では止めない。"""
    dll = _write(tmp_path, "MsdialLcImMsApi.dll", _pe(signed=False))
    assert app_control.project_save_blocked(dll, policy_state=policy_state) is False


def test_a_missing_assembly_does_not_claim_a_block(tmp_path):
    """判定材料が無いときに止めない——別の理由で落ちるならそちらで言うべき。"""
    assert app_control.project_save_blocked(tmp_path / "absent.dll", policy_state=1) is False
