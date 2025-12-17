"""
msa_osakt2_tool.py
==================

给 PP.py 的 "MSA 候选位点" 用的本地工具模块（不是 pip 库）。
核心能力：
1) 读取 Clustal/Clustal-Omega 的 .aln 对齐文件
2) 导出一个可读的 alignment view CSV（每列一个对齐位点）
3) 基于“多数派共识”给 OsAKT2（或自动识别的参考序列）生成候选突变位点表

设计目标：稳、可解释、出错也能给出清晰报错信息。
"""

from __future__ import annotations

import csv
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


# ----------------------------- Exceptions -----------------------------

class AlignmentParseError(RuntimeError):
    """对齐文件格式不符合预期（不是 Clustal / 内容损坏）"""


class ReferenceNotFoundError(RuntimeError):
    """找不到参考序列（例如 OsAKT2）"""


# ----------------------------- Data model -----------------------------

@dataclass
class MSA:
    names: List[str]                 # sequence ids in order
    seqs: List[str]                  # aligned sequences (same length)
    length: int                      # alignment columns

    def index_of(self, name: str) -> int:
        return self.names.index(name)

    def get(self, idx: int) -> str:
        return self.seqs[idx]


# ----------------------------- Helpers -----------------------------

# Clustal 里 gap 通常是 '-'; 有时也会出现 '.'（不同工具/格式的“空位/弱匹配”表现），稳妥起见当作 gap
_GAP_LIKE = set(["-", "."])

# 序列 chunk 允许的字符（含 gap / '.' / '*'）
_AA_RE = re.compile(r"^[A-Za-z\-\.\*]+$")


def _is_consensus_line(line: str) -> bool:
    """
    判定 Clustal 的共识行（consensus / match line）。

    特征：
    - 通常以空格开头（没有序列名）
    - 内容由 '*', ':', '.', 空格 组成（中间也会有空格）
    - 有时 split() 后会出现像 '::' 这种 token；我们用“去掉空格后字符集”判断最稳
    """
    s = line.rstrip("\n\r")
    if not s.strip():
        return False

    # 去掉所有空格/Tab 后，只剩匹配符号 -> 共识行
    compact = s.replace(" ", "").replace("\t", "")
    if compact and set(compact).issubset(set("*:.")):
        return True

    # 兼容极少数输出：共识行不一定以空格开头
    parts = s.split()
    if parts and len(parts) == 1 and set(parts[0]).issubset(set("*:.")):
        return True

    return False


def _parse_clustal(path: os.PathLike) -> MSA:
    """
    解析 .aln（CLUSTAL O / CLUSTAL W）格式。
    规则：
    - 跳过开头 'CLUSTAL' header
    - 块式拼接每条序列的对齐片段
    - 忽略共识行、空行
    - 忽略行尾的 position number
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"找不到对齐文件: {p}")

    lines = p.read_text(encoding="utf-8", errors="replace").splitlines()

    # 找 header
    header_ok = False
    for ln in lines[:10]:
        if ln.upper().startswith("CLUSTAL"):
            header_ok = True
            break
    if not header_ok:
        raise AlignmentParseError(
            f"文件看起来不是 CLUSTAL 格式（前几行找不到 'CLUSTAL'）：{p}"
        )

    order: List[str] = []
    chunks: Dict[str, List[str]] = {}

    for raw in lines:
        line = raw.rstrip("\n\r")
        if not line.strip():
            continue
        if line.upper().startswith("CLUSTAL"):
            continue
        if _is_consensus_line(line):
            continue

        # Clustal 的序列行：<name><spaces><chunk><spaces><optional number>
        parts = line.split()
        if len(parts) < 2:
            # 可能是奇怪的空白/注释，跳过
            continue
        name, chunk = parts[0], parts[1]

        # sanity: chunk 应该主要是 aa/gap
        if not _AA_RE.match(chunk):
            # 有些 clustal 输出 chunk 会带非字母字符，尽量容错：过滤掉非 AA/gap 字符
            chunk = re.sub(r"[^A-Za-z\-\.\*]", "", chunk)

        # 如果 chunk 仍为空，说明这行极可能是被误判的共识/杂行，直接跳过
        if not chunk:
            continue

        if name not in chunks:
            chunks[name] = []
            order.append(name)
        chunks[name].append(chunk)

    if not chunks:
        raise AlignmentParseError(f"未解析到任何序列行：{p}")

    names = order
    seqs = ["".join(chunks[n]) for n in names]
    length_set = {len(s) for s in seqs}
    if len(length_set) != 1:
        # 尝试指出最短/最长，方便定位问题
        lens = sorted((len(s), n) for n, s in zip(names, seqs))
        raise AlignmentParseError(
            "对齐后序列长度不一致，文件可能损坏或解析失败。"
            f" 最短: {lens[0]} 最长: {lens[-1]}"
        )
    length = length_set.pop()

    # 统一把 '.' 当 gap（更稳），避免后面统计出奇怪结果
    seqs = [s.replace(".", "-") for s in seqs]
    return MSA(names=names, seqs=seqs, length=length)


def _pick_reference(msa: MSA, prefer_patterns: Optional[Sequence[str]] = None) -> str:
    """
    自动挑参考序列。
    prefer_patterns: 按优先级匹配的关键词列表（大小写不敏感）。
    默认会优先尝试匹配 OsAKT2 / AKT2_ORY / OSAKT2 / 'AKT2' 等。
    """
    if prefer_patterns is None:
        prefer_patterns = [
            "OSAKT2", "OS-AKT2", "OS AKT2",  # 万一用户自己命名
            "AKT2_ORY", "AKT2-ORY", "AKT2ORY",  # 常见: AKT2_ORYSJ / AKT2_ORYSI...
            "AKT2",  # 兜底
        ]

    upper_names = [n.upper() for n in msa.names]

    # 逐个 pattern 找最先出现的
    for pat in prefer_patterns:
        up = pat.upper()
        for name_u, name in zip(upper_names, msa.names):
            if up in name_u:
                return name

    # 再兜底：如果只有 1 条，直接用
    if len(msa.names) == 1:
        return msa.names[0]

    # 最后兜底：用第一条（但这时你的候选位点更像“对第一条做建议”）
    return msa.names[0]


def _counts_to_string(counts: Dict[str, int], topn: int = 8) -> str:
    items = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    items = items[:topn]
    return ";".join(f"{k}:{v}" for k, v in items)


def _entropy_from_freqs(freqs: Sequence[float]) -> float:
    """Shannon entropy, base2."""
    ent = 0.0
    for f in freqs:
        if f > 0:
            ent -= f * math.log2(f)
    return ent


def _safe_mkdir_for_file(path: os.PathLike) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)


# ----------------------------- Public API -----------------------------

def export_alignment_view(
    aln_path: str,
    out_csv_path: str,
    reference_name: Optional[str] = None,
    max_sequences: Optional[int] = None,
) -> str:
    """
    把对齐展开成 CSV，方便肉眼查看 / Excel 过滤。

    CSV 列：
      aln_pos (1-based)
      ref_pos (1-based, ref gap -> empty)
      ref_aa
      consensus_aa (全体非gap的多数派)
      conservation (max_freq among non-gap, 0-1)
      gap_fraction (0-1)
      entropy (Shannon entropy on non-gap AAs)
      然后每条序列一列：<seq_name> = 该列字符

    返回：输出 CSV 路径（字符串）
    """
    msa = _parse_clustal(aln_path)

    if reference_name is None:
        reference_name = _pick_reference(msa)

    if reference_name not in msa.names:
        raise ReferenceNotFoundError(
            f"参考序列 '{reference_name}' 不在对齐里。可用序列名示例：{msa.names[:8]}"
        )

    ref_idx = msa.index_of(reference_name)

    names = msa.names
    seqs = msa.seqs

    if max_sequences is not None and max_sequences > 0 and len(names) > max_sequences:
        # 保留参考序列 + 前 max_sequences-1 条，避免 CSV 爆炸
        keep = [ref_idx] + [i for i in range(len(names)) if i != ref_idx][: max_sequences - 1]
        keep = sorted(set(keep))
        names = [names[i] for i in keep]
        seqs = [seqs[i] for i in keep]
        ref_idx = names.index(reference_name)
        msa = MSA(names=names, seqs=seqs, length=msa.length)

    ref_seq = msa.get(ref_idx)

    _safe_mkdir_for_file(out_csv_path)

    with open(out_csv_path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        header = [
            "aln_pos",
            "ref_pos",
            "ref_aa",
            "consensus_aa",
            "conservation",
            "gap_fraction",
            "entropy",
        ] + names
        w.writerow(header)

        ref_pos = 0
        for col in range(msa.length):
            col_chars = [s[col] for s in msa.seqs]
            gap_count = sum(1 for c in col_chars if c in _GAP_LIKE)
            gap_fraction = gap_count / len(col_chars)

            # 非 gap aa 统计
            nongap = [c.upper() for c in col_chars if c not in _GAP_LIKE]
            if nongap:
                aa_counts: Dict[str, int] = {}
                for c in nongap:
                    aa_counts[c] = aa_counts.get(c, 0) + 1
                total = len(nongap)
                consensus_aa, top = max(aa_counts.items(), key=lambda kv: kv[1])
                conservation = top / total
                freqs = [v / total for v in aa_counts.values()]
                entropy = _entropy_from_freqs(freqs)
            else:
                consensus_aa, conservation, entropy = "", 0.0, 0.0

            ref_char = ref_seq[col]
            if ref_char not in _GAP_LIKE:
                ref_pos += 1
                ref_pos_out = str(ref_pos)
            else:
                ref_pos_out = ""

            row = [
                str(col + 1),
                ref_pos_out,
                ref_char,
                consensus_aa,
                f"{conservation:.4f}",
                f"{gap_fraction:.4f}",
                f"{entropy:.4f}",
            ] + col_chars
            w.writerow(row)

    return out_csv_path


def suggest_candidates(
    aln_path: str,
    out_csv_path: str,
    reference_name: Optional[str] = None,
    min_non_gap_others: int = 5,
    max_gap_fraction_others: float = 0.2,
    min_conservation_others: float = 0.8,
    min_delta_conservation: float = 0.15,
    top_k: Optional[int] = None,
) -> str:
    """
    基于“其它序列的多数派共识”给参考序列生成候选点突变列表。

    逻辑（每列）：
    - 参考位点必须是非 gap
    - 其它序列（排除参考）非 gap 数 >= min_non_gap_others
    - 其它序列 gap_fraction <= max_gap_fraction_others
    - 其它序列的共识保守性 >= min_conservation_others
    - 参考 aa != 其它序列共识 aa
    - (可选) 与“全体”的保守性差值 >= min_delta_conservation：避免大家都很乱的位点

    输出 CSV 列：
      rank, ref_pos, ref_aa, suggested_aa, aln_pos,
      cons_others, conservation_others, n_non_gap_others, gap_fraction_others,
      conservation_all, gap_fraction_all, entropy_all,
      counts_others, note, reference_name

    返回：输出 CSV 路径（字符串）
    """
    msa = _parse_clustal(aln_path)

    if reference_name is None:
        reference_name = _pick_reference(msa)

    if reference_name not in msa.names:
        raise ReferenceNotFoundError(
            f"参考序列 '{reference_name}' 不在对齐里。可用序列名示例：{msa.names[:8]}"
        )

    ref_idx = msa.index_of(reference_name)
    ref_seq = msa.get(ref_idx)

    # others indices
    other_indices = [i for i in range(len(msa.names)) if i != ref_idx]
    if not other_indices:
        raise AlignmentParseError("对齐里只有参考序列一条，无法做“多数派共识”候选位点。")

    candidates: List[Tuple[float, Dict[str, str]]] = []

    ref_pos = 0
    for col in range(msa.length):
        ref_char = ref_seq[col]
        if ref_char in _GAP_LIKE:
            continue
        ref_pos += 1

        # 全体统计（用于辅助过滤/信息展示）
        col_chars_all = [s[col] for s in msa.seqs]
        gap_all = sum(1 for c in col_chars_all if c in _GAP_LIKE)
        gap_fraction_all = gap_all / len(col_chars_all)
        nongap_all = [c.upper() for c in col_chars_all if c not in _GAP_LIKE]
        if nongap_all:
            counts_all: Dict[str, int] = {}
            for c in nongap_all:
                counts_all[c] = counts_all.get(c, 0) + 1
            total_all = len(nongap_all)
            cons_all, top_all = max(counts_all.items(), key=lambda kv: kv[1])
            conservation_all = top_all / total_all
            freqs_all = [v / total_all for v in counts_all.values()]
            entropy_all = _entropy_from_freqs(freqs_all)
        else:
            cons_all, conservation_all, entropy_all = "", 0.0, 0.0

        # others 统计（真正用来做候选）
        col_chars_others = [msa.seqs[i][col] for i in other_indices]
        gap_others = sum(1 for c in col_chars_others if c in _GAP_LIKE)
        gap_fraction_others = gap_others / len(col_chars_others) if col_chars_others else 1.0
        nongap_others = [c.upper() for c in col_chars_others if c not in _GAP_LIKE]
        n_non_gap_others = len(nongap_others)

        if n_non_gap_others < min_non_gap_others:
            continue
        if gap_fraction_others > max_gap_fraction_others:
            continue

        counts_others: Dict[str, int] = {}
        for c in nongap_others:
            counts_others[c] = counts_others.get(c, 0) + 1
        total_others = n_non_gap_others
        cons_others, top_others = max(counts_others.items(), key=lambda kv: kv[1])
        conservation_others = top_others / total_others

        if conservation_others < min_conservation_others:
            continue
        if ref_char.upper() == cons_others:
            continue

        # 避免“全体也很乱”的位点（可调）
        if (conservation_others - conservation_all) < min_delta_conservation:
            # 如果全体也高度保守，还是可能是好位点 -> 放行
            if conservation_all < 0.95:
                continue

        # ranking：越保守越靠前，gap 越少越靠前
        score = (conservation_others * 1.2) - (gap_fraction_others * 0.5)

        note = (
            f"others_cons={cons_others}({conservation_others:.2f}), "
            f"ref={ref_char.upper()} differs; "
            f"gap_others={gap_fraction_others:.2f}"
        )
        row = {
            "ref_pos": str(ref_pos),
            "ref_aa": ref_char.upper(),
            "suggested_aa": cons_others,
            "aln_pos": str(col + 1),
            "cons_others": cons_others,
            "conservation_others": f"{conservation_others:.4f}",
            "n_non_gap_others": str(n_non_gap_others),
            "gap_fraction_others": f"{gap_fraction_others:.4f}",
            "conservation_all": f"{conservation_all:.4f}",
            "gap_fraction_all": f"{gap_fraction_all:.4f}",
            "entropy_all": f"{entropy_all:.4f}",
            "counts_others": _counts_to_string(counts_others),
            "note": note,
            "reference_name": reference_name,
        }
        candidates.append((score, row))

    # 排序 + 截断
    candidates.sort(key=lambda x: x[0], reverse=True)
    if top_k is not None and top_k > 0:
        candidates = candidates[:top_k]

    _safe_mkdir_for_file(out_csv_path)

    with open(out_csv_path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        header = [
            "rank",
            "ref_pos",
            "ref_aa",
            "suggested_aa",
            "aln_pos",
            "cons_others",
            "conservation_others",
            "n_non_gap_others",
            "gap_fraction_others",
            "conservation_all",
            "gap_fraction_all",
            "entropy_all",
            "counts_others",
            "note",
            "reference_name",
        ]
        w.writerow(header)

        for rank, (_, row) in enumerate(candidates, start=1):
            w.writerow([
                rank,
                row["ref_pos"],
                row["ref_aa"],
                row["suggested_aa"],
                row["aln_pos"],
                row["cons_others"],
                row["conservation_others"],
                row["n_non_gap_others"],
                row["gap_fraction_others"],
                row["conservation_all"],
                row["gap_fraction_all"],
                row["entropy_all"],
                row["counts_others"],
                row["note"],
                row["reference_name"],
            ])

    return out_csv_path


# ----------------------------- Optional CLI -----------------------------

def _main(argv: Optional[Sequence[str]] = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(
        description="OsAKT2 MSA helper: export alignment view and suggest candidate mutation sites."
    )
    ap.add_argument("aln", help="Clustal alignment file (.aln)")
    ap.add_argument("--ref", default=None, help="reference sequence name (default: auto)")
    ap.add_argument("--view_csv", default=None, help="output CSV for alignment view")
    ap.add_argument("--cand_csv", default=None, help="output CSV for candidate sites")
    ap.add_argument("--top", type=int, default=None, help="top K candidates")
    args = ap.parse_args(argv)

    aln = args.aln
    base = str(Path(aln).with_suffix(""))
    view_csv = args.view_csv or (base + "_alignment_view.csv")
    cand_csv = args.cand_csv or (base + "_candidate_sites.csv")

    export_alignment_view(aln, view_csv, reference_name=args.ref)
    suggest_candidates(aln, cand_csv, reference_name=args.ref, top_k=args.top)
    print("OK:", view_csv, cand_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
