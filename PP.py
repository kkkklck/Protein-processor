from typing import List, Dict, Tuple, Optional
from collections import defaultdict
import statistics
from pathlib import Path
import csv
import os, re, glob
import shutil
import subprocess
import math
import numpy as _np

OUTPUT_POLICIES = ("Minimal", "Standard", "Full")

try:
    import pandas as _pd
    import matplotlib.pyplot as _plt
except Exception:  # pragma: no cover - optional deps
    _pd = None
    _plt = None

try:
    from Bio.PDB import PDBParser as _PDBParser
except Exception:  # pragma: no cover - optional dep
    _PDBParser = None

try:
    from msa_osakt2_tool import export_alignment_view, suggest_candidates
except ImportError:  # pragma: no cover - optional dep
    export_alignment_view = None
    suggest_candidates = None

# 一字母 → 三字母氨基酸代码映射，只包含标准 20 个
ONE_TO_THREE = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}

# ===== HOLE / WSL 默认配置（只需要在自己电脑上改一次） =====
# 1. 在 WSL 里执行 `conda info --base` 得到 conda 的 base 路径，比如：
#    /home/k/miniforge3
# 2. 把 HOLE_WSL_CONDA_INIT 改成  <base>/etc/profile.d/conda.sh
# 3. 把 HOLE_WSL_CONDA_ENV 改成安装 hole 的那个环境名（例如 "hole_env"）

HOLE_WSL_CONDA_INIT = "$HOME/miniforge3/etc/profile.d/conda.sh"  # ← 根据实际路径改
HOLE_WSL_CONDA_ENV = "hole_env"  # ← 根据实际 env 名改
HOLE_WSL_EXE = "hole"  # env 里 HOLE 的命令名

# ===== Clustal-Omega / WSL 默认配置 =====
# lck这台机子现在的路径是 /usr/bin/clustalo
# 如果以后换机器，只需要改这里。
CLUSTALO_WSL_EXE = "/usr/bin/clustalo"


def _parse_clustal_alignment(aln_path: Path) -> Tuple[List[str], List[str]]:
    """解析 Clustal-Omega .aln，返回 (sequence_names, sequences)。"""

    aln_path = Path(aln_path)
    if not aln_path.is_file():
        raise FileNotFoundError(f"找不到 Clustal-Omega 对齐文件：{aln_path}")

    names: List[str] = []
    seq_chunks: Dict[str, List[str]] = defaultdict(list)

    with open(aln_path, encoding="utf-8", errors="ignore") as f:
        for raw_line in f:
            line = raw_line.rstrip("\n")
            striped = line.strip()
            if not striped or striped.upper().startswith("CLUSTAL"):
                continue
            # 共识行只含 *:. 等符号，直接跳过
            if striped.startswith(("*", ":", ".")):
                continue

            parts = line.split()
            if len(parts) < 2:
                continue

            name, seq_part = parts[0], parts[1]
            if name not in names:
                names.append(name)
            seq_chunks[name].append(seq_part)

    if not names:
        raise ValueError(f"在 {aln_path} 中没有解析到任何序列。")

    sequences = ["".join(seq_chunks[name]) for name in names]
    lengths = {len(seq) for seq in sequences}
    if len(lengths) != 1:
        raise ValueError(f"对齐序列长度不一致：{lengths}")

    return names, sequences


def _choose_reference(seq_names: List[str]) -> str:
    """选择参考序列：优先 OsAKT2，其次 AKT2_ORYS/ORYSJ，最后第一个。"""

    if not seq_names:
        raise ValueError("序列列表为空，无法选择参考序列。")

    preferred = ["OsAKT2"]
    fallbacks = ["AKT2_ORYS", "AKT2_ORYSJ"]

    for name in preferred + fallbacks:
        for candidate in seq_names:
            if candidate.upper().startswith(name.upper()):
                return candidate

    return seq_names[0]


def _export_alignment_view_fallback(aln_path: Path, view_csv_path: Path) -> None:
    names, sequences = _parse_clustal_alignment(aln_path)
    ref_name = _choose_reference(names)
    ref_idx = names.index(ref_name)

    alignment_len = len(sequences[0])
    ref_resnums: List[Optional[int]] = []
    ref_counter = 0
    for pos in range(alignment_len):
        ref_char = sequences[ref_idx][pos]
        if ref_char != "-":
            ref_counter += 1
            ref_resnums.append(ref_counter)
        else:
            ref_resnums.append(None)

    with open(view_csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Alignment_Pos", "Ref_Resnum", *names])

        for pos in range(alignment_len):
            ref_resnum = ref_resnums[pos] if ref_resnums[pos] is not None else ""
            row_chars = [seq[pos] for seq in sequences]
            writer.writerow([pos + 1, ref_resnum, *row_chars])


def _suggest_candidates_fallback(aln_path: Path, cand_csv_path: Path) -> None:
    names, sequences = _parse_clustal_alignment(aln_path)
    ref_name = _choose_reference(names)
    ref_idx = names.index(ref_name)

    alignment_len = len(sequences[0])
    ref_resnums: List[Optional[int]] = []
    ref_counter = 0
    for pos in range(alignment_len):
        ref_char = sequences[ref_idx][pos]
        if ref_char != "-":
            ref_counter += 1
            ref_resnums.append(ref_counter)
        else:
            ref_resnums.append(None)

    with open(cand_csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "Alignment_Pos",
                "Ref_Resnum",
                "Ref_AA",
                "Conserved_AA",
                "Support_Count",
                "Total_Others",
            ]
        )

        for pos in range(alignment_len):
            ref_char = sequences[ref_idx][pos]
            if ref_char == "-":
                continue

            others = [
                sequences[i][pos]
                for i in range(len(names))
                if i != ref_idx and sequences[i][pos] != "-"
            ]

            if len(others) < 2:
                continue

            conserved_set = set(others)
            if len(conserved_set) != 1:
                continue

            conserved = others[0]
            if conserved == ref_char:
                continue

            ref_resnum = ref_resnums[pos] if ref_resnums[pos] is not None else ""
            writer.writerow(
                [
                    pos + 1,
                    ref_resnum,
                    ref_char,
                    conserved,
                    len(others),
                    len(others),
                ]
            )


def _tables_dir(out_dir: str) -> str:
    path = os.path.join(out_dir, "tables")
    os.makedirs(path, exist_ok=True)
    return path


def table_path(out_dir: str, filename: str) -> str:
    """Return a table path under tables/ (creating the directory)."""

    return os.path.join(_tables_dir(out_dir), filename)


def find_table(out_dir: str, filename: str) -> str:
    """Find a CSV either under tables/ or directly in out_dir."""

    tables_candidate = os.path.join(out_dir, "tables", filename)
    if os.path.exists(tables_candidate):
        return tables_candidate
    legacy_candidate = os.path.join(out_dir, filename)
    return legacy_candidate


def _model_dir(out_dir: str, label: str) -> str:
    path = os.path.join(out_dir, label)
    os.makedirs(path, exist_ok=True)
    return path


def _is_better(src: str, dst: str) -> bool:
    """Heuristic: prefer newer or larger files."""

    try:
        src_stat = os.stat(src)
        dst_stat = os.stat(dst)
    except OSError:
        return True

    if src_stat.st_mtime != dst_stat.st_mtime:
        return src_stat.st_mtime > dst_stat.st_mtime
    return src_stat.st_size >= dst_stat.st_size


def organize_outputs(out_dir: str, labels: List[str], policy: str = "Standard") -> None:
    """Move legacy scattered files into per-model folders and deduplicate."""

    out_dir = (out_dir or "").strip()
    if not out_dir:
        return

    trash = os.path.join(out_dir, "_trash")
    os.makedirs(trash, exist_ok=True)

    for name in os.listdir(out_dir):
        src = os.path.join(out_dir, name)
        if not os.path.isfile(src):
            continue

        for label in labels:
            prefix = f"{label}_"
            if not name.startswith(prefix):
                continue

            dst_dir = _model_dir(out_dir, label)
            dst = os.path.join(dst_dir, name)

            if os.path.exists(dst):
                if _is_better(src, dst):
                    os.replace(src, dst)
                else:
                    os.replace(src, os.path.join(trash, name))
            else:
                os.replace(src, dst)
            break


def cleanup_minimal(out_dir: str, labels: List[str]) -> None:
    """Archive non-essential artifacts when Minimal output policy is selected."""

    archive = os.path.join(out_dir, "_archive")
    os.makedirs(archive, exist_ok=True)

    keep_dirs = {"gate_sites", "tables"}
    for entry in os.listdir(out_dir):
        if entry in keep_dirs or entry.startswith("."):
            continue
        path = os.path.join(out_dir, entry)
        if os.path.isdir(path) and entry in labels:
            shutil.move(path, os.path.join(archive, entry))
        elif os.path.isdir(path) and entry in {"_trash", "_archive"}:
            continue
        elif os.path.isdir(path):
            shutil.move(path, os.path.join(archive, entry))
        elif os.path.isfile(path):
            shutil.move(path, os.path.join(archive, entry))

# ===========================
# HOLE 管道相关工具函数
# ===========================


def hole_win_to_wsl(path: str) -> str:
    """把 Windows 风格路径转换为 WSL 路径。"""
    if path is None:
        raise ValueError("path 不可为空。")
    cleaned = path.strip().strip('"').strip("'")
    if not cleaned:
        raise ValueError("path 为空，无法转换。")
    if cleaned.startswith("/"):
        return cleaned
    if len(cleaned) >= 2 and cleaned[1] == ":":
        drive = cleaned[0].lower()
        tail = cleaned[2:].lstrip("\\/")
        tail = tail.replace("\\", "/")
        return f"/mnt/{drive}/{tail}"
    raise ValueError(f"无法识别的 Windows 路径：{path}")


def hole_write_input(
    base_dir_win: str,
    model: str,
    cpoint: Tuple[float, float, float],
    cvect: Tuple[float, float, float],
    sample: float = 0.25,
    endrad: float = 15.0,
    radius_filename: str = "simple.rad",
) -> str:
    """根据配置生成 HOLE 所需的 .inp 文件。"""

    base_dir_win = (base_dir_win or "").rstrip("\\/")
    if not base_dir_win:
        raise ValueError("base_dir_win 不能为空")
    base_dir_wsl = hole_win_to_wsl(base_dir_win)

    coord_wsl = f"{base_dir_wsl}/{model}.pdb"
    radius_wsl = f"{base_dir_wsl}/{radius_filename}"
    sphpdb_wsl = f"{base_dir_wsl}/{model}_hole_spheres.pdb"
    mappdb_wsl = f"{base_dir_wsl}/{model}_hole_profile.pdb"

    cpx, cpy, cpz = cpoint
    cvx, cvy, cvz = cvect

    lines = [
        f"coord  {coord_wsl}",
        f"radius {radius_wsl}",
        f"cpoint {cpx:.2f} {cpy:.2f} {cpz:.2f}",
        f"cvect  {cvx:.2f} {cvy:.2f} {cvz:.2f}",
        f"sample {sample:.2f}",
        f"endrad {endrad:.2f}",
        "ignore HOH",
        f"sphpdb {sphpdb_wsl}",
        f"mappdb {mappdb_wsl}",
        "",
    ]

    out_path = base_dir_win + f"\\{model}_hole.inp"
    with open(out_path, "w", encoding="ascii", errors="ignore") as f:
        f.write("\n".join(lines))
    return out_path


def hole_run_in_wsl(base_dir_win: str, model: str, hole_cmd: str = "") -> None:
    """
    在 WSL 环境中运行 HOLE。

    参数 hole_cmd：
      - 为空字符串 / "hole" / "auto"：走自动模式：
            . HOLE_WSL_CONDA_INIT
            conda activate HOLE_WSL_CONDA_ENV
            HOLE_WSL_EXE < inp > log
      - 其他非空：认为是高级用户手动指定的命令，原样执行。
    """

    base_dir_win = (base_dir_win or "").rstrip("\\/")
    if not base_dir_win:
        raise ValueError("base_dir_win 不能为空")
    base_dir_wsl = hole_win_to_wsl(base_dir_win)

    inp_name = f"{model}_hole.inp"
    log_name = f"{model}_hole.log"

    # 决定在 WSL 里真正要执行的 HOLE 命令
    cmd_raw = (hole_cmd or "").strip()
    if not cmd_raw or cmd_raw.lower() in {"hole", "auto"}:
        # —— 自动模式：使用上面配置的 conda 环境 ——
        # 注意：假设 HOLE_WSL_CONDA_INIT / HOLE_WSL_CONDA_ENV / HOLE_WSL_EXE
        # 已经被你改成正确值。
        cmd_in_shell = (
            f'. {HOLE_WSL_CONDA_INIT} && '
            f'conda activate {HOLE_WSL_CONDA_ENV} && '
            f'{HOLE_WSL_EXE}'
        )
    else:
        # —— 高级模式：直接使用用户输入的一整串命令 ——
        cmd_in_shell = cmd_raw

    inner_cmd = f'cd "{base_dir_wsl}" && {cmd_in_shell} < {inp_name} > {log_name}'
    full_cmd = ["wsl", "bash", "-lc", inner_cmd]
    subprocess.run(full_cmd, check=True)


def hole_parse_profile(log_path: str) -> List[Tuple[float, float]]:
    """从 HOLE 日志中提取 profile 数据。"""

    profile: List[Tuple[float, float]] = []
    with open(log_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if "mid-point" in line or "(sampled)" in line:
                parts = line.split()
                if len(parts) < 2:
                    continue
                try:
                    s = float(parts[0])
                    r = float(parts[1])
                except ValueError:
                    continue
                profile.append((s, r))
    if not profile:
        raise ValueError(f"日志中没有 profile 信息：{log_path}")
    return profile


def compute_gate_metrics(
    profile: List[Tuple[float, float]], threshold: float = 1.4
) -> Tuple[float, float, float, float, float]:
    """计算 profile 的最小半径以及 gate 区间。"""

    if not profile:
        raise ValueError("profile 不能为空。")

    # 先找出最小半径以及对应的位置
    s_min, r_min = min(profile, key=lambda item: item[1])

    # 为了避免日志里偶尔顺序被打乱，先按 s 排序
    sorted_profile = sorted(profile, key=lambda item: item[0])

    gate_segments: List[Tuple[float, float]] = []
    gate_start = gate_end = None
    for s, r in sorted_profile:
        if r < threshold:
            if gate_start is None:
                gate_start = s
            gate_end = s
        else:
            if gate_start is not None and gate_end is not None:
                gate_segments.append((gate_start, gate_end))
            gate_start = gate_end = None
    if gate_start is not None and gate_end is not None:
        gate_segments.append((gate_start, gate_end))

    if gate_segments:
        gate_start, gate_end = max(gate_segments, key=lambda seg: seg[1] - seg[0])
        gate_length = gate_end - gate_start
    else:
        gate_start = gate_end = float("nan")
        gate_length = float("nan")

    return r_min, s_min, gate_start, gate_end, gate_length


def hole_summarize_logs(log_paths: Dict[str, str], out_dir: str) -> None:
    """根据多个 HOLE 日志生成 CSV 汇总。"""

    if _pd is None:
        raise RuntimeError("需要安装 pandas 才能生成 CSV 汇总。")

    profiles_rows: List[Dict[str, float]] = []
    min_table_rows: List[Dict[str, float]] = []
    min_summary_rows: List[Dict[str, str]] = []

    for model, path in log_paths.items():
        prof = hole_parse_profile(path)
        for s, r in prof:
            profiles_rows.append({"Model": model, "s_A": s, "r_A": r})

        r_min, s_min, gate_start, gate_end, gate_length = compute_gate_metrics(prof)
        min_table_rows.append(
            {
                "Model": model,
                "s_min_A": s_min,
                "r_min_A": r_min,
                "gate_start_A": gate_start,
                "gate_end_A": gate_end,
                "gate_length_A": gate_length,
            }
        )
        min_summary_rows.append({
            "Model": model,
            "Summary": f"Minimum radius found: {r_min:8.3f} Å at s = {s_min:8.3f} Å",
        })

    os.makedirs(out_dir, exist_ok=True)
    _pd.DataFrame(profiles_rows).to_csv(os.path.join(out_dir, "hole_profile_samples.csv"), index=False)
    _pd.DataFrame(min_table_rows).to_csv(os.path.join(out_dir, "hole_min_table.csv"), index=False)
    _pd.DataFrame(min_summary_rows).to_csv(os.path.join(out_dir, "hole_min_summary.csv"), index=False)


def merge_all_metrics(hole_dir: str, sasa_dir: str, out_csv: str) -> None:
    """把 HOLE、SASA 和氢键的统计合成一张总表。"""

    if _pd is None:
        raise RuntimeError("需要安装 pandas 才能合并 CSV。")

    hole_csv = find_table(hole_dir, "hole_min_table.csv")
    sasa_csv = find_table(sasa_dir, "sasa_hbonds_summary.csv")

    if not os.path.exists(hole_csv):
        raise FileNotFoundError(f"未找到 HOLE 结果：{hole_csv}")
    if not os.path.exists(sasa_csv):
        raise FileNotFoundError(f"未找到 SASA/H-bond 结果：{sasa_csv}")

    df_hole = _pd.read_csv(hole_csv)
    df_sasa = _pd.read_csv(sasa_csv)
    df = df_hole.merge(df_sasa, on="Model", how="outer")
    os.makedirs(os.path.dirname(out_csv) or ".", exist_ok=True)
    df.to_csv(out_csv, index=False)

def _score_class(total_score: float) -> str:
    if total_score is None or (isinstance(total_score, float) and math.isnan(total_score)):
        return "unknown"
    if total_score >= 1.0:
        return "better_than_WT"
    if total_score >= 0.0:
        return "similar_to_WT"
    return "worse_than_WT"

def _fit_weights(df: _pd.DataFrame, standards_csv: str) -> Tuple[float, float, float]:
    standards = _pd.read_csv(standards_csv)
    if "Model" not in standards.columns or "y" not in standards.columns:
        raise ValueError("standards.csv 需要包含 'Model' 和 'y' 两列。")

    standards["y"] = _pd.to_numeric(standards["y"], errors="coerce")

    merged = df.merge(standards[["Model", "y"]], on="Model", how="inner")
    if merged.empty:
        raise ValueError("standards.csv 与 metrics_all.csv 没有重叠的 Model，无法拟合。")

    merged = merged.dropna(subset=["GateTightScore", "HydroScore", "y"])
    if len(merged) < 3:
        raise ValueError("用于拟合的标准点少于 3 个或存在空值，无法拟合。")

    X = _np.c_[
        merged["GateTightScore"].to_numpy(),
        merged["HydroScore"].to_numpy(),
        _np.ones(len(merged)),
    ]
    y = merged["y"].to_numpy()

    w1, w2, b = _np.linalg.lstsq(X, y, rcond=None)[0]
    scale = abs(w1) + abs(w2)
    if scale:
        w1, w2 = w1 / scale, w2 / scale
    return float(w1), float(w2), float(b)


def score_metrics_file(
    csv_path: str,
    wt_name: str = "WT",
    pdb_dir: str | None = None,
    standards_csv: str | None = None,
    default_weights: Tuple[float, float, float] = (0.7, 0.3, 0.0),
) -> str:
    """
    读 metrics_all.csv，按 WT 计算 GateTightScore / TotalScore / ScoreClass，
    若提供 standards_csv 则对 GateTightScore/HydroScore 做线性拟合得到权重；
    若提供 pdb_dir 且安装了 Biopython，则顺带拼上 pLDDT 置信度。
    返回最终写出的 metrics_scored.csv 路径。
    """

    if _pd is None:
        raise RuntimeError("需要 pandas 才能打分。")

    csv_path = os.path.abspath(csv_path)
    if not os.path.isfile(csv_path):
        raise FileNotFoundError(f"找不到 metrics_all.csv：{csv_path}")

    df = _pd.read_csv(csv_path)

    required_cols = ["Model", "r_min_A", "gate_length_A", "HBonds", "SASA_residue"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"metrics_all.csv 缺少这些列：{missing}")

    if wt_name not in df["Model"].values:
        raise ValueError(f"在 Model 列里找不到 WT='{wt_name}'，请确认命名。")

    wt = df[df["Model"] == wt_name].iloc[0]
    r0 = wt["r_min_A"]
    L0 = wt["gate_length_A"]
    s0 = wt["SASA_residue"]
    h0 = wt["HBonds"] if wt["HBonds"] != 0 else 1.0

    def _rel_delta(base: float, value: float) -> float:
        return (value - base) / base if base not in (0, None) else 0.0

    df["d_rmin_vs_WT"] = df["r_min_A"] - r0
    df["d_gateL_vs_WT"] = df["gate_length_A"] - L0
    df["d_SASAres_vs_WT"] = df["SASA_residue"] - s0
    df["d_HBonds_vs_WT"] = df["HBonds"] - h0

    def _gate_score(row: _pd.Series) -> float:
        delta_r = _rel_delta(r0, row["r_min_A"])
        delta_L = _rel_delta(L0, row["gate_length_A"])
        return (-delta_r) + 0.5 * delta_L

    def _hydro_score(row: _pd.Series) -> float:
        delta_s = _rel_delta(s0, row["SASA_residue"])
        delta_h = _rel_delta(h0, row["HBonds"])
        return (-0.5 * delta_s) + 0.5 * delta_h

    df["GateTightScore"] = df.apply(_gate_score, axis=1)
    df["HydroScore"] = df.apply(_hydro_score, axis=1)

    if standards_csv:
        w_gate, w_hydro, bias = _fit_weights(df, standards_csv)
    else:
        w_gate, w_hydro, bias = default_weights

    df["TotalScore"] = df["GateTightScore"] * w_gate + df["HydroScore"] * w_hydro + bias
    wt_total_series = df.loc[df["Model"] == wt_name, "TotalScore"]
    if not wt_total_series.empty:
        wt_total = wt_total_series.iloc[0]
        if not _pd.isna(wt_total):
            df["TotalScore"] = df["TotalScore"] - wt_total
    df["ScoreClass"] = df["TotalScore"].apply(_score_class)

    if pdb_dir and _PDBParser is not None:
        conf_df = plddt_summary_for_models(pdb_dir, df["Model"].astype(str).tolist())
        df = df.merge(conf_df, on="Model", how="left")

    out_path = os.path.join(os.path.dirname(csv_path), "metrics_scored.csv")
    df.to_csv(out_path, index=False, encoding="utf-8-sig")
    return out_path

def make_stage3_table(out_dir: str, pick_models: List[str] | None = None) -> str:
    """
    生成 Stage3 决策表模板：
    - 读取 metrics_all.csv，并合并 metrics_scored.csv 中的新增列
    - 追加两列定性占位符：Patch_Electrostatics、Contacts_Qualitative
    - 若 gate_sites 图已生成，则填入对应的图片路径
    返回写出的 stage3_table.csv 路径。
    """

    if _pd is None:
        raise RuntimeError("需要 pandas 才能生成 Stage3 表格。")

    out_dir = (out_dir or "").strip()
    if not out_dir:
        raise ValueError("out_dir 不能为空")

    metrics_all = find_table(out_dir, "metrics_all.csv")
    if not os.path.isfile(metrics_all):
        raise FileNotFoundError(f"找不到 metrics_all.csv：{metrics_all}")

    df = _pd.read_csv(metrics_all)

    if pick_models:
        pick_models = [m for m in pick_models if m in df["Model"].astype(str).tolist()]

    scored_path = find_table(out_dir, "metrics_scored.csv")
    if os.path.isfile(scored_path):
        scored_df = _pd.read_csv(scored_path)
        if "Model" in scored_df.columns:
            new_cols = [c for c in scored_df.columns if c != "Model" and c not in df.columns]
            if new_cols:
                df = df.merge(scored_df[["Model", *new_cols]], on="Model", how="left")

    if pick_models:
        df = df[df["Model"].isin(pick_models)]
        if pick_models:
            df["Model"] = _pd.Categorical(df["Model"], pick_models)
            df = df.sort_values("Model")
            df["Model"] = df["Model"].astype(str)

    sites_dir = os.path.join(out_dir, "gate_sites")

    def _img_path(label: str, suffix: str) -> str:
        path = os.path.join(sites_dir, f"{label}_{suffix}.png")
        return path if os.path.exists(path) else ""

    df["Patch_Electrostatics"] = ""
    df["Contacts_Qualitative"] = ""
    df["sites_contacts_img"] = df["Model"].apply(lambda m: _img_path(str(m), "sites_contacts"))
    df["sites_coulombic_img"] = df["Model"].apply(lambda m: _img_path(str(m), "sites_coulombic"))

    os.makedirs(_tables_dir(out_dir), exist_ok=True)
    out_path = table_path(out_dir, "stage3_table.csv")
    df.to_csv(out_path, index=False, encoding="utf-8-sig")
    return out_path

def run_osakt2_msa(
        fasta_path: str,
        clustalo_cmd: str = "clustalo",
) -> Tuple[str, str, str]:
    """
    给定一个多序列 FASTA，调用 Clustal-Omega 生成 .aln，
    再用 msa_osakt2_tool 导出 alignment_osakt2_view.csv 和 candidate_sites_auto_v0.1.csv。

    返回 (aln_path, view_csv, cand_csv) 三个路径字符串。
    """

    fasta = Path((fasta_path or "").strip().strip('"').strip("'"))
    if not fasta.is_file():
        raise FileNotFoundError(f"找不到 FASTA 文件：{fasta}")

    out_dir = fasta.parent
    aln_path = out_dir / f"{fasta.stem}_OsAKT2.aln"

    cmd = [
        clustalo_cmd,
        "-i",
        str(fasta),
        "-o",
        str(aln_path),
        "--outfmt=clu",
        "--force",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            "Clustal-Omega 运行失败：\n"
            f"STDOUT:\n{result.stdout}\n\nSTDERR:\n{result.stderr}"
        )

    view_csv = out_dir / "alignment_osakt2_view.csv"
    cand_csv = out_dir / "candidate_sites_auto_v0.1.csv"

    if export_alignment_view and suggest_candidates:
        export_alignment_view(aln_path, view_csv)
        suggest_candidates(aln_path, cand_csv)
    else:
        _export_alignment_view_fallback(aln_path, view_csv)
        _suggest_candidates_fallback(aln_path, cand_csv)

    return str(aln_path), str(view_csv), str(cand_csv)


def run_osakt2_msa_wsl(fasta_path_win: str) -> Tuple[str, str, str]:
    """
    在 WSL 里用 Clustal-Omega 对 fasta_path_win 对应的多序列做比对，
    然后在 Windows 里用 msa_osakt2_tool 导出：
      1) XXX_OsAKT2.aln
      2) alignment_osakt2_view.csv
      3) candidate_sites_auto_v0.1.csv

    返回值（三个 Windows 路径，全是 str）：
        (aln_path_win, view_csv_win, cand_csv_win)
    """

    fasta_path = Path((fasta_path_win or "").strip().strip('"').strip("'"))
    if not fasta_path.is_file():
        raise FileNotFoundError(f"找不到 FASTA 文件：{fasta_path}")

    base_dir_win = str(fasta_path.parent)
    fasta_name = fasta_path.name
    stem = fasta_path.stem
    aln_name = f"{stem}_OsAKT2.aln"

    base_dir_wsl = hole_win_to_wsl(base_dir_win)
    exe = (CLUSTALO_WSL_EXE or "clustalo").strip()

    inner_cmd = (
        f'cd "{base_dir_wsl}" && '
        f'"{exe}" -i "{fasta_name}" -o "{aln_name}" '
        f"--outfmt=clu --force"
    )

    full_cmd = ["wsl", "bash", "-lc", inner_cmd]
    result = subprocess.run(
        full_cmd,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            "Clustal-Omega(WLS) 运行失败，请检查 CLUSTALO_WSL_EXE 或 WSL 配置：\n"
            f"命令: {' '.join(full_cmd)}\n\n"
            f"STDOUT:\n{result.stdout}\n\nSTDERR:\n{result.stderr}"
        )

    aln_path_win = str(fasta_path.with_name(aln_name))
    out_dir = fasta_path.parent

    view_csv = out_dir / "alignment_osakt2_view.csv"
    cand_csv = out_dir / "candidate_sites_auto_v0.1.csv"

    if export_alignment_view and suggest_candidates:
        export_alignment_view(Path(aln_path_win), view_csv)
        suggest_candidates(Path(aln_path_win), cand_csv)
    else:
        _export_alignment_view_fallback(Path(aln_path_win), view_csv)
        _suggest_candidates_fallback(Path(aln_path_win), cand_csv)

    return aln_path_win, str(view_csv), str(cand_csv)
def _extract_plddt_from_pdb(pdb_path: Path):
    """返回该 pdb 的逐残基 pLDDT 列表：[{chain, resi, resn, plddt}, ...]"""

    if _PDBParser is None:
        raise RuntimeError("需要 Biopython 才能解析 pLDDT（pip install biopython）")

    parser = _PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.stem, str(pdb_path))

    rows = []
    for model in structure:
        for chain in model:
            per_res = defaultdict(list)
            for res in chain:
                if res.id[0] != " ":  # 跳过水/异配体
                    continue
                for atom in res.get_atoms():
                    b = atom.get_bfactor()
                    if b is not None:
                        key = (chain.id, res.id)
                        per_res[key].append(float(b))

            for res in chain:
                if res.id[0] != " ":
                    continue
                key = (chain.id, res.id)
                if key in per_res:
                    vals = per_res[key]
                    rows.append(
                        {
                            "chain": chain.id,
                            "resi": res.id[1],
                            "resn": res.get_resname(),
                            "plddt": round(statistics.mean(vals), 2),
                        }
                    )
    return rows

def _summarize_plddt_rows(rows):
    """把逐残基列表压成一个 summary dict。"""

    if not rows:
        return {}
    vals = [r["plddt"] for r in rows]
    return {
        "plddt_mean": round(statistics.mean(vals), 2),
        "plddt_median": round(statistics.median(vals), 2),
        "plddt_lt50": sum(v < 50 for v in vals),
        "plddt_50_70": sum(50 <= v < 70 for v in vals),
        "plddt_ge70": sum(v >= 70 for v in vals),
    }

def plddt_summary_for_models(pdb_dir: str, models: list[str]) -> "_pd.DataFrame":
    """
    假设每个模型的 PDB 在 pdb_dir 下叫 <Model>.pdb，
    返回一个 DataFrame：Model + pLDDT 汇总列。
    """

    if _pd is None:
        raise RuntimeError("需要 pandas 才能生成 pLDDT 汇总表。")
    if _PDBParser is None:
        raise RuntimeError("需要 Biopython 才能解析 pLDDT。")

    pdb_dir = Path(pdb_dir)
    rows = []
    for m in models:
        pdb_path = pdb_dir / f"{m}.pdb"
        if not pdb_path.is_file():
            continue
        detail_rows = _extract_plddt_from_pdb(pdb_path)
        summary = _summarize_plddt_rows(detail_rows)
        if not summary:
            continue
        mean = summary["plddt_mean"]
        if mean >= 90:
            conf_class = "very_high_conf"
        elif mean >= 70:
            conf_class = "confident"
        elif mean >= 50:
            conf_class = "low_conf"
        else:
            conf_class = "very_low_conf"

        rows.append(
            {
                "Model": m,
                **summary,
                "ConfidenceClass": conf_class,
            }
        )

    if not rows:
        return _pd.DataFrame(
            columns=[
                "Model",
                "plddt_mean",
                "plddt_median",
                "plddt_lt50",
                "plddt_50_70",
                "plddt_ge70",
                "ConfidenceClass",
            ]
        )

    return _pd.DataFrame(rows)

def short_model_label(model: str) -> str:
    """
    尽量把模型名缩短：
      - 取第一个 '_' 前的部分：E174A_single → E174A
      - OsAKT2_WT → WT（你要的话可以单独再特判）
    """

    m = str(model)
    if "_" in m:
        head = m.split("_", 1)[0]
        if head.lower().startswith("osakt"):
            parts = m.split("_")
            if len(parts) >= 2:
                return parts[1]
        return head
    return m

def hole_plot_profiles(
        log_paths: Dict[str, str],
        out_dir: str,
        out_name: str = "hole_profiles.png",
) -> str:
    """把多条 profile 曲线绘制成图像。"""

    if _plt is None:
        raise RuntimeError("需要安装 matplotlib 才能绘图。")

    profiles = {model: hole_parse_profile(path) for model, path in log_paths.items()}

    os.makedirs(out_dir, exist_ok=True)

    n_models = len(profiles)
    fig_width = max(7.0, n_models * 0.6)
    fig, ax = _plt.subplots(figsize=(fig_width, 5))

    for model, prof in profiles.items():
        s_vals = [s for s, _ in prof]
        r_vals = [r for _, r in prof]
        label = short_model_label(model)
        ax.plot(s_vals, r_vals, label=label, linewidth=1.6)

    ax.axhline(0.0, linestyle="--", linewidth=0.8)
    ax.set_xlabel("s (Å)")
    ax.set_ylabel("Hole radius (Å)")

    handles, labels = _plt.gca().get_legend_handles_labels()
    _plt.legend(
        handles,
        labels,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        borderaxespad=0.0,
        fontsize=8,
        title="Model",
        title_fontsize=9,
    )

    fig.tight_layout(rect=[0, 0, 0.78, 1])

    out_path = os.path.join(out_dir, out_name)
    fig.savefig(out_path, dpi=300)
    _plt.close(fig)
    return out_path


def plot_basic_hole_metrics(hole_dir: str, out_dir: str | None = None):
    """
    读取 hole_min_table.csv，画基础对比图：
      - r_min_A 的柱状图
      - gate_length_A 的柱状图（如果存在）
    返回生成的图片路径列表。
    """

    if _pd is None or _plt is None:
        raise RuntimeError("需要安装 pandas 和 matplotlib 才能绘图。")

    hole_dir = (hole_dir or "").strip()
    if not hole_dir:
        raise ValueError("hole_dir 不能为空")

    csv_path = os.path.join(hole_dir, "hole_min_table.csv")
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"未找到 HOLE 汇总表：{csv_path}")

    df = _pd.read_csv(csv_path)
    if "Model" not in df.columns or "r_min_A" not in df.columns:
        raise ValueError("hole_min_table.csv 缺少 Model 或 r_min_A 列。")

    if out_dir is None:
        out_dir = hole_dir
    os.makedirs(out_dir, exist_ok=True)

    models = df["Model"].astype(str).tolist()
    labels = [short_model_label(m) for m in models]

    x = list(range(len(models)))
    out_paths: List[str] = []

    _plt.figure(figsize=(max(6, len(models) * 0.6), 4.5))
    _plt.bar(x, df["r_min_A"].tolist())
    _plt.xticks(x, labels, rotation=45, ha="right")
    _plt.xlabel("Model")
    _plt.ylabel("最小孔径 r_min (Å)")
    _plt.tight_layout()
    out1 = os.path.join(out_dir, "hole_r_min_bar.png")
    _plt.savefig(out1, dpi=300)
    _plt.close()
    out_paths.append(out1)

    if "gate_length_A" in df.columns:
        _plt.figure(figsize=(max(6, len(models) * 0.6), 4.5))
        _plt.bar(x, df["gate_length_A"].tolist())
        _plt.xticks(x, labels, rotation=45, ha="right")
        _plt.xlabel("Model")
        _plt.ylabel("gate 长度 (Å)")
        _plt.tight_layout()
        out2 = os.path.join(out_dir, "hole_gate_length_bar.png")
        _plt.savefig(out2, dpi=300)
        _plt.close()
        out_paths.append(out2)

    return out_paths

def normalize_path_for_chimerax(path: str) -> str:
    """
    把 Windows 路径变成 ChimeraX 好用的样子：
    去掉首尾引号和空格，反斜杠换成 /
    """
    p = path.strip().strip('"').strip("'")
    return p.replace("\\", "/")


def build_axis_cxc(
    wt_pdb_path: str,
    chain_id: str,
    residue_expr: str,
    out_dir: str,
    label: str = "axis",
) -> str:
    """根据 WT + 残基表达式生成“找轴”脚本。"""

    wt_pdb_path = (wt_pdb_path or "").strip()
    residue_expr = (residue_expr or "").strip()
    out_dir = (out_dir or "").strip()
    if not wt_pdb_path:
        raise ValueError("wt_pdb_path 不能为空")
    if not residue_expr:
        raise ValueError("residue_expr 不能为空")
    if not out_dir:
        raise ValueError("out_dir 不能为空")

    wt_cx = normalize_path_for_chimerax(wt_pdb_path)
    out_dir_cx = normalize_path_for_chimerax(out_dir)

    safe_label = (label or "axis").strip() or "axis"
    log_name = f"axis_{safe_label}.log"
    log_path_cx = f"{out_dir_cx}/{log_name}"

    lines = [
        f"# === 自动找轴工具: {safe_label} ===",
        "close all",
        f"open \"{wt_cx}\"",
        f"select #1/{chain_id or 'A'}:{residue_expr}",
        "view sel",
        "clip near 8; clip far 60",
        "log clear",
        "measure center sel mark true radius 0.6 color yellow",
        f"log save \"{log_path_cx}\" executableLinks false",
        "getcrd sel",
        "# 运行完后，到输出目录里找 axis_*.log 并交给 Python 解析。",
    ]

    os.makedirs(out_dir, exist_ok=True)
    axis_cxc_path = os.path.join(out_dir, f"{safe_label}_axis.cxc")
    with open(axis_cxc_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    return axis_cxc_path


def parse_axis_log(log_path: str) -> Tuple[float, float, float]:
    """从 axis log 中解析出质心坐标。"""

    if not os.path.exists(log_path):
        raise FileNotFoundError(f"找不到 log 文件：{log_path}")

    pattern = re.compile(r"Center of mass[^=]*=\s*\(([^)]+)\)")
    with open(log_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if "Center of mass" not in line:
                continue
            match = pattern.search(line)
            if not match:
                continue
            coord_str = match.group(1)
            parts = [p.strip() for p in coord_str.split(",")]
            if len(parts) != 3:
                continue
            try:
                x, y, z = map(float, parts)
            except ValueError:
                continue
            return x, y, z

    raise ValueError(f"在 log 里没找到质心坐标: {log_path}")


def _split_multi_value(value: str) -> List[str]:
    """把用户输入里用逗号、空格或分号分隔的值拆成列表。"""
    if not value:
        return []
    # 允许中文逗号（有些输入法会产生）
    value = value.replace("，", ",")
    parts = re.split(r"[;,\s]+", value)
    return [p for p in (part.strip() for part in parts) if p]


def build_mutation_cxc(
        wt_pdb_path: str,
        mut_label: str,
        chain: str,
        residue: str,
        new_aa: str,
        out_dir: str,
) -> str:
    """
    自动生成突变体的 cxc，用 ChimeraX 一键构建突变体。
    out_dir/MUT/ 下生成:
        OsAKT_xxx.cxc
        OsAKT_xxx.pdb (由 ChimeraX 生成)
    """
    # 规范一下路径和参数
    out_dir_cx = normalize_path_for_chimerax(out_dir)
    wt_cx = normalize_path_for_chimerax(wt_pdb_path)

    mut_dir = os.path.join(out_dir, "MUT")
    os.makedirs(mut_dir, exist_ok=True)
    mut_dir_cx = normalize_path_for_chimerax(mut_dir)

    residues = _split_multi_value(str(residue).strip())
    aas_raw = [item.upper() for item in _split_multi_value(new_aa)]
    chains = _split_multi_value(chain) or ["A"]

    if not residues:
        raise ValueError("至少需要提供一个残基号。")
    if not aas_raw:
        raise ValueError("至少需要提供一个目标氨基酸。")

    if len(chains) == 1:
        chains *= len(residues)
    elif len(chains) != len(residues):
        raise ValueError("链 ID 数量需要和残基号数量一致（或只填一个链 ID）。")

    if len(aas_raw) == 1:
        aas_raw *= len(residues)
    elif len(aas_raw) != len(residues):
        raise ValueError("目标氨基酸数量需要和残基号数量一致（或只填一个氨基酸）。")

    mutation_steps = []
    for ch, res, aa_raw in zip(chains, residues, aas_raw):
        aa_raw = aa_raw.strip()
        if len(aa_raw) == 1:
            if aa_raw not in ONE_TO_THREE:
                raise ValueError(f"不认识的氨基酸一字母代码：{aa_raw}")
            aa = ONE_TO_THREE[aa_raw]
        elif len(aa_raw) == 3:
            aa = aa_raw
        else:
            raise ValueError(
                f"氨基酸名请用一字母或三字母，比如 N / ASN（现在是：{aa_raw}）"
            )
        mutation_steps.append({
            "chain": ch or "A",
            "residue": res,
            "aa": aa,
        })

    cxc_path = os.path.join(mut_dir, f"{mut_label}.cxc")

    mut_lines = ["# ChimeraX swapaa 语法：swapaa <残基选择> <目标氨基酸>"]
    for idx, step in enumerate(mutation_steps, start=1):
        desc = f"{step['chain']}:{step['residue']} -> {step['aa']}"
        mut_lines.append(f"# Step {idx}: {desc}")
        mut_lines.append(
            "swapaa #1/{chain}:{residue} {aa}".format(
                chain=step["chain"],
                residue=step["residue"],
                aa=step["aa"],
            )
        )

    script = f"""# === 自动突变：{mut_label} ===
    close all

    open "{wt_cx}"

    {os.linesep.join(mut_lines)}

    save "{mut_dir_cx}/{mut_label}.pdb"

    close all
    """

    with open(cxc_path, "w", encoding="utf-8") as f:
        f.write(script)

    return cxc_path


def build_cxc_script(
    wt_pdb_path: str,
    mutant_list: List[Dict[str, str]],
    chain_id: str,
    residue_expr: str,
    out_dir: str,
    features: Dict[str, bool],
    policy: str = "Standard",
    roi_expr: Optional[str] = None,
    roi_view_name: str = "ROI",
) -> str:
    """
    根据配置拼出 ChimeraX .cxc 脚本文本。
    residue_expr: 目标残基表达式，例如 "298,299,300" 或 "298-305"
    roi_expr: ROI 残基表达式，用于局部视角 5/6，例如 "283,286,291,298-300"
    features: {
        "full_coulombic": bool,
        "contacts": bool,
        "hbonds": bool,
        "sasa": bool,
        "sites_contacts": bool,
        "sites_coulombic": bool,
    }
    """
    out_dir_cx = normalize_path_for_chimerax(out_dir)
    wt_pdb_cx = normalize_path_for_chimerax(wt_pdb_path)

    # 清理残基表达式（去空格）
    residue_expr = (residue_expr or "").replace(" ", "")
    roi_expr = (roi_expr or "").replace(" ", "")
    chain_id = (chain_id or "A").strip()

    if (features.get("sites_contacts") or features.get("sites_coulombic")) and not roi_expr:
        raise ValueError("roi_expr 不能为空：已勾选 5/6，但没有 ROI 残基表达式。")

    if features.get("sites_contacts") or features.get("sites_coulombic"):
        os.makedirs(os.path.join(out_dir, "gate_sites"), exist_ok=True)

    lines: List[str] = []
    add = lines.append

    # Header
    add("# === Auto-generated ChimeraX script ===")
    add("# 由 Python 工具生成，自动完成 open / 对齐 / 出图 等步骤。")
    add("")
    add("# WT PDB: %s" % wt_pdb_cx)
    if mutant_list:
        labels = ", ".join(m.get("label", "?") for m in mutant_list)
        add("# Mutants: %s" % labels)
    add("")

    # 打开所有模型（突变体在前，WT 在最后）
    add("# ——打开所有 PDB——")
    for i, m in enumerate(mutant_list, start=1):
        pdb_cx = normalize_path_for_chimerax(m["pdb"])
        m["model_id"] = i
        add('open "%s"' % pdb_cx)
    wt_id = len(mutant_list) + 1
    add('open "%s"' % wt_pdb_cx)
    last_id = wt_id
    add("")

    # 对齐突变体到 WT
    if mutant_list:
        add("# ——把全部突变体对齐到 WT——")
        if len(mutant_list) == 1:
            # 只有一个突变体：mm #1 to #2 这种写法
            mid = mutant_list[0]["model_id"]
            add("mm #%d to #%d showAlignment false" % (mid, wt_id))
        else:
            # 突变体在前，WT 在最后，模型 ID 一定是 1..N 连续
            max_id = wt_id - 1
            # 对应你操作文档里的 mm #1-4 to #5
            add("mm #1-%d to #%d showAlignment false" % (max_id, wt_id))
        add("")


    # 相机与画幅
    add("# ——画幅 / 相机 / 灯光——")
    add("windowsize 2200 1400")
    add("camera ortho")
    add("lighting soft")
    add("")

    def _cx_dir(label: str) -> str:
        return normalize_path_for_chimerax(_model_dir(out_dir, label))

    # FULL 静电势
    if features.get("full_coulombic"):
        add("# ===================== FULL 视角 + 静电势 =====================")
        # WT
        add("# ——WT FULL——")
        add("hide #1-%d" % last_id)
        add("show #%d" % wt_id)
        add("surface #%d" % wt_id)
        add("coulombic #%d range -10,10" % wt_id)
        add('save "%s/WT_coulombic.png" width 2200 supersample 3' % _cx_dir("WT"))
        add("surface hide")
        add("")
        # 每个突变体
        for m in mutant_list:
            label = m.get("label", "MUT")
            mid = m["model_id"]
            add("# ——%s FULL——" % label)
            add("hide #1-%d" % last_id)
            add("show #%d" % mid)
            add("surface #%d" % mid)
            add("coulombic #%d range -10,10" % mid)
            add('save "%s/%s_coulombic.png" width 2200 supersample 3' % (_cx_dir(label), label))
            add("surface hide")
            add("")

    # 是否需要近景视角（局部分析）
    need_local_view = any(
        features.get(flag)
        for flag in ("contacts", "hbonds", "sasa")
    )
    if need_local_view and residue_expr:
        add("# ===================== 目标残基近景视角 =====================")
        add("# 以 WT 为基准，把视角对准指定残基附近")
        add("hide #1-%d" % last_id)
        add("show #%d" % wt_id)
        add("select #%d/%s:%s" % (wt_id, chain_id, residue_expr))
        add("view sel")
        add("clip near 8; clip far 60")
        add("view name TARGET")
        add("delete pseudobonds")
        add("surface hide")
        add("")

    # 接触图
    if features.get("contacts") and residue_expr:
        add("# ===================== 接触图（≤4 Å） =====================")
        add("surface hide")
        all_models = mutant_list + [{"label": "WT", "model_id": wt_id}]
        for m in all_models:
            label = m.get("label", "MUT")
            mid = m["model_id"]
            add("# ——%s——" % label)
            add("hide #1-%d" % last_id)
            add("show #%d" % mid)
            add("view TARGET")
            add("delete pseudobonds")
            add("select #%d/%s:%s" % (mid, chain_id, residue_expr))
            add("contacts sel distanceOnly 4 reveal true intermodel false")
            add("size #%d pseudobondradius 0.35" % mid)
            add("color yellow pseudobonds")
            add('save "%s/%s_contacts.png" width 2200 supersample 3' % (_cx_dir(label), label))
            add("delete pseudobonds")
            add("")

    # 氢键图
    if features.get("hbonds") and residue_expr:
        add("# ===================== 氢键 + 文本日志 =====================")
        add("surface hide")
        all_models = mutant_list + [{"label": "WT", "model_id": wt_id}]
        for m in all_models:
            label = m.get("label", "MUT")
            mid = m["model_id"]
            add("# ——%s——" % label)
            add("hide #1-%d" % last_id)
            add("show #%d" % mid)
            add("view TARGET")
            add("delete pseudobonds")
            add("select #%d/%s:%s" % (mid, chain_id, residue_expr))
            add("hbonds sel reveal true interModel false showDist true")
            add("size #%d pseudobondradius 0.35" % mid)
            add("color yellow pseudobonds")
            if policy in ("Standard", "Full"):
                add('save "%s/%s_hbonds.png" width 2200 supersample 3' % (_cx_dir(label), label))
            add(
                'hbonds #%d/%s:%s interModel false intraModel true '
                'log true saveFile "%s/%s_hbonds.txt"'
                % (mid, chain_id, residue_expr, _cx_dir(label), label)
            )
            add("delete pseudobonds")
            add("")

    # SASA
    if features.get("sasa") and residue_expr:
        add("# ===================== SASA（溶剂可及面积） =====================")
        all_models = mutant_list + [{"label": "WT", "model_id": wt_id}]
        for m in all_models:
            label = m.get("label", "MUT")
            mid = m["model_id"]
            add("# ——%s——" % label)
            add("log clear")
            add("measure sasa #%d/%s:%s" % (mid, chain_id, residue_expr))
            add("info residues #%d/%s:%s attribute area" % (mid, chain_id, residue_expr))
            add('log save "%s/%s_sasa.html" executableLinks false' % (_cx_dir(label), label))
            add("")

    # ROI 局部批处理
    if (features.get("sites_contacts") or features.get("sites_coulombic")) and roi_expr:
        sites_dir_cx = f"{out_dir_cx}/gate_sites"

        add("# ===================== ROI 统一视角批处理 =====================")
        add("hide #1-%d" % last_id)
        add("show #%d" % wt_id)
        add(f"show #{wt_id}/{chain_id}:{roi_expr}")
        add(f"view #{wt_id}/{chain_id}:{roi_expr}")
        add("scale 1.8")
        add("clip near 5")
        add("clip far 50")
        add(f"view name \"{roi_view_name}\"")
        add("delete pseudobonds")
        add("surface hide")
        add("")

        site_models = [{"label": "WT", "model_id": wt_id}, *mutant_list]

        for m in site_models:
            label = m.get("label", "MUT")
            mid = m["model_id"]
            add("# ——%s ROI——" % label)
            add("hide #1-%d" % last_id)
            add(f"show #{mid}")
            add(f"view \"{roi_view_name}\"")
            add("hide atoms")
            add(f"cartoon #{mid}")
            add(f"color #{mid} white")

            if features.get("sites_contacts"):
                add("surface hide")
                add(f"show #{mid}/{chain_id}:{roi_expr}")
                add(f"color #{mid}/{chain_id}:{roi_expr} yellow")
                add(f"style #{mid}/{chain_id}:{roi_expr} stick")
                add(f"contacts #{mid}/{chain_id}:{roi_expr} distanceOnly 4 reveal true")
                add(
                    f'save "{sites_dir_cx}/{label}_sites_contacts.png" '
                    "width 1600 height 1000 supersample 3"
                )
                add("delete pseudobonds")

            if features.get("sites_coulombic"):
                add(f"surface #{mid}")
                add(f"coulombic protein surfaces #{mid}")
                add(
                    f'save "{sites_dir_cx}/{label}_sites_coulombic.png" '
                    "width 1600 height 1000 supersample 3"
                )
                add("surface hide")

            add("")

    add("# === 脚本结束 ===")
    return "\n".join(lines) + "\n"

def parse_sasa_html(path: str) -> Dict[str, float]:
    """解析 ChimeraX 输出的 *_sasa.html，得到残基 SASA。"""

    if not os.path.exists(path):
        raise FileNotFoundError(f"找不到 SASA 文件：{path}")

    with open(path, encoding="utf-8", errors="ignore") as f:
        text = f.read()

    text_no_html = re.sub(r"<[^>]+>", " ", text)
    pattern = re.compile(r"[#/]*([\w:]+).*?area\s+([0-9.]+)", re.IGNORECASE)
    sasa: Dict[str, float] = {}

    for raw_line in text_no_html.splitlines():
        line = raw_line.strip()
        if not line or "area" not in line.lower():
            continue
        match = pattern.search(line)
        if not match:
            continue
        resid = match.group(1)
        try:
            area = float(match.group(2))
        except ValueError:
            continue
        sasa[resid] = area

    if not sasa:
        raise ValueError(f"在 {path} 中没有解析到 SASA 数据。")

    return sasa

def parse_hbonds_txt(path: str) -> int:
    """统计氢键日志（txt）中的氢键数量。"""

    if not os.path.exists(path):
        raise FileNotFoundError(f"找不到氢键日志：{path}")

    count = 0
    with open(path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            striped = line.strip()
            if not striped or striped.startswith("#"):
                continue
            count += 1
    return count

def summarize_sasa_hbonds(out_dir: str) -> Tuple[str, str]:
    """扫描 *_sasa.html / *_hbonds.txt，输出汇总 CSV。"""

    if _pd is None:
        raise RuntimeError("需要安装 pandas 才能汇总 SASA/H-bonds。")

    out_dir = (out_dir or "").strip()
    if not out_dir:
        raise ValueError("out_dir 不能为空")

    sasa_files = sorted(
        glob.glob(os.path.join(out_dir, "**", "*_sasa.html"), recursive=True)
    )
    if not sasa_files:
        raise ValueError(f"在 {out_dir}（含子目录）没找到任何 *_sasa.html 文件。")

    hbonds_files = sorted(
        glob.glob(os.path.join(out_dir, "**", "*_hbonds.txt"), recursive=True)
    )
    if not hbonds_files:
        raise ValueError(f"在 {out_dir}（含子目录）没找到任何 *_hbonds.txt 文件。")

    summary_rows: List[Dict[str, float]] = []
    detail_rows: List[Dict[str, str]] = []
    all_resids = set()
    hbonds_map = {
        os.path.basename(path).replace("_hbonds.txt", ""): path for path in hbonds_files
    }

    for sasa_path in sasa_files:
        label = os.path.basename(sasa_path).replace("_sasa.html", "")
        sasa_dict = parse_sasa_html(sasa_path)
        hb_path = hbonds_map.get(label)
        hbonds = parse_hbonds_txt(hb_path) if hb_path and os.path.exists(hb_path) else 0
        total_sasa = sum(sasa_dict.values())

        row: Dict[str, float] = {
            "Model": label,
            "HBonds": hbonds,
            "Total_SASA": total_sasa,
        }

        for resid, area in sasa_dict.items():
            col = f"SASA_{resid}"
            row[col] = area
            all_resids.add(resid)
            detail_rows.append({
                "Model": label,
                "Residue": resid,
                "SASA": area,
            })

        summary_rows.append(row)

    df_summary = _pd.DataFrame(summary_rows)
    resid_cols = [f"SASA_{resid}" for resid in sorted(all_resids)]
    for col in resid_cols:
        if col not in df_summary.columns:
            df_summary[col] = _pd.NA
    ordered_cols = ["Model", "HBonds", "Total_SASA"] + resid_cols
    df_summary = df_summary[ordered_cols].sort_values("Model")

    df_detail = _pd.DataFrame(detail_rows).sort_values(["Residue", "Model"])

    summary_csv = table_path(out_dir, "sasa_hbonds_summary.csv")
    detail_csv = table_path(out_dir, "sasa_per_residue.csv")

    df_summary.to_csv(summary_csv, index=False)
    df_detail.to_csv(detail_csv, index=False)

    return summary_csv, detail_csv
