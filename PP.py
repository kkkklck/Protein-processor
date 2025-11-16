from typing import List, Dict, Tuple
import os, re, glob
import shutil
import subprocess

try:
    import pandas as _pd
    import matplotlib.pyplot as _plt
except Exception:  # pragma: no cover - optional deps
    _pd = None
    _plt = None

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

# ===== HOLE / WSL 默认配置（你只需要在自己电脑上改一次） =====
# 1. 在 WSL 里执行 `conda info --base` 得到 conda 的 base 路径，比如：
#    /home/k/miniforge3
# 2. 把 HOLE_WSL_CONDA_INIT 改成  <base>/etc/profile.d/conda.sh
# 3. 把 HOLE_WSL_CONDA_ENV 改成你安装 hole 的那个环境名（例如 "hole_env"）

HOLE_WSL_CONDA_INIT = "$HOME/miniforge3/etc/profile.d/conda.sh"  # ← 根据实际路径改
HOLE_WSL_CONDA_ENV = "hole_env"  # ← 根据实际 env 名改
HOLE_WSL_EXE = "hole"  # env 里 HOLE 的命令名


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
        s_min, r_min = min(prof, key=lambda item: item[1])
        min_table_rows.append({"Model": model, "s_min_A": s_min, "r_min_A": r_min})
        min_summary_rows.append({
            "Model": model,
            "Summary": f"Minimum radius found: {r_min:8.3f} Å at s = {s_min:8.3f} Å",
        })

    os.makedirs(out_dir, exist_ok=True)
    _pd.DataFrame(profiles_rows).to_csv(os.path.join(out_dir, "hole_profile_samples.csv"), index=False)
    _pd.DataFrame(min_table_rows).to_csv(os.path.join(out_dir, "hole_min_table.csv"), index=False)
    _pd.DataFrame(min_summary_rows).to_csv(os.path.join(out_dir, "hole_min_summary.csv"), index=False)


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
    _plt.figure(figsize=(6, 4))
    for model, prof in profiles.items():
        s_vals = [s for s, _ in prof]
        r_vals = [r for _, r in prof]
        _plt.plot(s_vals, r_vals, label=model)

    _plt.axhline(0.0, linestyle="--", linewidth=0.8)
    _plt.xlabel("s (Å)")
    _plt.ylabel("Hole radius (Å)")
    _plt.legend()
    _plt.tight_layout()

    out_path = os.path.join(out_dir, out_name)
    _plt.savefig(out_path, dpi=300)
    _plt.close()
    return out_path

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
) -> str:
    """
    根据配置拼出 ChimeraX .cxc 脚本文本。
    residue_expr: 目标残基表达式，例如 "298,299,300" 或 "298-305"
    features: {
        "full_coulombic": bool,
        "contacts": bool,
        "hbonds": bool,
        "sasa": bool,
    }
    """
    out_dir_cx = normalize_path_for_chimerax(out_dir)
    wt_pdb_cx = normalize_path_for_chimerax(wt_pdb_path)

    # 清理残基表达式（去空格）
    residue_expr = (residue_expr or "").replace(" ", "")
    chain_id = (chain_id or "A").strip()

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
        ids_str = ",".join("#%d" % m["model_id"] for m in mutant_list)
        add("mm %s to #%d showAlignment false" % (ids_str, wt_id))
        add("")

    # 相机与画幅
    add("# ——画幅 / 相机 / 灯光——")
    add("windowsize 2200 1400")
    add("camera ortho")
    add("lighting soft")
    add("")

    # FULL 静电势
    if features.get("full_coulombic"):
        add("# ===================== FULL 视角 + 静电势 =====================")
        # WT
        add("# ——WT FULL——")
        add("hide #1-%d" % last_id)
        add("show #%d" % wt_id)
        add("surface #%d" % wt_id)
        add("coulombic #%d range -10,10" % wt_id)
        add('save "%s/WT_coulombic.png" width 2200 supersample 3' % out_dir_cx)
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
            add('save "%s/%s_coulombic.png" width 2200 supersample 3' % (out_dir_cx, label))
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
        add("")

    # 接触图
    if features.get("contacts") and residue_expr:
        add("# ===================== 接触图（≤4 Å） =====================")
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
            add('save "%s/%s_contacts.png" width 2200 supersample 3' % (out_dir_cx, label))
            add("delete pseudobonds")
            add("")

    # 氢键图
    if features.get("hbonds") and residue_expr:
        add("# ===================== 氢键 + 文本日志 =====================")
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
            add('save "%s/%s_hbonds.png" width 2200 supersample 3' % (out_dir_cx, label))
            add(
                'hbonds #%d/%s:%s interModel false intraModel true '
                'log true saveFile "%s/%s_hbonds.txt"'
                % (mid, chain_id, residue_expr, out_dir_cx, label)
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
            add('log save "%s/%s_sasa.html" executableLinks false' % (out_dir_cx, label))
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

    sasa_files = sorted(glob.glob(os.path.join(out_dir, "*_sasa.html")))
    if not sasa_files:
        raise ValueError(f"在 {out_dir} 没找到任何 *_sasa.html 文件。")

    summary_rows: List[Dict[str, float]] = []
    detail_rows: List[Dict[str, str]] = []
    all_resids = set()

    for sasa_path in sasa_files:
        label = os.path.basename(sasa_path).replace("_sasa.html", "")
        sasa_dict = parse_sasa_html(sasa_path)
        hb_path = os.path.join(out_dir, f"{label}_hbonds.txt")
        hbonds = parse_hbonds_txt(hb_path) if os.path.exists(hb_path) else 0
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

    summary_csv = os.path.join(out_dir, "sasa_hbonds_summary.csv")
    detail_csv = os.path.join(out_dir, "sasa_per_residue.csv")

    df_summary.to_csv(summary_csv, index=False)
    df_detail.to_csv(detail_csv, index=False)

    return summary_csv, detail_csv
