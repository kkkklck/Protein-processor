from typing import List, Dict
import os, re
import shutil

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


def normalize_path_for_chimerax(path: str) -> str:
    """
    把 Windows 路径变成 ChimeraX 好用的样子：
    去掉首尾引号和空格，反斜杠换成 /
    """
    p = path.strip().strip('"').strip("'")
    return p.replace("\\", "/")


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
