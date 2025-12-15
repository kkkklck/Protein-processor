import csv
import os
import shutil
import tkinter as tk
from tkinter import filedialog, messagebox

from PP import (
    build_axis_cxc,
    build_cxc_script,
    build_mutation_cxc,
    parse_axis_log,
    hole_write_input,
    hole_run_in_wsl,
    hole_summarize_logs,
    hole_plot_profiles,
    summarize_sasa_hbonds,
    plot_basic_hole_metrics,
    merge_all_metrics,
    score_metrics_file,
    make_stage3_table,
    run_osakt2_msa_wsl,
)


def create_gui():
    root = tk.Tk()
    root.title("ChimeraX .cxc 自动生成器（GUI 版）")

    # 整体窗口稍微宽一点
    root.geometry("900x600")

    # ===== 外层可滚动区域（大滚动条） =====
    canvas = tk.Canvas(root, highlightthickness=0)
    vbar = tk.Scrollbar(root, orient="vertical", command=canvas.yview)
    canvas.configure(yscrollcommand=vbar.set)

    vbar.pack(side="right", fill="y")
    canvas.pack(side="left", fill="both", expand=True)

    # 真正放控件的地方，全都往这个 content 里塞
    content = tk.Frame(canvas)
    content_window = canvas.create_window((0, 0), window=content, anchor="nw")

    # 内容大小变化时，更新滚动区域
    def _on_content_config(event):
        canvas.configure(scrollregion=canvas.bbox("all"))

    content.bind("<Configure>", _on_content_config)

    # 画布大小变化时，让 content 宽度跟着变，避免只占左侧一小条
    def _on_canvas_config(event):
        canvas.itemconfig(content_window, width=event.width)

    canvas.bind("<Configure>", _on_canvas_config)

    # 绑定鼠标滚轮（Windows）
    def _on_mousewheel(event):
        canvas.yview_scroll(int(-event.delta / 120), "units")

    canvas.bind_all("<MouseWheel>", _on_mousewheel)

    # ===== WT 区 =====
    wt_frame = tk.LabelFrame(content, text="WT PDB", padx=8, pady=8)
    wt_frame.pack(fill="x", padx=10, pady=5)

    wt_path_var = tk.StringVar()

    tk.Label(wt_frame, text="WT PDB 路径：").grid(row=0, column=0, sticky="w")
    wt_entry = tk.Entry(wt_frame, textvariable=wt_path_var, width=70)
    wt_entry.grid(row=0, column=1, sticky="w")

    def browse_wt():
        path = filedialog.askopenfilename(
            title="选择 WT PDB 文件",
            filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")]
        )
        if path:
            wt_path_var.set(path)

    tk.Button(wt_frame, text="浏览", command=browse_wt).grid(row=0, column=2, padx=5)

    # ===== 模式切换 =====
    research_container = tk.Frame(content)
    mutate_container = tk.Frame(content)
    hole_container = tk.Frame(content)

    mode_var = tk.StringVar(value="research")  # "research" / "mutate" / "hole"

    def update_mode():
        mode = mode_var.get()
        for frame in (research_container, mutate_container, hole_container):
            frame.pack_forget()
        if mode == "research":
            research_container.pack(fill="both", expand=True)
        elif mode == "mutate":
            mutate_container.pack(fill="both", expand=True)
        else:
            hole_container.pack(fill="both", expand=True)

    mode_frame = tk.Frame(wt_frame)
    mode_frame.grid(row=0, column=3, padx=10, sticky="w")
    tk.Label(mode_frame, text="模式：").pack(side="left")
    tk.Radiobutton(
        mode_frame,
        text="研究",
        value="research",
        variable=mode_var,
        command=update_mode
    ).pack(side="left")
    tk.Radiobutton(
        mode_frame,
        text="突变",
        value="mutate",
        variable=mode_var,
        command=update_mode
    ).pack(side="left", padx=(4, 0))
    tk.Radiobutton(
        mode_frame,
        text="HOLE",
        value="hole",
        variable=mode_var,
        command=update_mode
    ).pack(side="left", padx=(4, 0))

    # ===== 研究模式：突变体区 =====
    mutants_outer = tk.LabelFrame(research_container, text="突变体 PDB（可选，多条）", padx=8, pady=8)
    mutants_outer.pack(fill="both", padx=10, pady=5, expand=True)

    scroll = ScrollableFrame(mutants_outer)
    scroll.pack(fill="both", expand=True)

    mutants_frame = scroll.scrollable_frame

    mutant_rows = []

    def add_mutant_row(default_label=None):
        idx = len(mutant_rows) + 1
        row_frame = tk.Frame(mutants_frame)
        row_frame.pack(fill="x", pady=2)

        label_var = tk.StringVar(value=default_label or f"MUT{idx}")
        pdb_var = tk.StringVar()

        tk.Label(row_frame, text=f"{idx}. 标签：").grid(row=0, column=0, sticky="w")
        tk.Entry(row_frame, textvariable=label_var, width=10).grid(row=0, column=1, sticky="w", padx=(0, 10))

        tk.Label(row_frame, text="PDB：").grid(row=0, column=2, sticky="w")
        tk.Entry(row_frame, textvariable=pdb_var, width=50).grid(row=0, column=3, sticky="w")

        def browse_mutant():
            path = filedialog.askopenfilename(
                title="选择突变体 PDB 文件",
                filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")]
            )
            if path:
                pdb_var.set(path)

        tk.Button(row_frame, text="浏览", command=browse_mutant).grid(row=0, column=4, padx=5)

        row_dict = {
            "label_var": label_var,
            "pdb_var": pdb_var,
            "frame": row_frame,
        }

        def delete_row():
            if row_dict in mutant_rows:
                mutant_rows.remove(row_dict)
            row_frame.destroy()

        del_btn = tk.Button(row_frame, text="删除", fg="red", command=delete_row)
        del_btn.grid(row=0, column=5, padx=5)
        row_dict["del_btn"] = del_btn

        mutant_rows.append(row_dict)
    # 默认先给两行（方便你现在 DMI / DMT）
    add_mutant_row("DMI")
    add_mutant_row("DMT")

    # 注意：按钮放在滚动区外面，这样永远在列表最下面
    tk.Button(
        mutants_outer,
        text="添加突变体",
        command=add_mutant_row
    ).pack(anchor="w", padx=10, pady=4)

    # ===== 功能选择 =====
    feature_frame = tk.LabelFrame(research_container, text="ChimeraX 自动化", padx=8, pady=8)
    feature_frame.pack(fill="x", padx=10, pady=5)

    full_var = tk.IntVar(value=1)
    contacts_var = tk.IntVar(value=1)
    hbonds_var = tk.IntVar(value=1)
    sasa_var = tk.IntVar(value=1)
    sites_contacts_var = tk.IntVar(value=1)
    sites_coulombic_var = tk.IntVar(value=1)

    tk.Checkbutton(
        feature_frame, text="1. FULL 静电势图", variable=full_var
    ).grid(row=0, column=0, sticky="w")
    tk.Checkbutton(
        feature_frame, text="2. 近景接触图（≤4 Å）", variable=contacts_var
    ).grid(row=0, column=1, sticky="w")
    tk.Checkbutton(
        feature_frame, text="3. 近景氢键图 + 文本日志", variable=hbonds_var
    ).grid(row=1, column=0, sticky="w")
    tk.Checkbutton(
        feature_frame, text="4. 目标残基 SASA", variable=sasa_var
    ).grid(row=1, column=1, sticky="w")
    tk.Checkbutton(
        feature_frame,
        text="5. ROI 局部接触图（自定义）",
        variable=sites_contacts_var,
    ).grid(row=2, column=0, sticky="w", pady=(4, 0))
    tk.Checkbutton(
        feature_frame,
        text="6. ROI 静电势图（自定义）",
        variable=sites_coulombic_var,
    ).grid(row=2, column=1, sticky="w", pady=(4, 0))

    tk.Label(
        feature_frame,
        text="说明：2 / 3 / 4 需要指定“目标残基”；5 / 6 必须指定 ROI 残基表达式；只勾 1 时可以都不填。",
        fg="#555"
    ).grid(row=3, column=0, columnspan=2, sticky="w", pady=(4, 0))

    # ===== 目标残基 & 链 =====
    target_frame = tk.LabelFrame(research_container, text="目标残基（可选）", padx=8, pady=8)
    target_frame.pack(fill="x", padx=10, pady=5)

    chain_var = tk.StringVar(value="A")
    residue_expr_var = tk.StringVar()
    roi_expr_var = tk.StringVar()

    tk.Label(target_frame, text="链 ID：").grid(row=0, column=0, sticky="w")
    tk.Entry(target_frame, textvariable=chain_var, width=5).grid(row=0, column=1, sticky="w")

    tk.Label(target_frame, text="残基表达式：").grid(row=0, column=2, sticky="w", padx=(15, 0))
    tk.Entry(target_frame, textvariable=residue_expr_var, width=40).grid(row=0, column=3, sticky="w")

    tk.Label(
        target_frame,
        text="例：298,299,300 或 298-305。",
        fg="#555"
    ).grid(row=1, column=0, columnspan=4, sticky="w", pady=(4, 0))

    tk.Label(target_frame, text="ROI 残基表达式（用于 5/6）：").grid(row=2, column=0, sticky="w", pady=(8, 0))
    tk.Entry(target_frame, textvariable=roi_expr_var, width=40).grid(row=2, column=1, columnspan=3, sticky="w", pady=(8, 0))
    tk.Label(
        target_frame,
        text="例：283,286,291,298-300 或 45-60,120,155-170。",
        fg="#555",
    ).grid(row=3, column=0, columnspan=4, sticky="w", pady=(4, 0))

    # ===== 输出设置 =====
    out_frame = tk.LabelFrame(research_container, text="输出位置", padx=8, pady=8)
    out_frame.pack(fill="x", padx=10, pady=5)

    out_dir_var = tk.StringVar(value=os.path.join("D:\\", "demo"))
    cxc_path_var = tk.StringVar(value=os.path.join("D:\\", "demo", "auto_chimerax.cxc"))

    tk.Label(out_frame, text="图片 / 文本输出目录：").grid(row=0, column=0, sticky="w")
    tk.Entry(out_frame, textvariable=out_dir_var, width=60).grid(row=0, column=1, sticky="w")

    def browse_out_dir():
        path = filedialog.askdirectory(title="选择输出目录")
        if path:
            out_dir_var.set(path)

    tk.Button(out_frame, text="浏览", command=browse_out_dir).grid(row=0, column=2, padx=5)

    tk.Label(out_frame, text=".cxc 脚本保存为：").grid(row=1, column=0, sticky="w", pady=(6, 0))
    tk.Entry(out_frame, textvariable=cxc_path_var, width=60).grid(row=1, column=1, sticky="w", pady=(6, 0))

    def browse_cxc_file():
        path = filedialog.asksaveasfilename(
            title="保存 .cxc 文件",
            defaultextension=".cxc",
            filetypes=[("ChimeraX script", "*.cxc"), ("All files", "*.*")]
        )
        if path:
            cxc_path_var.set(path)

    tk.Button(out_frame, text="浏览", command=browse_cxc_file).grid(row=1, column=2, padx=5, pady=(6, 0))

    # ===== MSA 自动候选 =====
    msa_frame = tk.LabelFrame(
        research_container,
        text="MSA 自动候选小工具（Clustal-Omega + OsAKT2）",
        padx=8,
        pady=8,
    )
    msa_frame.pack(fill="x", padx=10, pady=5)

    msa_fasta_var = tk.StringVar()

    tk.Label(msa_frame, text="FASTA 文件：").grid(row=0, column=0, sticky="w")
    tk.Entry(msa_frame, textvariable=msa_fasta_var, width=60).grid(
        row=0, column=1, sticky="w"
    )

    def browse_msa_fasta():
        path = filedialog.askopenfilename(
            title="选择多序列 FASTA 文件",
            filetypes=[
                ("FASTA", "*.fasta *.fa *.faa *.fna *.fas"),
                ("All files", "*.*"),
            ],
        )
        if path:
            msa_fasta_var.set(path)

    tk.Button(msa_frame, text="浏览", command=browse_msa_fasta).grid(
        row=0, column=2, padx=5
    )

    def on_run_msa():
        fasta = msa_fasta_var.get().strip()
        if not fasta:
            messagebox.showerror("缺少 FASTA", "请先选择多序列 FASTA 文件。")
            return

        try:
            aln_path, view_csv, cand_csv = run_osakt2_msa_wsl(fasta)
        except Exception as e:
            messagebox.showerror(
                "MSA 失败",
                f"在 WSL 中运行 Clustal-Omega 或后处理时出错：\n{e}",
            )
            return

        msg = (
            "多序列比对 + 自动候选筛选 已完成：\n\n"
            f"对齐文件：\n{aln_path}\n\n"
            f"可视化表：\n{view_csv}\n\n"
            f"候选位点（自动）：\n{cand_csv}\n\n"
            "用 Excel 打开候选表，挑选你感兴趣的位点，再手动整理到 candidate_sites_v0.1.xlsx 即可。"
        )
        messagebox.showinfo("MSA 完成", msg)

    tk.Button(
        msa_frame,
        text="一键跑 MSA + 自动候选",
        command=on_run_msa,
        width=20,
    ).grid(row=1, column=0, columnspan=3, padx=5, pady=(6, 0))

    # ===== HOLE 管道配置 =====
    hole_frame = tk.LabelFrame(hole_container, text="HOLE 管道配置", padx=8, pady=8)
    hole_frame.pack(fill="x", padx=10, pady=5)
    hole_frame.grid_columnconfigure(1, weight=1)
    hole_frame.grid_columnconfigure(3, weight=1)

    hole_base_dir_var = tk.StringVar(value=r"D:\demo\hole")
    tk.Label(hole_frame, text="HOLE 工作目录：").grid(row=0, column=0, sticky="w")
    tk.Entry(hole_frame, textvariable=hole_base_dir_var, width=60).grid(row=0, column=1, columnspan=2, sticky="we")

    def browse_hole_dir():
        path = filedialog.askdirectory(title="选择 HOLE 工作目录")
        if path:
            hole_base_dir_var.set(path)

    tk.Button(hole_frame, text="浏览", command=browse_hole_dir).grid(row=0, column=3, padx=5)

    hole_models_var = tk.StringVar(value="WT,DMI,DMT,GT,ND")
    tk.Label(hole_frame, text="模型列表（逗号分隔）：").grid(row=1, column=0, sticky="w", pady=(6, 0))
    tk.Entry(hole_frame, textvariable=hole_models_var, width=60).grid(row=1, column=1, columnspan=3, sticky="we", pady=(6, 0))

    hole_cpx_var = tk.StringVar(value="38.30")
    hole_cpy_var = tk.StringVar(value="7.55")
    hole_cpz_var = tk.StringVar(value="-22.15")
    hole_cvx_var = tk.StringVar(value="0.0")
    hole_cvy_var = tk.StringVar(value="0.0")
    hole_cvz_var = tk.StringVar(value="1.0")

    tk.Label(hole_frame, text="cpoint (Å)：").grid(row=2, column=0, sticky="w", pady=(6, 0))
    tk.Entry(hole_frame, textvariable=hole_cpx_var, width=8).grid(row=2, column=1, sticky="w", pady=(6, 0))
    tk.Entry(hole_frame, textvariable=hole_cpy_var, width=8).grid(row=2, column=2, sticky="w", pady=(6, 0))
    tk.Entry(hole_frame, textvariable=hole_cpz_var, width=8).grid(row=2, column=3, sticky="w", pady=(6, 0))

    tk.Label(hole_frame, text="cvect：").grid(row=3, column=0, sticky="w", pady=(6, 0))
    tk.Entry(hole_frame, textvariable=hole_cvx_var, width=8).grid(row=3, column=1, sticky="w", pady=(6, 0))
    tk.Entry(hole_frame, textvariable=hole_cvy_var, width=8).grid(row=3, column=2, sticky="w", pady=(6, 0))
    tk.Entry(hole_frame, textvariable=hole_cvz_var, width=8).grid(row=3, column=3, sticky="w", pady=(6, 0))

    hole_sample_var = tk.StringVar(value="0.25")
    hole_endrad_var = tk.StringVar(value="15.0")
    tk.Label(hole_frame, text="sample (Å)：").grid(row=4, column=0, sticky="w", pady=(6, 0))
    tk.Entry(hole_frame, textvariable=hole_sample_var, width=8).grid(row=4, column=1, sticky="w", pady=(6, 0))
    tk.Label(hole_frame, text="endrad (Å)：").grid(row=4, column=2, sticky="w", pady=(6, 0))
    tk.Entry(hole_frame, textvariable=hole_endrad_var, width=8).grid(row=4, column=3, sticky="w", pady=(6, 0))

    hole_radius_var = tk.StringVar(value="simple.rad")
    tk.Label(hole_frame, text="radius 文件名：").grid(row=5, column=0, sticky="w", pady=(6, 0))
    tk.Entry(hole_frame, textvariable=hole_radius_var, width=20).grid(row=5, column=1, sticky="w", pady=(6, 0))

    # 默认留空，让后台自动选择如何调用 HOLE
    hole_cmd_var = tk.StringVar(value="")
    tk.Label(hole_frame, text="HOLE 命令（可留空）：").grid(row=6, column=0, sticky="w", pady=(6, 0))
    tk.Entry(hole_frame, textvariable=hole_cmd_var, width=24).grid(row=6, column=1, sticky="w", pady=(6, 0))

    hole_run_var = tk.IntVar(value=0)
    hole_parse_var = tk.IntVar(value=1)
    tk.Checkbutton(
        hole_frame,
        text="自动在 WSL 中调用 HOLE",
        variable=hole_run_var,
    ).grid(row=7, column=0, columnspan=2, sticky="w", pady=(8, 0))
    tk.Checkbutton(
        hole_frame,
        text="根据 *_hole.log 生成 CSV & 曲线图",
        variable=hole_parse_var,
    ).grid(row=7, column=2, columnspan=2, sticky="w", pady=(8, 0))

    def on_prepare_hole_pdb():
        """从 WT + 突变体列表，把对应 PDB 复制到 HOLE 工作目录，命名为 <模型名>.pdb。"""
        base_dir = hole_base_dir_var.get().strip()
        if not base_dir:
            messagebox.showerror("缺少 HOLE 目录", "请先在 HOLE 区设置“HOLE 工作目录”。")
            return

        models_str = hole_models_var.get().strip()
        if not models_str:
            messagebox.showerror("缺少模型列表", "请先在 HOLE 区设置模型列表，比如 WT,DMI,DMT,GT,ND。")
            return

        models = [m.strip() for m in models_str.split(",") if m.strip()]
        if not models:
            messagebox.showerror("无模型", "模型列表似乎是空的。")
            return

        mapping = {}

        if "WT" in models:
            wt_path = wt_path_var.get().strip()
            if not wt_path:
                messagebox.showerror("缺少 WT PDB", "模型列表包含 WT，但上面没有选择 WT PDB。")
                return
            if not os.path.exists(wt_path):
                messagebox.showerror("WT PDB 不存在", f"找不到 WT PDB 文件：\n{wt_path}")
                return
            mapping["WT"] = wt_path

        label_to_pdb = {}
        for idx, row in enumerate(mutant_rows, start=1):
            label = (row["label_var"].get() or "").strip()
            pdb = (row["pdb_var"].get() or "").strip()
            if not label or not pdb:
                continue
            label_to_pdb[label] = pdb

        for m in models:
            if m == "WT":
                continue
            pdb_path = label_to_pdb.get(m)
            if not pdb_path:
                messagebox.showerror(
                    "缺少突变体 PDB",
                    f"模型列表包含 {m}，但在“突变体 PDB”区没有找到同名标签的行，"
                    "请确保突变体标签和模型名一模一样。"
                )
                return
            if not os.path.exists(pdb_path):
                messagebox.showerror("PDB 文件不存在", f"{m} 对应的 PDB 不存在：\n{pdb_path}")
                return
            mapping[m] = pdb_path

        os.makedirs(base_dir, exist_ok=True)

        for label, src in mapping.items():
            dst = os.path.join(base_dir, f"{label}.pdb")
            try:
                shutil.copy2(src, dst)
            except Exception as e:
                messagebox.showerror("复制失败", f"复制 {label} 时出错：\n{e}")
                return

        msg = (
            "已把 WT + 突变体 PDB 复制到 HOLE 工作目录：\n"
            f"{base_dir}\n\n"
            "文件名统一为：<模型名>.pdb。"
        )
        messagebox.showinfo("准备完成", msg)

    tk.Button(
        hole_frame,
        text="从 WT+突变体准备 HOLE PDB",
        command=on_prepare_hole_pdb,
    ).grid(row=8, column=0, columnspan=4, sticky="w", pady=(8, 0))

    # ===== 自动推荐 cpoint / cvect 小工具 =====
    axis_frame = tk.LabelFrame(hole_frame, text="自动推荐 cpoint / cvect", padx=8, pady=8)
    axis_frame.grid(row=9, column=0, columnspan=4, sticky="we", pady=(10, 0))

    axis_chain_var = tk.StringVar(value="A")
    axis_res_expr_var = tk.StringVar(value="")
    axis_log_var = tk.StringVar(value="axis_axis.log")

    tk.Label(axis_frame, text="链 ID：").grid(row=0, column=0, sticky="w")
    tk.Entry(axis_frame, textvariable=axis_chain_var, width=5).grid(row=0, column=1, sticky="w")

    tk.Label(axis_frame, text="找轴残基表达式：").grid(row=0, column=2, sticky="w")
    tk.Entry(axis_frame, textvariable=axis_res_expr_var, width=20).grid(row=0, column=3, sticky="w")

    tk.Label(axis_frame, text="axis log 文件名：").grid(row=1, column=0, sticky="w", pady=(4, 0))
    tk.Entry(axis_frame, textvariable=axis_log_var, width=20).grid(row=1, column=1, sticky="w", pady=(4, 0))

    def on_generate_axis_cxc():
        wt_pdb = wt_path_var.get().strip()
        base_dir = hole_base_dir_var.get().strip()
        chain = axis_chain_var.get().strip() or "A"
        res_expr = axis_res_expr_var.get().strip()

        if not wt_pdb:
            messagebox.showerror("缺少 WT", "请先在上面选择 WT PDB。")
            return
        if not base_dir:
            messagebox.showerror("缺少目录", "请先设置 HOLE 工作目录。")
            return
        if not res_expr:
            messagebox.showerror("缺少残基", "请填写用于找轴的残基表达式（例如 298-300）。")
            return

        try:
            axis_cxc = build_axis_cxc(
                wt_pdb_path=wt_pdb,
                chain_id=chain,
                residue_expr=res_expr,
                out_dir=base_dir,
                label="axis",
            )
        except Exception as e:
            messagebox.showerror("生成失败", f"创建找轴脚本时出错：\n{e}")
            return

        msg = (
            "脚本路径：\n"
            f"{axis_cxc}\n\n"
            "接下来在 ChimeraX 中运行：\n"
            f"runscript {axis_cxc}\n\n"
            "跑完后会在 HOLE 目录生成 axis_axis.log。"
        )
        messagebox.showinfo("已生成找轴脚本", msg)

    def on_fill_cpoint_from_axis_log():
        base_dir = hole_base_dir_var.get().strip()
        log_name = axis_log_var.get().strip() or "axis_axis.log"
        if not base_dir:
            messagebox.showerror("缺少目录", "请先设置 HOLE 工作目录。")
            return

        log_path = os.path.join(base_dir, log_name)
        if not os.path.exists(log_path):
            messagebox.showerror("找不到 log", f"没有发现日志文件：\n{log_path}")
            return

        try:
            x, y, z = parse_axis_log(log_path)
        except Exception as e:
            messagebox.showerror("解析失败", f"解析质心坐标时出错：\n{e}")
            return

        hole_cpx_var.set(f"{x:.2f}")
        hole_cpy_var.set(f"{y:.2f}")
        hole_cpz_var.set(f"{z:.2f}")

        hole_cvx_var.set("0.0")
        hole_cvy_var.set("0.0")
        hole_cvz_var.set("1.0")

        msg = (
            f"推荐 cpoint = ({x:.2f}, {y:.2f}, {z:.2f})\n"
            "cvect 已设为 (0.0, 0.0, 1.0)。"
        )
        messagebox.showinfo("已填入推荐轴", msg)

    tk.Button(axis_frame, text="生成找轴 .cxc", command=on_generate_axis_cxc).grid(
        row=2, column=0, columnspan=2, sticky="w", pady=(6, 0)
    )
    tk.Button(
        axis_frame,
        text="从 axis log 填入 cpoint/cvect",
        command=on_fill_cpoint_from_axis_log,
    ).grid(row=2, column=2, columnspan=2, sticky="w", pady=(6, 0))

    # ===== 生成按钮 =====
    def on_generate():
        wt_pdb = wt_path_var.get().strip()
        if not wt_pdb:
            messagebox.showerror("缺少 WT", "请先选择 WT 的 PDB 文件。")
            return

        mode = mode_var.get()

        if mode == "mutate":
            mut_out_dir = mut_out_dir_var.get().strip()
            if not mut_out_dir:
                messagebox.showerror("缺少输出目录", "请选择突变体输出目录。")
                return

            mutations = []
            for idx, row in enumerate(mutation_rows, start=1):
                residue = row["residue_var"].get().strip()
                new_aa = row["new_aa_var"].get().strip().upper()
                if not residue and not new_aa:
                    continue
                if not residue or not new_aa:
                    messagebox.showerror("缺少参数", "突变行需要同时填写残基号和要改成的氨基酸。")
                    return
                label = row["label_var"].get().strip() or f"MUT{idx}"
                chain = row["chain_var"].get().strip() or "A"
                mutations.append({
                    "label": label,
                    "chain": chain,
                    "residue": residue,
                    "new_aa": new_aa,
                })

            if not mutations:
                messagebox.showerror("没有突变", "请至少填写一行突变信息。")
                return

            generated_paths = []
            try:
                os.makedirs(mut_out_dir, exist_ok=True)
                for mut in mutations:
                    cxc_path = build_mutation_cxc(
                        wt_pdb_path=wt_pdb,
                        mut_label=mut["label"],
                        chain=mut["chain"],
                        residue=mut["residue"],
                        new_aa=mut["new_aa"],
                        out_dir=mut_out_dir,
                    )
                    generated_paths.append(cxc_path)
            except Exception as e:
                messagebox.showerror("生成失败", f"创建突变脚本时出错：\n{e}")
                return

            preview = "\n".join(generated_paths)
            msg = (
                f"已生成 {len(generated_paths)} 个 swapaa 脚本：\n{preview}\n\n"
                "在 ChimeraX 中执行：\n"
                f"runscript {generated_paths[0]}"
            )
            messagebox.showinfo("完成", msg)
            return

        if mode == "hole":
            base_dir = hole_base_dir_var.get().strip()
            if not base_dir:
                messagebox.showerror("缺少目录", "请选择 HOLE 工作目录。")
                return

            models = [m.strip() for m in hole_models_var.get().split(",") if m.strip()]
            if not models:
                messagebox.showerror("缺少模型", "请至少填写一个模型名。")
                return

            try:
                cpoint = (
                    float(hole_cpx_var.get()),
                    float(hole_cpy_var.get()),
                    float(hole_cpz_var.get()),
                )
                cvect = (
                    float(hole_cvx_var.get()),
                    float(hole_cvy_var.get()),
                    float(hole_cvz_var.get()),
                )
                sample = float(hole_sample_var.get())
                endrad = float(hole_endrad_var.get())
            except ValueError:
                messagebox.showerror("参数错误", "cpoint/cvect/sample/endrad 请输入数字。")
                return

            radius_file = hole_radius_var.get().strip() or "simple.rad"
            hole_cmd = hole_cmd_var.get().strip()
            run_flag = bool(hole_run_var.get())
            parse_flag = bool(hole_parse_var.get())

            for m in models:
                model_dir = os.path.join(base_dir, f"{m}-HOLE")
                os.makedirs(model_dir, exist_ok=True)

                src_pdb = os.path.join(base_dir, f"{m}.pdb")
                if not os.path.exists(src_pdb):
                    messagebox.showerror(
                        "缺少 PDB",
                        f"找不到 {m} 的 PDB 文件：\n{src_pdb}\n请先点击“从 WT+突变体准备 HOLE PDB”。",
                    )
                    return

                dst_pdb = os.path.join(model_dir, f"{m}.pdb")
                try:
                    shutil.copy2(src_pdb, dst_pdb)
                except Exception as e:
                    messagebox.showerror("复制失败", f"复制 {m} 的 PDB 到子目录时出错：\n{e}")
                    return

                radius_src = os.path.join(base_dir, radius_file)
                if os.path.isfile(radius_src):
                    try:
                        shutil.copy2(radius_src, os.path.join(model_dir, radius_file))
                    except Exception:
                        pass

                try:
                    hole_write_input(
                        base_dir_win=model_dir,
                        model=m,
                        cpoint=cpoint,
                        cvect=cvect,
                        sample=sample,
                        endrad=endrad,
                        radius_filename=radius_file,
                    )
                except Exception as e:
                    messagebox.showerror("生成失败", f"写入 {m} 的 HOLE 输入时出错：\n{e}")
                    return

            if run_flag:
                for m in models:
                    try:
                        hole_run_in_wsl(
                            base_dir_win=os.path.join(base_dir, f"{m}-HOLE"),
                            model=m,
                            hole_cmd=hole_cmd,
                        )
                    except Exception as e:
                        messagebox.showerror("HOLE 运行失败", f"{m} 出错：{e}")
                        return

            summary_msg = []
            if parse_flag:
                logs = {}
                for m in models:
                    log_path = os.path.join(base_dir, f"{m}-HOLE", f"{m}_hole.log")
                    if os.path.exists(log_path):
                        logs[m] = log_path
                    else:
                        summary_msg.append(f"警告：没找到 {m}_hole.log，跳过。")

                if logs:
                    try:
                        hole_summarize_logs(logs, base_dir)
                        fig_path = hole_plot_profiles(logs, base_dir)
                        summary_msg.append(
                            "已生成：hole_profile_samples.csv, hole_min_table.csv, "
                            "hole_min_summary.csv, "
                            f"{os.path.basename(fig_path)}"
                        )
                    except Exception as e:
                        messagebox.showerror("后处理失败", f"解析 log 或绘图时出错：\n{e}")
                        return

            msg_lines = [
                "已为下列模型生成 HOLE 输入文件：",
                ", ".join(models),
            ]
            if run_flag:
                msg_lines.append("已尝试在 WSL 中运行 HOLE。")
            if summary_msg:
                msg_lines.append("")
                msg_lines.extend(summary_msg)

            messagebox.showinfo("HOLE 管道完成", "\n".join(msg_lines))
            return

        # ===== 研究模式 =====
        mutants = []
        for idx, row in enumerate(mutant_rows, start=1):
            label = row["label_var"].get().strip() or f"MUT{idx}"
            pdb = row["pdb_var"].get().strip()
            if not pdb:
                continue  # 路径空就忽略这一行
            mutants.append({"label": label, "pdb": pdb})

        features = {
            "full_coulombic": bool(full_var.get()),
            "contacts": bool(contacts_var.get()),
            "hbonds": bool(hbonds_var.get()),
            "sasa": bool(sasa_var.get()),
            "sites_contacts": bool(sites_contacts_var.get()),
            "sites_coulombic": bool(sites_coulombic_var.get()),
        }
        if not any(features.values()):
            messagebox.showerror("没选功能", "至少勾选一个要自动执行的功能。")
            return

        chain_id = chain_var.get().strip() or "A"
        residue_expr = residue_expr_var.get().strip()
        roi_expr = roi_expr_var.get().strip()

        need_residue = any(features[k] for k in ("contacts", "hbonds", "sasa"))
        if need_residue and not residue_expr:
            messagebox.showerror(
                "缺少残基表达式",
                "你勾选了 2/3/4 中的至少一项，必须填写目标残基表达式。"
            )
            return

        need_roi = any(features[k] for k in ("sites_contacts", "sites_coulombic"))
        if need_roi and not roi_expr:
            messagebox.showerror("缺少 ROI", "你勾选了 5/6，必须填写 ROI 残基表达式或导入位点表。")
            return

        out_dir = out_dir_var.get().strip()
        if not out_dir:
            messagebox.showerror("缺少输出目录", "请指定图片 / 文本的输出目录。")
            return

        cxc_path = cxc_path_var.get().strip()
        if not cxc_path:
            messagebox.showerror("缺少 .cxc 路径", "请指定要保存的 .cxc 脚本文件名。")
            return

        try:
            script_text = build_cxc_script(
                wt_pdb_path=wt_pdb,
                mutant_list=mutants,
                chain_id=chain_id,
                residue_expr=residue_expr,
                out_dir=out_dir,
                features=features,
                roi_expr=roi_expr,
            )
        except Exception as e:
            messagebox.showerror("生成失败", f"生成脚本时出错：\n{e}")
            return

        try:
            os.makedirs(os.path.dirname(cxc_path), exist_ok=True)
            with open(cxc_path, "w", encoding="utf-8") as f:
                f.write(script_text)
        except OSError as e:
            messagebox.showerror("写文件失败", f"无法写入 .cxc 文件：\n{e}")
            return

        msg = (
            f"已生成 ChimeraX 脚本：\n{cxc_path}\n\n"
            "在 ChimeraX 命令行中执行：\n"
            f"runscript {cxc_path}"
        )
        messagebox.showinfo("完成", msg)

    def on_plot_hole_metrics():
        base_dir = hole_base_dir_var.get().strip()
        if not base_dir:
            messagebox.showerror("缺少目录", "请选择 HOLE 工作目录。")
            return

        try:
            paths = plot_basic_hole_metrics(base_dir)
        except Exception as e:
            messagebox.showerror("绘图失败", f"生成 HOLE 对比图时出错：\n{e}")
            return

        messagebox.showinfo(
            "绘图完成",
            "已在 HOLE 工作目录生成基础对比图：\n" + "\n".join(paths)
        )

    btn_frame = tk.Frame(research_container)
    btn_frame.pack(fill="x", padx=10, pady=10)

    tk.Button(btn_frame, text="生成 .cxc", command=on_generate, width=15).pack(side="left")

    fit_enabled_var = tk.BooleanVar(value=True)
    standard_rows = []

    fit_frame = tk.LabelFrame(research_container, text="尺子拟合（可选）")
    fit_frame.pack(fill="x", padx=10, pady=(6, 8))

    tk.Checkbutton(
        fit_frame,
        text="启用标准集拟合（不勾选则默认权重）",
        variable=fit_enabled_var,
    ).pack(anchor="w")

    rows_container = tk.Frame(fit_frame)
    rows_container.pack(fill="x", pady=4)

    def collect_models_from_gui():
        models = ["WT"]
        for row in mutant_rows:
            label = (row["label_var"].get() or "").strip()
            if label:
                models.append(label)
        seen = set()
        uniq = []
        for m in models:
            if m not in seen:
                seen.add(m)
                uniq.append(m)
        return uniq

    def rebuild_standards_rows():
        for w in rows_container.winfo_children():
            w.destroy()
        standard_rows.clear()

        models = collect_models_from_gui()

        tk.Label(rows_container, text="作为标准").grid(row=0, column=0, padx=6, sticky="w")
        tk.Label(rows_container, text="Model").grid(row=0, column=1, padx=6, sticky="w")
        tk.Label(
            rows_container,
            text="y（目标分，例：better=1/similar=0/worse=-1）",
        ).grid(row=0, column=2, padx=6, sticky="w")

        for i, m in enumerate(models, start=1):
            use_var = tk.BooleanVar(value=(m == "WT"))
            y_var = tk.StringVar(value=("0" if m == "WT" else ""))

            tk.Checkbutton(rows_container, variable=use_var).grid(
                row=i, column=0, padx=6, sticky="w"
            )
            tk.Label(rows_container, text=m).grid(row=i, column=1, padx=6, sticky="w")
            tk.Entry(rows_container, textvariable=y_var, width=18).grid(
                row=i, column=2, padx=6, sticky="w"
            )

            standard_rows.append({"model": m, "use_var": use_var, "y_var": y_var})

    tk.Button(fit_frame, text="刷新模型列表", command=rebuild_standards_rows).pack(
        anchor="w", pady=(4, 0)
    )

    rebuild_standards_rows()

    def on_summarize_sasa_hbonds():
        out_dir = out_dir_var.get().strip()
        if not out_dir:
            messagebox.showerror("缺少输出目录", "请先在上面设置“图片 / 文本输出目录”。")
            return

        try:
            summary_csv, detail_csv = summarize_sasa_hbonds(out_dir)
        except Exception as e:
            messagebox.showerror("汇总失败", f"解析 SASA/H-bonds 日志时出错：\n{e}")
            return

        labels = ["WT"]
        for idx, row in enumerate(mutant_rows, start=1):
            label = (row["label_var"].get() or "").strip()
            if label:
                labels.append(label)
        labels = sorted(set(labels))

        for label in labels:
            subdir = os.path.join(out_dir, label)
            os.makedirs(subdir, exist_ok=True)
            prefix = f"{label}_"
            for name in os.listdir(out_dir):
                if not name.startswith(prefix):
                    continue
                src = os.path.join(out_dir, name)
                if not os.path.isfile(src):
                    continue
                new_name = name[len(prefix):] if len(name) > len(prefix) else name
                dst = os.path.join(subdir, new_name)
                try:
                    shutil.move(src, dst)
                except Exception:
                    continue

        msg = (
            "已生成：\n"
            f"{summary_csv}\n"
            f"{detail_csv}\n\n"
            "可以用 Excel 打开做对比分析。"
        )
        messagebox.showinfo("汇总完成", msg)

    def build_standards_csv_from_gui(sasa_dir: str) -> str | None:
        if not fit_enabled_var.get():
            return None

        points = []
        for r in standard_rows:
            if not r["use_var"].get():
                continue
            raw = (r["y_var"].get() or "").strip()
            if raw == "":
                continue
            try:
                y = float(raw)
            except Exception:
                continue
            points.append((r["model"], y))

        if len(points) < 3:
            messagebox.showwarning("标准点不足", "用于拟合的标准点少于 3 个，本次将使用默认权重评分。")
            return None

        os.makedirs(sasa_dir, exist_ok=True)
        std_path = os.path.join(sasa_dir, "standards_gui.csv")
        with open(std_path, "w", newline="", encoding="utf-8-sig") as f:
            writer = csv.writer(f)
            writer.writerow(["Model", "y"])
            for model, y in points:
                writer.writerow([model, y])
        return std_path


    def on_merge_all_metrics():
        hole_dir = hole_base_dir_var.get().strip()
        sasa_dir = out_dir_var.get().strip()

        if not hole_dir:
            messagebox.showerror("缺少 HOLE 目录", "请先在 HOLE 模式里设置 HOLE 工作目录。")
            return
        if not sasa_dir:
            messagebox.showerror("缺少 SASA 目录", "请先在研究模式里设置“图片 / 文本输出目录”。")
            return

        metrics_all = os.path.join(sasa_dir, "metrics_all.csv")

        std_csv = build_standards_csv_from_gui(sasa_dir)

        try:
            merge_all_metrics(hole_dir=hole_dir, sasa_dir=sasa_dir, out_csv=metrics_all)
        except Exception as e:
            messagebox.showerror("合并失败", f"merge_all_metrics 出错：\n{e}")
            return

        std_used = std_csv
        try:
            scored_path = score_metrics_file(
                metrics_all,
                wt_name="WT",
                pdb_dir=hole_dir,
                standards_csv=std_used,
            )
        except Exception as e:
            if std_used:
                messagebox.showwarning(
                    "拟合失败，已回退默认权重",
                    f"标准集：{std_used}\n\n原因：{e}\n\n将使用默认权重继续评分。",
                )
                try:
                    scored_path = score_metrics_file(
                        metrics_all,
                        wt_name="WT",
                        pdb_dir=hole_dir,
                        standards_csv=None,
                    )
                    std_used = None
                except Exception as e2:
                    messagebox.showerror("评分失败", f"回退默认权重仍失败：\n{e2}")
                    return
            else:
                messagebox.showerror("评分失败", f"score_metrics 出错：\n{e}")
                return
        used = "GUI标准集（standards_gui.csv）" if std_used else "默认权重"

        msg = (
            "已生成 HOLE + SASA 总表：\n"
            f"{metrics_all}\n\n"
            "并已写出结构评分 + 置信度表：\n"
            f"{scored_path}\n\n"
            "metrics_scored.csv 里包含：\n"
            "- r_min / gate_length / HBonds / SASA_residue\n"
            "- GateTightScore / TotalScore / ScoreClass\n"
            "- pLDDT 均值 / 中位数 / 低-中-高残基数 + ConfidenceClass\n\n"
            f"本次尺子来源：{used}"
        )
        messagebox.showinfo("合并 + 评分完成", msg)

    def on_make_stage3_table():
        sasa_dir = out_dir_var.get().strip()

        if not sasa_dir:
            messagebox.showerror("缺少输出目录", "请先在研究模式里设置“图片 / 文本输出目录”。")
            return

        try:
            table_path = make_stage3_table(sasa_dir)
        except Exception as e:
            messagebox.showerror("生成失败", f"make_stage3_table 出错：\n{e}")
            return

        msg = (
            "已生成决策表模板：\n"
            f"{table_path}\n\n"
            "Patch_Electrostatics / Contacts_Qualitative 两列留空，方便看图填写。"
        )
        messagebox.showinfo("决策表已生成", msg)

    tk.Button(
            btn_frame,
            text="汇总 SASA / H-bonds",
            command=on_summarize_sasa_hbonds,
            width=18,
    ).pack(side="left", padx=8)

    tk.Button(
            btn_frame,
            text="合并 HOLE + SASA 指标",
        command=on_merge_all_metrics,
        width=22,
    ).pack(side="left", padx=8)

    tk.Button(
        btn_frame,
        text="生成决策表模板",
        command=on_make_stage3_table,
        width=22,
    ).pack(side="left", padx=8)

    tk.Label(
        btn_frame,
        text="先在 ChimeraX 里 runscript 跑完，再点汇总，就会吐出 CSV 指标表。",
        fg="#555"
    ).pack(side="left", padx=10)
    tk.Label(
        btn_frame,
        text="窗口不会自动退出。",
        fg="#555"
    ).pack(side="left", padx=10)

    # ===== 突变模式：突变体构建 =====
    mut_builder_frame = tk.LabelFrame(mutate_container, text="突变体构建（swapaa）", padx=8, pady=8)
    mut_builder_frame.pack(fill="x", padx=10, pady=5)

    tk.Label(
        mut_builder_frame,
        text="一次配置多个突变，生成对应的 swapaa .cxc（输出在所选目录/MUT/ 下）。",
        fg="#555"
    ).grid(row=0, column=0, columnspan=6, sticky="w")
    tk.Label(
        mut_builder_frame,
        text="提示：链 ID / 残基号 / 目标氨基酸都可以用逗号或空格分隔多个（数量要匹配）。",
        fg="#777"
    ).grid(row=1, column=0, columnspan=6, sticky="w", pady=(4, 0))

    mutation_rows = []
    mutations_scroll = ScrollableFrame(mut_builder_frame, height=200)
    mutations_scroll.grid(row=2, column=0, columnspan=6, sticky="nsew", pady=(8, 4))
    mutations_inner = mutations_scroll.scrollable_frame
    mut_builder_frame.grid_rowconfigure(2, weight=1)
    mut_builder_frame.grid_columnconfigure(0, weight=1)

    def add_mutation_row(default_label=None):
        idx = len(mutation_rows) + 1
        row = tk.Frame(mutations_inner)
        row.pack(fill="x", pady=2)

        label_var = tk.StringVar(value=default_label or f"MUT{idx}")
        chain_var = tk.StringVar(value="A")
        residue_var = tk.StringVar()
        new_aa_var = tk.StringVar()

        tk.Label(row, text="标签：").grid(row=0, column=0, sticky="w")
        tk.Entry(row, textvariable=label_var, width=10).grid(row=0, column=1, sticky="w", padx=(0, 10))

        tk.Label(row, text="链(可多个)：").grid(row=0, column=2, sticky="w")
        tk.Entry(row, textvariable=chain_var, width=5).grid(row=0, column=3, sticky="w", padx=(0, 10))

        tk.Label(row, text="残基号(可多个)：").grid(row=0, column=4, sticky="w")
        tk.Entry(row, textvariable=residue_var, width=8).grid(row=0, column=5, sticky="w", padx=(0, 10))

        tk.Label(row, text="改成(可多个)：").grid(row=0, column=6, sticky="w")
        tk.Entry(row, textvariable=new_aa_var, width=12).grid(row=0, column=7, sticky="w", padx=(0, 10))

        def delete_row():
            if row_dict in mutation_rows:
                mutation_rows.remove(row_dict)
            row.destroy()

        del_btn = tk.Button(row, text="删除", fg="red", command=delete_row)
        del_btn.grid(row=0, column=8, padx=(0, 5))

        row_dict = {
            "label_var": label_var,
            "chain_var": chain_var,
            "residue_var": residue_var,
            "new_aa_var": new_aa_var,
            "frame": row,
        }
        mutation_rows.append(row_dict)

    add_mutation_row()

    tk.Button(mut_builder_frame, text="添加突变", command=add_mutation_row).grid(
        row=3, column=0, columnspan=2, sticky="w", pady=(4, 0)
    )

    mut_out_dir_var = tk.StringVar(value=os.path.join("D:\\", "demo"))
    tk.Label(mut_builder_frame, text="输出目录：").grid(row=4, column=0, sticky="w", pady=(10, 0))
    tk.Entry(mut_builder_frame, textvariable=mut_out_dir_var, width=50).grid(
        row=4, column=1, columnspan=4, sticky="w", pady=(10, 0)
    )

    def browse_mut_out_dir():
        path = filedialog.askdirectory(title="选择输出目录")
        if path:
            mut_out_dir_var.set(path)

    tk.Button(mut_builder_frame, text="浏览", command=browse_mut_out_dir).grid(
        row=4, column=5, padx=5, pady=(10, 0)
    )

    mut_btn_frame = tk.Frame(mutate_container)
    mut_btn_frame.pack(fill="x", padx=10, pady=10)
    tk.Button(mut_btn_frame, text="生成突变脚本", command=on_generate, width=18).pack(side="left")
    tk.Label(
        mut_btn_frame,
        text="勾选突变模式后，只会生成 swapaa 的 .cxc 文件。",
        fg="#555"
    ).pack(side="left", padx=10)

    hole_btn_frame = tk.Frame(hole_container)
    hole_btn_frame.pack(fill="x", padx=10, pady=10)
    tk.Button(hole_btn_frame, text="执行 HOLE 管道", command=on_generate, width=18).pack(side="left")
    tk.Button(
        hole_btn_frame,
        text="画 HOLE 对比图",
        command=on_plot_hole_metrics,
        width=18,
    ).pack(side="left", padx=8)

    # 默认展示研究模式
    update_mode()

    return root

class ScrollableFrame(tk.Frame):
    def __init__(self, container, *args, height=280, **kwargs):
        super().__init__(container, *args, **kwargs)

        canvas = tk.Canvas(self, height=height)
        scrollbar = tk.Scrollbar(self, orient="vertical", command=canvas.yview)
        self.scrollable_frame = tk.Frame(canvas)

        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(
                scrollregion=canvas.bbox("all")
            )
        )

        canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")



if __name__ == "__main__":
    app = create_gui()
    app.mainloop()
