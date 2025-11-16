import os
import shutil
import tkinter as tk
from tkinter import filedialog, messagebox

from PP import build_cxc_script, build_mutation_cxc


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

    mode_var = tk.StringVar(value="research")

    def update_mode():
        if mode_var.get() == "research":
            mutate_container.pack_forget()
            research_container.pack(fill="both", expand=True)
        else:
            research_container.pack_forget()
            mutate_container.pack(fill="both", expand=True)

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

    tk.Button(mutants_frame, text="添加突变体", command=add_mutant_row).pack(anchor="w", pady=4)

    # ===== 功能选择 =====
    feature_frame = tk.LabelFrame(research_container, text="要让 ChimeraX 自动干的事", padx=8, pady=8)
    feature_frame.pack(fill="x", padx=10, pady=5)

    full_var = tk.IntVar(value=1)
    contacts_var = tk.IntVar(value=1)
    hbonds_var = tk.IntVar(value=1)
    sasa_var = tk.IntVar(value=1)

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

    tk.Label(
        feature_frame,
        text="说明：2 / 3 / 4 需要你指定“目标残基”，只勾 1 的话可以不填残基。",
        fg="#555"
    ).grid(row=2, column=0, columnspan=2, sticky="w", pady=(4, 0))

    # ===== 目标残基 & 链 =====
    target_frame = tk.LabelFrame(research_container, text="目标残基（可选，不一定是三联体）", padx=8, pady=8)
    target_frame.pack(fill="x", padx=10, pady=5)

    chain_var = tk.StringVar(value="A")
    residue_expr_var = tk.StringVar()

    tk.Label(target_frame, text="链 ID：").grid(row=0, column=0, sticky="w")
    tk.Entry(target_frame, textvariable=chain_var, width=5).grid(row=0, column=1, sticky="w")

    tk.Label(target_frame, text="残基表达式：").grid(row=0, column=2, sticky="w", padx=(15, 0))
    tk.Entry(target_frame, textvariable=residue_expr_var, width=40).grid(row=0, column=3, sticky="w")

    tk.Label(
        target_frame,
        text="例：298,299,300 或 298-305，后续可以自己扩展更复杂写法。",
        fg="#555"
    ).grid(row=1, column=0, columnspan=4, sticky="w", pady=(4, 0))

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

    # ===== 生成按钮 =====
    def on_generate():
        wt_pdb = wt_path_var.get().strip()
        if not wt_pdb:
            messagebox.showerror("缺少 WT", "请先选择 WT 的 PDB 文件。")
            return

        if mode_var.get() == "mutate":
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
            messagebox.showinfo(
                "完成",
                f"已生成 {len(generated_paths)} 个 swapaa 脚本：\n{preview}\n\n"
                "在 ChimeraX 中执行：\n"
                f"runscript {generated_paths[0]}"
            )
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
        }
        if not any(features.values()):
            messagebox.showerror("没选功能", "至少勾选一个要自动执行的功能。")
            return

        chain_id = chain_var.get().strip() or "A"
        residue_expr = residue_expr_var.get().strip()

        need_residue = any(features[k] for k in ("contacts", "hbonds", "sasa"))
        if need_residue and not residue_expr:
            messagebox.showerror(
                "缺少残基表达式",
                "你勾选了 2/3/4 中的至少一项，必须填写目标残基表达式。"
            )
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

        messagebox.showinfo(
            "完成",
            f"已生成 ChimeraX 脚本：\n{cxc_path}\n\n"
            "在 ChimeraX 命令行中执行：\n"
            f"runscript {cxc_path}"
        )

    btn_frame = tk.Frame(research_container)
    btn_frame.pack(fill="x", padx=10, pady=10)

    tk.Button(btn_frame, text="生成 .cxc", command=on_generate, width=15).pack(side="left")
    tk.Label(
        btn_frame,
        text="窗口不会自动退出，你可以多次修改参数反复生成不同的 .cxc。",
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

    mutation_rows = []
    mutations_scroll = ScrollableFrame(mut_builder_frame, height=200)
    mutations_scroll.grid(row=1, column=0, columnspan=6, sticky="nsew", pady=(8, 4))
    mutations_inner = mutations_scroll.scrollable_frame
    mut_builder_frame.grid_rowconfigure(1, weight=1)
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

        tk.Label(row, text="链：").grid(row=0, column=2, sticky="w")
        tk.Entry(row, textvariable=chain_var, width=5).grid(row=0, column=3, sticky="w", padx=(0, 10))

        tk.Label(row, text="残基号：").grid(row=0, column=4, sticky="w")
        tk.Entry(row, textvariable=residue_var, width=8).grid(row=0, column=5, sticky="w", padx=(0, 10))

        tk.Label(row, text="改成：").grid(row=0, column=6, sticky="w")
        tk.Entry(row, textvariable=new_aa_var, width=8).grid(row=0, column=7, sticky="w", padx=(0, 10))

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
        row=2, column=0, columnspan=2, sticky="w", pady=(4, 0)
    )

    mut_out_dir_var = tk.StringVar(value=os.path.join("D:\\", "demo"))
    tk.Label(mut_builder_frame, text="输出目录：").grid(row=3, column=0, sticky="w", pady=(10, 0))
    tk.Entry(mut_builder_frame, textvariable=mut_out_dir_var, width=50).grid(
        row=3, column=1, columnspan=4, sticky="w", pady=(10, 0)
    )

    def browse_mut_out_dir():
        path = filedialog.askdirectory(title="选择输出目录")
        if path:
            mut_out_dir_var.set(path)

    tk.Button(mut_builder_frame, text="浏览", command=browse_mut_out_dir).grid(
        row=3, column=5, padx=5, pady=(10, 0)
    )

    mut_btn_frame = tk.Frame(mutate_container)
    mut_btn_frame.pack(fill="x", padx=10, pady=10)
    tk.Button(mut_btn_frame, text="生成突变脚本", command=on_generate, width=18).pack(side="left")
    tk.Label(
        mut_btn_frame,
        text="勾选突变模式后，只会生成 swapaa 的 .cxc 文件。",
        fg="#555"
    ).pack(side="left", padx=10)

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
