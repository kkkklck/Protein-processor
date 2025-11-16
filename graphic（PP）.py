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

    # ===== WT 区 =====
    wt_frame = tk.LabelFrame(root, text="WT PDB", padx=8, pady=8)
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

    # ===== 突变体区 =====
    mutants_outer = tk.LabelFrame(root, text="突变体 PDB（可选，多条）", padx=8, pady=8)
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
    feature_frame = tk.LabelFrame(root, text="要让 ChimeraX 自动干的事", padx=8, pady=8)
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
    target_frame = tk.LabelFrame(root, text="目标残基（可选，不一定是三联体）", padx=8, pady=8)
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

    # ===== 突变体构建 (swapaa) =====
    mut_builder_frame = tk.LabelFrame(root, text="突变体构建（swapaa）", padx=8, pady=8)
    mut_builder_frame.pack(fill="x", padx=10, pady=5)

    gen_mut_var = tk.IntVar(value=0)
    mut_chain_var = tk.StringVar(value="A")
    mut_residue_var = tk.StringVar()
    mut_new_aa_var = tk.StringVar(value="ASP")

    tk.Checkbutton(
        mut_builder_frame,
        text="生成突变体构建 cxc（swapaa）",
        variable=gen_mut_var
    ).grid(row=0, column=0, columnspan=4, sticky="w")

    tk.Label(mut_builder_frame, text="链：").grid(row=1, column=0, sticky="w", pady=(6, 0))
    tk.Entry(mut_builder_frame, textvariable=mut_chain_var, width=5).grid(row=1, column=1, sticky="w", pady=(6, 0))

    tk.Label(mut_builder_frame, text="残基号：").grid(row=1, column=2, sticky="w", padx=(10, 0), pady=(6, 0))
    tk.Entry(mut_builder_frame, textvariable=mut_residue_var, width=10).grid(row=1, column=3, sticky="w", pady=(6, 0))

    tk.Label(mut_builder_frame, text="改成：").grid(row=2, column=0, sticky="w", pady=(6, 0))
    tk.Entry(mut_builder_frame, textvariable=mut_new_aa_var, width=8).grid(row=2, column=1, sticky="w", pady=(6, 0))

    tk.Label(
        mut_builder_frame,
        text="会把 WT 复制到输出目录/MUT/ 下，并生成同名 .cxc 用于 ChimeraX swapaa。",
        fg="#555"
    ).grid(row=3, column=0, columnspan=4, sticky="w", pady=(4, 0))

    # ===== 输出设置 =====
    out_frame = tk.LabelFrame(root, text="输出位置", padx=8, pady=8)
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

        # 收集突变体
        mutants = []
        for idx, row in enumerate(mutant_rows, start=1):
            label = row["label_var"].get().strip() or f"MUT{idx}"
            pdb = row["pdb_var"].get().strip()
            if not pdb:
                continue  # 路径空就忽略这一行
            mutants.append({"label": label, "pdb": pdb})

        # 至少得有一个功能被勾选
        features = {
            "full_coulombic": bool(full_var.get()),
            "contacts": bool(contacts_var.get()),
            "hbonds": bool(hbonds_var.get()),
            "sasa": bool(sasa_var.get()),
        }
        if not any(features.values()):
            messagebox.showerror("没选功能", "至少勾选一个要自动执行的功能。")
            return

        # 目标残基（只在需要时强制）
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

        # 尝试生成脚本
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

        mut_builder_outputs = []
        if gen_mut_var.get():
            if not mutant_rows:
                messagebox.showerror("没有突变体", "请至少保留一行突变体标签，用于生成 swapaa 脚本。")
                return

            mut_chain = mut_chain_var.get().strip() or "A"
            mut_residue = mut_residue_var.get().strip()
            mut_new_aa = mut_new_aa_var.get().strip().upper()

            if not mut_residue:
                messagebox.showerror("缺少残基号", "请填写要突变的残基编号。")
                return
            if not mut_new_aa:
                messagebox.showerror("缺少新氨基酸", "请填写要替换成的氨基酸三字母代码。")
                return

            mut_dir = os.path.join(out_dir, "MUT")

            try:
                os.makedirs(mut_dir, exist_ok=True)
                for idx, row in enumerate(mutant_rows, start=1):
                    mut_label = row["label_var"].get().strip() or f"MUT{idx}"
                    wt_copy_path = os.path.join(mut_dir, f"{mut_label}.pdb")
                    shutil.copy2(wt_pdb, wt_copy_path)
                    cxc_file = build_mutation_cxc(
                        wt_pdb,
                        mut_label=mut_label,
                        chain=mut_chain,
                        residue=mut_residue,
                        new_aa=mut_new_aa,
                        out_dir=out_dir,
                    )
                    mut_builder_outputs.append(f"- {mut_label}: {cxc_file}")
            except Exception as e:
                messagebox.showerror("突变体构建失败", f"生成突变体 cxc 时出错：\n{e}")
                return

        msg = (
            f"已生成 ChimeraX 脚本：\n{cxc_path}\n\n"
            f"在 ChimeraX 命令行中执行：\nrunscript {cxc_path}"
        )
        if mut_builder_outputs:
            msg += "\n\n突变体构建 cxc：\n" + "\n".join(mut_builder_outputs)

        messagebox.showinfo("完成", msg)

    btn_frame = tk.Frame(root)
    btn_frame.pack(fill="x", padx=10, pady=10)

    tk.Button(btn_frame, text="生成 .cxc", command=on_generate, width=15).pack(side="left")
    tk.Label(
        btn_frame,
        text="窗口不会自动退出，你可以多次修改参数反复生成不同的 .cxc。",
        fg="#555"
    ).pack(side="left", padx=10)

    return root

class ScrollableFrame(tk.Frame):
    def __init__(self, container, *args, **kwargs):
        super().__init__(container, *args, **kwargs)

        canvas = tk.Canvas(self, height=280)
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
