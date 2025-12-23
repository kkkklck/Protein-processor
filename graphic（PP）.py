A. 一键“安全清理”（主按钮）

目标：不碰关键成果，只扫“明显没用/过期”的东西。
做法：复用你 delet.py 的扫描/执行内核，但加上 PP 的“保护规则”：

默认保护（永不清）

tables/、gate_sites/

*_tables/、*_gate_sites/（你 PP 里已经把这俩当作 minimal 的核心保留项）

PP

汇总结果：比如 sasa_hbonds_summary.csv、sasa_per_residue.csv（PP 的汇总函数会写这俩）

PP

默认清理（优先命中）

__pycache__/、*.pyc、*.tmp、*.bak、*~、.DS_Store、Thumbs.db

_trash_*、_archive 里超过 N 天的旧回收（防止回收站长成宇宙）

以及：out_dir 内“早于 N 天”的非保护文件（N 可给下拉：3/7/30）

动作默认用你现在最稳的：移入 _trash_YYYYMMDD_HHMMSS，可反悔。

delet

B. “高级清理器”（副按钮/工具菜单）

就是把你 CleanerApp(tk.Tk) 改成 CleanerDialog(tk.Toplevel)（因为 PP 主程序已经有 Tk root 了，不能再起第二个 Tk）。
并且把默认参数“本土化预填”：

目标文件夹 = 当前 PP 的 out_dir

勾上“跳过 _trash_*”（你本来就推荐）

delet

排除规则预填：tables, gate_sites, *_tables, *_gate_sites, *.csv（csv 也一般别乱动）

给 codex 的落地改造清单（直接照着做）
1）把 delet.py 拆成“内核 + UI”

把这些函数提到一个可复用模块里（可以新建 pp_cleaner_core.py，或直接塞进 PP.py 里也行）：

worker_scan(...)（你已有，带 skip_trash）

delet

worker_execute(...)（你已有，支持 trash/delete，保留相对路径移动）

delet

make_trash_dir(...)

delet

再加一个 PP 专用的过滤钩子：should_protect(path_rel_parts) -> bool

伪代码：

PROTECT_TOP = {"tables", "gate_sites"}
PROTECT_SUFFIX = ("_tables", "_gate_sites")

def is_protected(rel_parts: tuple[str, ...]) -> bool:
    if not rel_parts:
        return False
    if rel_parts[0] in PROTECT_TOP:
        return True
    # 任何层级目录名命中 *_tables / *_gate_sites 都保护
    for p in rel_parts:
        if p.endswith(PROTECT_SUFFIX):
            return True
    return False


然后在扫描时：

rel = p.relative_to(root)
if is_protected(rel.parts):
    continue

2）在 PP.py 加一个“快捷清理”函数
def quick_clean_pp(out_dir: str, days: int = 7, action: str = "trash"):
    # 扫描 out_dir 下所有文件
    # 命中条件：mtime < now - days 且不在 protected
    # 执行：move to _trash_...（默认）/ delete（危险，隐藏在高级里）

3）在 graphic（PP）.py 接一个按钮

放在输出目录那一行附近最顺（用户天然会把“输出管理”当成同一块）。
按钮文案建议：🧹 清理旧输出
旁边一个下拉：保留最近：3天 / 7天（默认）/ 30天

点击后：

起线程跑 quick_clean_pp

进度和结果写到 log（你们现在已经有日志系统/窗口了）
