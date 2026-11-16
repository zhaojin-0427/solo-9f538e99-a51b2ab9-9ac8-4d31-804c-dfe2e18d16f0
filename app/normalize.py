"""变异解析、校验、左对齐标准化、等价归并与单倍型合成。

纯函数模块：不依赖 Flask / SQLite，便于单独测试。
所有非法输入抛出 ValidationError，消息面向用户展示；调用方保证整批失败。
变异组合冲突抛出 HaplotypeConflictError（ValidationError 子类），整次拒绝。
"""
from __future__ import annotations

import csv
import io

BASES = set("ACGTN")
MAX_REF_LEN = 4096
MAX_ROWS = 20

SNV = "SNV"
INS = "INS"
DEL = "DEL"

KIND_LABELS = {SNV: "单碱基替换", INS: "插入", DEL: "缺失"}


class ValidationError(Exception):
    """输入非法（文件级或行级）；消息面向用户展示。"""


# ---------------------------------------------------------------- 输入解析

def parse_fasta(text: str):
    """解析单条染色体 FASTA，返回 (chrom, seq)。序列统一大写。"""
    lines = [ln.strip() for ln in text.replace("\r\n", "\n").replace("\r", "\n").split("\n")]
    lines = [ln for ln in lines if ln]
    if not lines:
        raise ValidationError("FASTA 文件为空")
    if not lines[0].startswith(">"):
        raise ValidationError("FASTA 首行必须是以 '>' 开头的标题行")
    if any(ln.startswith(">") for ln in lines[1:]):
        raise ValidationError("只接受单条染色体 FASTA（发现多个标题行）")
    tokens = lines[0][1:].strip().split()
    if not tokens:
        raise ValidationError("FASTA 标题行缺少染色体名称")
    chrom = tokens[0]
    seq = "".join(lines[1:]).replace(" ", "").upper()
    if not seq:
        raise ValidationError("FASTA 序列为空")
    if len(seq) > MAX_REF_LEN:
        raise ValidationError(f"参考序列长度 {len(seq)} 超过上限 {MAX_REF_LEN}")
    bad = sorted(set(seq) - BASES)
    if bad:
        raise ValidationError(f"参考序列含非法字符 {''.join(bad)}（仅允许 A/C/G/T/N）")
    return chrom, seq


REQUIRED_COLS = ("id", "pos", "ref", "alt")


def parse_csv(text: str):
    """解析变异 CSV，返回行字典列表（id/pos/ref/alt/line）。"""
    try:
        rows = [r for r in csv.reader(io.StringIO(text))]
    except csv.Error as exc:
        raise ValidationError(f"CSV 解析失败: {exc}") from exc
    rows = [r for r in rows if any(cell.strip() for cell in r)]
    if not rows:
        raise ValidationError("CSV 文件为空")
    header = [c.strip().lower() for c in rows[0]]
    missing = [c for c in REQUIRED_COLS if c not in header]
    if missing:
        raise ValidationError("CSV 缺少必需列: " + ", ".join(c.upper() for c in missing))
    idx = {c: header.index(c) for c in REQUIRED_COLS}
    data = rows[1:]
    if not data:
        raise ValidationError("CSV 没有数据行")
    if len(data) > MAX_ROWS:
        raise ValidationError(f"CSV 数据行数 {len(data)} 超过上限 {MAX_ROWS}")

    out, seen = [], set()
    for n, r in enumerate(data, start=2):  # n = 含表头的文件行号
        def col(name, r=r):
            i = idx[name]
            return r[i].strip() if i < len(r) else ""

        vid = col("id")
        if not vid:
            raise ValidationError(f"第 {n} 行: ID 为空")
        if any(ch.isspace() for ch in vid) or "," in vid or ";" in vid:
            raise ValidationError(f"第 {n} 行: ID '{vid}' 含空白/逗号/分号")
        if vid in seen:
            raise ValidationError(f"第 {n} 行: ID '{vid}' 重复")
        seen.add(vid)

        pos_s = col("pos")
        if not pos_s.isdigit() or int(pos_s) < 1:
            raise ValidationError(f"第 {n} 行: POS '{pos_s}' 不是从 1 开始的正整数")
        pos = int(pos_s)

        ref, alt = col("ref").upper(), col("alt").upper()
        if not ref or not alt:
            raise ValidationError(f"第 {n} 行: REF/ALT 不能为空")
        for label, v in (("REF", ref), ("ALT", alt)):
            bad = sorted(set(v) - BASES)
            if bad:
                raise ValidationError(f"第 {n} 行: {label} 含非法字符 {''.join(bad)}")
        out.append({"id": vid, "pos": pos, "ref": ref, "alt": alt, "line": n})
    return out


# ---------------------------------------------------------------- 行级校验

def classify(row) -> str:
    """判定变异类型；不支持的写法抛 ValidationError。"""
    ref, alt, n = row["ref"], row["alt"], row["line"]
    if ref == alt:
        raise ValidationError(f"第 {n} 行: REF 与 ALT 相同")
    if len(ref) == 1 and len(alt) == 1:
        return SNV
    if len(ref) == 1 and len(alt) > 1 and alt.startswith(ref):
        return INS
    if len(alt) == 1 and len(ref) > 1 and ref.startswith(alt):
        return DEL
    raise ValidationError(
        f"第 {n} 行: 仅支持单碱基替换，或带锚定碱基的纯插入/缺失"
        f"（REF={ref}, ALT={alt}）"
    )


def check_ref(seq: str, row) -> None:
    """逐行核对 REF 与参考序列是否一致。"""
    pos, ref, n = row["pos"], row["ref"], row["line"]
    if pos > len(seq):
        raise ValidationError(f"第 {n} 行: POS {pos} 超出参考序列长度 {len(seq)}")
    end = pos - 1 + len(ref)
    if end > len(seq):
        raise ValidationError(
            f"第 {n} 行: REF 超出参考序列末端（POS {pos} 起 {len(ref)} 碱基，参考仅 {len(seq)}）"
        )
    actual = seq[pos - 1:end]
    if actual != ref:
        raise ValidationError(
            f"第 {n} 行: REF '{ref}' 与参考序列不符（该处参考为 '{actual}'）"
        )


# ---------------------------------------------------------------- 标准化

def normalize_variant(seq: str, row):
    """左对齐标准化，返回 (pos, ref, alt)（1 基 POS，含锚定碱基的 VCF 写法）。"""
    kind, pos, ref, alt = row["kind"], row["pos"], row["ref"], row["alt"]
    if kind == SNV:
        return pos, ref, alt
    if kind == INS:
        ins = alt[1:]
        p = pos - 1  # 锚定碱基的 0 基下标
        # 插入序列末碱基等于锚定碱基时，整体左移一位并旋转插入序列
        # （等价推导：seq[p] + S == S' + seq[p]  ⟺  S' = seq[p] + S[:-1]）
        while p >= 1 and seq[p] == ins[-1]:
            ins = seq[p] + ins[:-1]
            p -= 1
        return p + 1, seq[p], seq[p] + ins
    # DEL：s 为被删片段的 0 基起始下标（锚定在其左侧）
    dele = ref[1:]
    s = pos
    # 左移后锚定下标为 s-1，必须 >= 1（即 POS>=1 处仍有锚定碱基）
    while s >= 2 and seq[s - 1] == dele[-1]:
        dele = seq[s - 1] + dele[:-1]
        s -= 1
    return s, seq[s - 1] + dele, seq[s - 1]


def apply_variant(seq: str, row) -> str:
    """单条变异作用后的完整序列（等价判定依据）。"""
    i = row["pos"] - 1
    return seq[:i] + row["alt"] + seq[i + len(row["ref"]):]


# ---------------------------------------------------------------- 报告计算

def render_vcf(chrom: str, ref_len: int, groups) -> str:
    lines = [
        "##fileformat=VCFv4.2",
        "##source=variant-normalizer",
        f"##contig=<ID={chrom},length={ref_len}>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    for g in groups:
        lines.append("\t".join([
            chrom, str(g["pos"]), ";".join(g["ids"]), g["ref"], g["alt"], ".", ".", ".",
        ]))
    return "\n".join(lines) + "\n"


def compute_report(fasta_text: str, csv_text: str) -> dict:
    """整批计算：任一非法行 → ValidationError（整批失败，无部分结果）。"""
    chrom, seq = parse_fasta(fasta_text)
    rows = parse_csv(csv_text)
    for row in rows:  # 先全量校验，再计算，保证不产生半份结果
        row["kind"] = classify(row)
        check_ref(seq, row)

    groups = {}  # 完整突变序列 -> 归并组
    for row in rows:
        npos, nref, nalt = normalize_variant(seq, row)
        row["norm"] = (npos, nref, nalt)
        key = apply_variant(seq, row)  # 等价 = 单条变异作用后的完整序列相同
        g = groups.get(key)
        if g is None:
            g = {"pos": npos, "ref": nref, "alt": nalt, "kind": row["kind"], "ids": []}
            groups[key] = g
        g["ids"].append(row["id"])
        row["_key"] = key

    ordered = sorted(groups.values(), key=lambda g: (g["pos"], g["ref"], g["alt"]))
    for g in ordered:
        g["ids"].sort()
        g["kind_label"] = KIND_LABELS[g["kind"]]
    group_index = {id(g): i for i, g in enumerate(ordered)}

    out_rows = [{
        "id": row["id"], "line": row["line"],
        "pos": row["pos"], "ref": row["ref"], "alt": row["alt"],
        "kind": row["kind"], "kind_label": KIND_LABELS[row["kind"]],
        "norm_pos": row["norm"][0], "norm_ref": row["norm"][1], "norm_alt": row["norm"][2],
        "group": group_index[id(groups[row["_key"]])],
    } for row in rows]

    out_groups = [{k: g[k] for k in ("pos", "ref", "alt", "kind", "kind_label", "ids")}
                  for g in ordered]
    return {
        "chrom": chrom,
        "ref_len": len(seq),
        "groups": out_groups,
        "rows": out_rows,
        "vcf": render_vcf(chrom, len(seq), ordered),
    }


# ---------------------------------------------------------------- 单倍型合成

class HaplotypeConflictError(ValidationError):
    """所选变异组合在参考坐标上冲突；携带全部涉及源 ID（去重排序）。"""

    def __init__(self, conflicts):
        self.conflicts = list(conflicts)
        self.ids = sorted({i for c in self.conflicts for i in c["ids"]})
        detail = "；".join(c["message"] for c in self.conflicts)
        super().__init__(
            f"变异组合冲突，整次拒绝（涉及源 ID：{', '.join(self.ids)}）：{detail}"
        )


def _group_edit(g):
    """标准记录 → 参考坐标编辑 (start, end, replacement)。

    0 基半开区间 [start, end)：SNV 替换一个碱基；DEL 删除锚点右侧片段；
    INS 为零宽编辑（start == end == 锚定碱基的 1 基 POS，即插到锚点右侧）。
    """
    kind, pos, ref, alt = g["kind"], g["pos"], g["ref"], g["alt"]
    if kind == SNV:
        return (pos - 1, pos, alt)
    if kind == INS:
        return (pos, pos, alt[1:])
    return (pos, pos + len(ref) - 1, "")  # DEL


def _describe(g):
    ids = "、".join(g["ids"])
    kind, pos, ref, alt = g["kind"], g["pos"], g["ref"], g["alt"]
    if kind == SNV:
        return f"单碱基替换（源 ID：{ids}；POS {pos}：{ref}→{alt}）"
    if kind == INS:
        return f"插入（源 ID：{ids}；POS {pos} 锚定碱基后插入 {alt[1:]}）"
    return (f"缺失（源 ID：{ids}；POS {pos}：{ref}→{alt}，"
            f"删除第 {pos + 1}–{pos + len(ref) - 1} 碱基）")


def _pair_conflict(g1, g2):
    """两条标准记录是否冲突；返回冲突描述 dict 或 None。

    规则（全部按原参考坐标判定）：
    - 替换/缺失的碱基区间不得重叠（被删除的碱基不能再被改动）；
    - 插入锚点碱基不得落在另一缺失删除的碱基中；
    - 不同插入不得占同一断点；
    - 共享 VCF 锚点但实际编辑不冲突（插入/替换与缺失同锚点等）可共存，
      首尾相接的缺失同样可共存（删除区间并集确定）。
    """
    e1, e2 = _group_edit(g1), _group_edit(g2)

    def clash(reason):
        return {"message": f"{_describe(g1)} 与 {_describe(g2)}：{reason}",
                "ids": sorted(set(g1["ids"]) | set(g2["ids"]))}

    # 替换/缺失的碱基足迹（插入为零宽编辑，无足迹）
    f1 = e1 if e1[0] != e1[1] else None
    f2 = e2 if e2[0] != e2[1] else None
    if f1 and f2 and f1[0] < f2[1] and f2[0] < f1[1]:
        return clash("被删除/替换的碱基被再次改动")
    # 插入锚点落在被删碱基上：0 基插入位置 q（锚定碱基的 1 基 POS 即 q），
    # 缺失区间 [s,e) 删除 1 基碱基 s+1..e，锚点被删 ⟺ s < q ≤ e
    for gi, ei, gj, ej in ((g1, e1, g2, e2), (g2, e2, g1, e1)):
        if gi["kind"] == INS and gj["kind"] == DEL and ej[0] < ei[0] <= ej[1]:
            return clash("插入锚点碱基被另一缺失删除")
    if g1["kind"] == g2["kind"] == INS and e1[0] == e2[0]:
        return clash("不同插入占同一断点")
    return None


def synthesize_haplotype(seq, payload, selected_ids):
    """把所选源 ID 对应的标准记录按原参考坐标组合到参考序列上。

    - 映射到同一标准位点（归并组）的源 ID 只应用一次；
    - 全部编辑按参考坐标同时组合（替换与插缺互不位移）；
    - 任一冲突 → HaplotypeConflictError（整次拒绝，消息列出涉及源 ID）。
    返回单倍型完整序列。
    """
    ids = list(dict.fromkeys(selected_ids))  # 去重保序
    if not ids:
        raise ValidationError("未选择任何源 ID")
    rows_by_id = {r["id"]: r for r in payload["rows"]}
    unknown = [i for i in ids if i not in rows_by_id]
    if unknown:
        raise ValidationError("报告中不存在源 ID: " + ", ".join(unknown))

    groups = payload["groups"]
    chosen = {}  # 归并组序号 -> 标准记录（同组只应用一次）
    for i in ids:
        gi = rows_by_id[i]["group"]
        chosen.setdefault(gi, groups[gi])
    selected = [chosen[k] for k in sorted(chosen)]

    conflicts = []
    for a in range(len(selected)):
        for b in range(a + 1, len(selected)):
            c = _pair_conflict(selected[a], selected[b])
            if c:
                conflicts.append(c)
    if conflicts:
        raise HaplotypeConflictError(conflicts)

    # 同一断点处：零宽插入排在同起点编辑之前（先插后删/替），结果唯一确定
    edits = sorted((_group_edit(g) for g in selected), key=lambda e: (e[0], e[1]))
    out, cur = [], 0
    for s, e, rep in edits:
        out.append(seq[cur:s])
        out.append(rep)
        cur = e
    out.append(seq[cur:])
    return "".join(out)
