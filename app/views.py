"""HTTP 路由：批次创建、输入更正、报告查看、单倍型合成。"""
from __future__ import annotations

import json

from flask import (Blueprint, Response, abort, flash, redirect,
                   render_template, request, url_for)

from .db import get_db
from .normalize import (HaplotypeConflictError, ValidationError, parse_fasta,
                        synthesize_haplotype)

bp = Blueprint("main", __name__)

KIND_LABEL = {"fasta": "FASTA", "csv": "CSV"}


# ---------------------------------------------------------------- 辅助

def _batch_or_404(db, bid):
    batch = db.execute("SELECT * FROM batches WHERE id=?", (bid,)).fetchone()
    if batch is None:
        abort(404)
    return batch


def _current_versions(db, bid):
    rows = db.execute(
        "SELECT kind, MAX(version) AS v FROM inputs WHERE batch_id=? GROUP BY kind",
        (bid,),
    ).fetchall()
    d = {r["kind"]: r["v"] for r in rows}
    return d.get("fasta"), d.get("csv")


def _enqueue(db, bid):
    """按当前最新版本组合排队一个计算任务（相同版本组合的未完成任务不重复排）。"""
    fv, cv = _current_versions(db, bid)
    dup = db.execute(
        "SELECT 1 FROM jobs WHERE batch_id=? AND fasta_version=? AND csv_version=?"
        " AND status IN ('pending','running')",
        (bid, fv, cv),
    ).fetchone()
    if dup is None:
        db.execute(
            "INSERT INTO jobs(batch_id, fasta_version, csv_version, status, created_at)"
            " VALUES (?,?,?,'pending', datetime('now'))",
            (bid, fv, cv),
        )


def _read_upload(field):
    """读取上传文件并解码为文本；失败时 flash 并返回 None。"""
    f = request.files.get(field)
    if f is None or not f.filename:
        flash("请选择要上传的文件")
        return None
    try:
        return f.read().decode("utf-8-sig")
    except UnicodeDecodeError:
        flash("文件必须是 UTF-8 编码的文本")
        return None


def _decode_haps(rows):
    """haplotypes 行 → dict，并解析 selected_ids JSON。"""
    return [dict(r, ids=json.loads(r["selected_ids"])) for r in rows]


def _render_fasta(hap, selected_ids, width=60):
    """单倍型记录 → FASTA 文本（头部注明依据报告与所选源 ID）。"""
    header = (f">{hap['chrom']}_hap{hap['id']} "
              f"report={hap['report_id']} ids={','.join(selected_ids)}")
    seq = hap["sequence"]
    lines = [header] + [seq[i:i + width] for i in range(0, len(seq), width)]
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------- 页面

@bp.get("/")
def index():
    db = get_db()
    batches = db.execute("SELECT * FROM batches ORDER BY id DESC").fetchall()
    return render_template("index.html", batches=batches)


@bp.post("/batches")
def create_batch():
    name = request.form.get("name", "").strip() or "未命名批次"
    fasta_text = _read_upload("fasta")
    csv_text = _read_upload("csv")
    if fasta_text is None or csv_text is None:
        return redirect(url_for("main.index"))

    db = get_db()
    cur = db.execute(
        "INSERT INTO batches(name, created_at) VALUES (?, datetime('now'))", (name,)
    )
    bid = cur.lastrowid
    for kind, text in (("fasta", fasta_text), ("csv", csv_text)):
        db.execute(
            "INSERT INTO inputs(batch_id, kind, version, content, created_at)"
            " VALUES (?,?,1,?, datetime('now'))",
            (bid, kind, text),
        )
    _enqueue(db, bid)
    db.commit()
    flash(f"批次 #{bid} 已创建（FASTA v1 / CSV v1），计算已排队")
    return redirect(url_for("main.batch_detail", bid=bid))


@bp.post("/batches/<int:bid>/inputs/<kind>")
def correct_input(bid, kind):
    """整份更正某一份输入；计算进行中也可提交。"""
    if kind not in KIND_LABEL:
        abort(404)
    db = get_db()
    _batch_or_404(db, bid)
    text = _read_upload("file")
    if text is None:
        return redirect(url_for("main.batch_detail", bid=bid))

    cur = db.execute(
        "SELECT MAX(version) AS v FROM inputs WHERE batch_id=? AND kind=?",
        (bid, kind),
    ).fetchone()["v"]
    db.execute(
        "INSERT INTO inputs(batch_id, kind, version, content, created_at)"
        " VALUES (?,?,?,?, datetime('now'))",
        (bid, kind, cur + 1, text),
    )
    _enqueue(db, bid)
    db.commit()
    flash(f"{KIND_LABEL[kind]} 第 {cur + 1} 版已提交，已按最新版本组合排队重新计算")
    return redirect(url_for("main.batch_detail", bid=bid))


@bp.get("/batches/<int:bid>")
def batch_detail(bid):
    db = get_db()
    batch = _batch_or_404(db, bid)
    fv, cv = _current_versions(db, bid)
    report = db.execute(
        "SELECT * FROM reports WHERE batch_id=? AND fasta_version=? AND csv_version=?"
        " ORDER BY id DESC LIMIT 1",
        (bid, fv, cv),
    ).fetchone()
    computing = db.execute(
        "SELECT COUNT(*) AS n FROM jobs WHERE batch_id=? AND fasta_version=?"
        " AND csv_version=? AND status IN ('pending','running')",
        (bid, fv, cv),
    ).fetchone()["n"] > 0
    reports = db.execute(
        "SELECT * FROM reports WHERE batch_id=? ORDER BY id DESC LIMIT 50", (bid,)
    ).fetchall()
    jobs = db.execute(
        "SELECT * FROM jobs WHERE batch_id=? ORDER BY id DESC LIMIT 20", (bid,)
    ).fetchall()
    haplotypes = _decode_haps(db.execute(
        "SELECT * FROM haplotypes WHERE batch_id=? ORDER BY id DESC LIMIT 50",
        (bid,),
    ).fetchall())
    payload = (json.loads(report["payload"])
               if report and report["status"] == "ok" else None)
    return render_template(
        "batch.html", batch=batch, fv=fv, cv=cv, report=report, payload=payload,
        computing=computing, reports=reports, jobs=jobs, haplotypes=haplotypes,
    )


@bp.get("/batches/<int:bid>/reports/<int:rid>")
def report_detail(bid, rid):
    db = get_db()
    batch = _batch_or_404(db, bid)
    report = db.execute(
        "SELECT * FROM reports WHERE id=? AND batch_id=?", (rid, bid)
    ).fetchone()
    if report is None:
        abort(404)
    fv, cv = _current_versions(db, bid)
    is_current = (report["fasta_version"] == fv and report["csv_version"] == cv)
    payload = (json.loads(report["payload"])
               if report["status"] == "ok" else None)
    haplotypes = _decode_haps(db.execute(
        "SELECT * FROM haplotypes WHERE report_id=? ORDER BY id DESC", (rid,)
    ).fetchall())
    return render_template(
        "report.html", batch=batch, report=report, payload=payload,
        is_current=is_current, fv=fv, cv=cv, haplotypes=haplotypes,
    )


@bp.get("/batches/<int:bid>/reports/<int:rid>/vcf")
def report_vcf(bid, rid):
    db = get_db()
    report = db.execute(
        "SELECT * FROM reports WHERE id=? AND batch_id=?", (rid, bid)
    ).fetchone()
    if report is None or report["status"] != "ok":
        abort(404)
    vcf = json.loads(report["payload"])["vcf"]
    return Response(
        vcf,
        content_type="text/plain; charset=utf-8",
        headers={"Content-Disposition":
                 f"attachment; filename=batch{bid}_report{rid}.vcf"},
    )


# ---------------------------------------------------------------- 单倍型合成

@bp.post("/batches/<int:bid>/reports/<int:rid>/haplotypes")
def create_haplotype(bid, rid):
    """勾选源 ID 合成单倍型；冲突整次拒绝（留失败记录，不留部分序列）。

    依据报告固定为其 (fasta_version, csv_version)：提交前若输入已更正，
    结果照常保存可查，但因版本不再匹配而标为"已过期"，不会标为当前。
    """
    db = get_db()
    _batch_or_404(db, bid)
    report = db.execute(
        "SELECT * FROM reports WHERE id=? AND batch_id=?", (rid, bid)
    ).fetchone()
    if report is None or report["status"] != "ok":
        abort(404)

    ids = sorted(set(request.form.getlist("ids")))
    payload = json.loads(report["payload"])
    fasta = db.execute(
        "SELECT content FROM inputs WHERE batch_id=? AND kind='fasta' AND version=?",
        (bid, report["fasta_version"]),
    ).fetchone()
    _, seq = parse_fasta(fasta["content"])

    try:
        hap_seq = synthesize_haplotype(seq, payload, ids)
        status, error, sequence = "ok", None, hap_seq
    except HaplotypeConflictError as exc:
        status, error, sequence = "failed", str(exc), None
    except ValidationError as exc:
        # 空选择 / 未知源 ID：表单错误，直接提示，不产生合成记录
        flash(str(exc))
        return redirect(url_for("main.report_detail", bid=bid, rid=rid))

    cur = db.execute(
        "INSERT INTO haplotypes(batch_id, report_id, fasta_version, csv_version,"
        " status, selected_ids, error, chrom, seq_length, sequence, created_at)"
        " VALUES (?,?,?,?,?,?,?,?,?,?, datetime('now'))",
        (bid, rid, report["fasta_version"], report["csv_version"], status,
         json.dumps(ids, ensure_ascii=False), error, payload["chrom"],
         len(sequence) if sequence is not None else None, sequence),
    )
    db.commit()
    hid = cur.lastrowid
    if status == "ok":
        flash(f"单倍型 #{hid} 已生成：{len(sequence)} bp，"
              f"依据报告 #{rid}，所选 {len(ids)} 个源 ID")
    else:
        flash(f"单倍型 #{hid} 合成失败：{error}")
    return redirect(url_for("main.haplotype_detail", bid=bid, hid=hid))


@bp.get("/batches/<int:bid>/haplotypes/<int:hid>")
def haplotype_detail(bid, hid):
    db = get_db()
    batch = _batch_or_404(db, bid)
    hap = db.execute(
        "SELECT * FROM haplotypes WHERE id=? AND batch_id=?", (hid, bid)
    ).fetchone()
    if hap is None:
        abort(404)
    fv, cv = _current_versions(db, bid)
    is_current = (hap["status"] == "ok" and hap["fasta_version"] == fv
                  and hap["csv_version"] == cv)
    selected_ids = json.loads(hap["selected_ids"])
    fasta_text = (_render_fasta(hap, selected_ids)
                  if hap["status"] == "ok" else None)
    return render_template(
        "haplotype.html", batch=batch, hap=hap, selected_ids=selected_ids,
        is_current=is_current, fv=fv, cv=cv, fasta_text=fasta_text,
    )


@bp.get("/batches/<int:bid>/haplotypes/<int:hid>/fasta")
def haplotype_fasta(bid, hid):
    db = get_db()
    hap = db.execute(
        "SELECT * FROM haplotypes WHERE id=? AND batch_id=?", (hid, bid)
    ).fetchone()
    if hap is None or hap["status"] != "ok":
        abort(404)
    text = _render_fasta(hap, json.loads(hap["selected_ids"]))
    return Response(
        text,
        content_type="text/plain; charset=utf-8",
        headers={"Content-Disposition":
                 f"attachment; filename=batch{bid}_haplotype{hid}.fa"},
    )
