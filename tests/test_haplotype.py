"""单倍型合成测试：参考坐标组合、冲突规则、整次拒绝、版本派生"当前"语义。"""
import io
import json

import pytest

from app import create_app
from app.db import connect
from app.normalize import (HaplotypeConflictError, ValidationError,
                           compute_report, synthesize_haplotype)
from app.worker import run_pending_jobs

REF = "ACGTACGT"          # 1:A 2:C 3:G 4:T 5:A 6:C 7:G 8:T
FASTA = f">chr1\n{REF}\n"


def csv_of(*rows):
    return "ID,POS,REF,ALT\n" + "\n".join(",".join(map(str, r)) for r in rows) + "\n"


def payload_of(csv_text, fasta=FASTA):
    return compute_report(fasta, csv_text)


# ---------------------------------------------------------------- 纯函数：组合

def test_same_locus_records_applied_once():
    # i1/i2 是同聚物中同一插入的两种写法，归并到同一标准位点
    payload = payload_of(csv_of(("i1", 3, "A", "AA"), ("i2", 7, "A", "AA")),
                         fasta=">c\n" + "A" * 10 + "\n")
    assert len(payload["groups"]) == 1
    hap = synthesize_haplotype("A" * 10, payload, ["i1", "i2"])
    assert hap == "A" * 11            # 只应用一次，而非两次


def test_shared_anchor_ins_and_del_coexist():
    # 插入与缺失共享锚点 C2，但实际编辑不冲突：插 T 到 C2 后、删 G3
    payload = payload_of(csv_of(("i1", 2, "C", "CT"), ("d1", 2, "CG", "C")))
    assert synthesize_haplotype(REF, payload, ["i1", "d1"]) == "ACTTACGT"


def test_shared_anchor_snv_and_del_coexist():
    # 替换改动锚点碱基本身，缺失删其右侧碱基：实际编辑不重叠
    payload = payload_of(csv_of(("s1", 2, "C", "T"), ("d1", 2, "CG", "C")))
    assert synthesize_haplotype(REF, payload, ["s1", "d1"]) == "ATTACGT"


def test_snv_and_ins_same_base_coexist():
    payload = payload_of(csv_of(("s1", 2, "C", "T"), ("i1", 2, "C", "CG")))
    assert synthesize_haplotype(REF, payload, ["s1", "i1"]) == "ATGGTACGT"


def test_adjacent_deletions_coexist():
    # 首尾相接的缺失：删除区间不重叠，并集确定
    payload = payload_of(csv_of(("d1", 2, "CG", "C"), ("d2", 3, "GT", "G")))
    assert synthesize_haplotype(REF, payload, ["d1", "d2"]) == "ACACGT"


def test_combination_is_order_independent():
    payload = payload_of(csv_of(
        ("s1", 1, "A", "T"),      # 替换首碱基
        ("i1", 5, "A", "AC"),     # A5 后插 C
        ("d1", 7, "GT", "G"),     # 删 T8
    ))
    expect = "TCGTACCG"
    assert synthesize_haplotype(REF, payload, ["s1", "i1", "d1"]) == expect
    assert synthesize_haplotype(REF, payload, ["d1", "s1", "i1"]) == expect
    assert synthesize_haplotype(REF, payload, ["i1", "d1", "s1"]) == expect


def test_dedup_group_combined_with_other_groups():
    payload = payload_of(csv_of(
        ("i1", 3, "A", "AA"), ("i2", 7, "A", "AA"),   # 同一标准位点
        ("s1", 5, "A", "C"),
    ), fasta=">c\n" + "A" * 10 + "\n")
    hap = synthesize_haplotype("A" * 10, payload, ["i1", "i2", "s1"])
    # 插入左对齐到第 1 位后（插一个 A），第 5 位替换为 C → C 落在第 6 位
    assert hap == "A" * 5 + "C" + "A" * 5


# ---------------------------------------------------------------- 纯函数：冲突

def test_snv_on_deleted_base_conflicts():
    payload = payload_of(csv_of(("d1", 2, "CGT", "C"), ("s1", 3, "G", "A")))
    with pytest.raises(HaplotypeConflictError) as ei:
        synthesize_haplotype(REF, payload, ["d1", "s1"])
    assert ei.value.ids == ["d1", "s1"]
    assert "被删除" in str(ei.value) and "d1" in str(ei.value) and "s1" in str(ei.value)


def test_ins_anchor_on_deleted_base_conflicts():
    payload = payload_of(csv_of(
        ("d1", 2, "CGT", "C"),     # 删除 G3、T4
        ("i1", 3, "G", "GA"),      # 锚点 G3 被删除
        ("i2", 4, "T", "TA"),      # 锚点 T4 被删除
        ("ok", 2, "C", "CA"),      # 与缺失共享锚点：可共存
    ))
    with pytest.raises(HaplotypeConflictError) as ei:
        synthesize_haplotype(REF, payload, ["d1", "i1"])
    assert ei.value.ids == ["d1", "i1"]
    with pytest.raises(HaplotypeConflictError):
        synthesize_haplotype(REF, payload, ["d1", "i2"])
    # 共享锚点的插入不冲突
    assert synthesize_haplotype(REF, payload, ["d1", "ok"]) == "ACAACGT"


def test_two_insertions_same_breakpoint_conflict():
    payload = payload_of(csv_of(("i1", 5, "A", "AG"), ("i2", 5, "A", "AT")))
    assert len(payload["groups"]) == 2          # 插入序列不同 → 不同标准位点
    with pytest.raises(HaplotypeConflictError) as ei:
        synthesize_haplotype(REF, payload, ["i1", "i2"])
    assert ei.value.ids == ["i1", "i2"]
    assert "断点" in str(ei.value)


def test_overlapping_deletions_conflict():
    payload = payload_of(csv_of(("d1", 2, "CGT", "C"), ("d2", 3, "GTA", "G")))
    with pytest.raises(HaplotypeConflictError) as ei:
        synthesize_haplotype(REF, payload, ["d1", "d2"])
    assert ei.value.ids == ["d1", "d2"]


def test_conflict_error_lists_all_involved_ids():
    # d1 与 s1 冲突、d1 与 i1 冲突：涉及源 ID 取并集
    payload = payload_of(csv_of(
        ("d1", 2, "CGT", "C"), ("s1", 3, "G", "A"), ("i1", 4, "T", "TA"),
    ))
    with pytest.raises(HaplotypeConflictError) as ei:
        synthesize_haplotype(REF, payload, ["d1", "s1", "i1"])
    assert ei.value.ids == ["d1", "i1", "s1"]
    assert "涉及源 ID" in str(ei.value)


def test_unselected_conflicting_ids_do_not_block():
    # 只勾选互不冲突的子集即可合成；未勾选的冲突记录不参与
    payload = payload_of(csv_of(("d1", 2, "CGT", "C"), ("s1", 3, "G", "A"),
                                ("s2", 8, "T", "A")))
    assert synthesize_haplotype(REF, payload, ["d1", "s2"]) == "ACACGA"


def test_empty_selection_rejected():
    payload = payload_of(csv_of(("s1", 2, "C", "T")))
    with pytest.raises(ValidationError, match="未选择"):
        synthesize_haplotype(REF, payload, [])


def test_unknown_id_rejected():
    payload = payload_of(csv_of(("s1", 2, "C", "T")))
    with pytest.raises(ValidationError, match="不存在"):
        synthesize_haplotype(REF, payload, ["s1", "ghost"])


# ---------------------------------------------------------------- 应用层

APP_FASTA = ">chr1\nACGTACGTACGT\n"   # 12 bp
APP_CSV = "ID,POS,REF,ALT\ns1,2,C,T\nd1,5,AC,A\ni1,9,A,AG\n"
APP_HAP = "ATGTAGTAGCGT"              # s1 替换、d1 删 C6、i1 在 A9 后插 G


@pytest.fixture()
def app(tmp_path):
    return create_app({
        "TESTING": True,
        "DATABASE": str(tmp_path / "t.db"),
        "WORKER_ENABLED": False,
        "COMPUTE_DELAY_SECONDS": 0,
    })


@pytest.fixture()
def client(app):
    return app.test_client()


def _upload(client, fasta=APP_FASTA, csv_text=APP_CSV):
    return client.post("/batches", data={
        "name": "t",
        "fasta": (io.BytesIO(fasta.encode()), "ref.fa"),
        "csv": (io.BytesIO(csv_text.encode()), "vars.csv"),
    })


def _correct_csv(client, text, bid=1):
    return client.post(f"/batches/{bid}/inputs/csv",
                       data={"file": (io.BytesIO(text.encode()), "new.csv")})


def _synthesize(client, ids, bid=1, rid=1):
    return client.post(f"/batches/{bid}/reports/{rid}/haplotypes",
                       data={"ids": ids}, follow_redirects=False)


def _haps(app, bid=1):
    conn = connect(app.config["DATABASE"])
    rows = conn.execute(
        "SELECT * FROM haplotypes WHERE batch_id=? ORDER BY id", (bid,)).fetchall()
    conn.close()
    return rows


def _run(app):
    run_pending_jobs(app.config["DATABASE"])


def test_haplotype_happy_path_and_fasta_download(app, client):
    _upload(client)
    _run(app)

    page = client.get("/batches/1/reports/1").data.decode()
    assert 'name="ids"' in page                 # 当前成功报告页可勾选源 ID

    resp = _synthesize(client, ["s1", "d1", "i1"])
    assert resp.status_code == 302
    assert resp.headers["Location"].endswith("/batches/1/haplotypes/1")

    haps = _haps(app)
    assert len(haps) == 1 and haps[0]["status"] == "ok"
    assert haps[0]["report_id"] == 1
    assert json.loads(haps[0]["selected_ids"]) == ["d1", "i1", "s1"]
    assert haps[0]["sequence"] == APP_HAP and haps[0]["seq_length"] == 12

    page = client.get("/batches/1/haplotypes/1").data.decode()
    assert APP_HAP in page and "当前" in page
    assert "d1, i1, s1" in page and "报告" in page

    resp = client.get("/batches/1/haplotypes/1/fasta")
    assert resp.status_code == 200
    body = resp.data.decode()
    assert body.startswith(">chr1_hap1 report=1 ids=d1,i1,s1\n")
    assert APP_HAP in body.replace("\n", "")

    batch_page = client.get("/batches/1").data.decode()
    assert "单倍型" in batch_page and "当前" in batch_page


def test_same_locus_applied_once_via_web(app, client):
    _upload(client, fasta=">c\n" + "A" * 10 + "\n",
            csv_text="ID,POS,REF,ALT\ni1,3,A,AA\ni2,7,A,AA\n")
    _run(app)
    _synthesize(client, ["i1", "i2"])
    haps = _haps(app)
    assert haps[0]["sequence"] == "A" * 11      # 同一标准位点只应用一次


def test_conflict_rejected_lists_ids_and_keeps_no_partial_sequence(app, client):
    _upload(client, csv_text="ID,POS,REF,ALT\nd1,2,CGT,C\ns1,3,G,A\n")
    _run(app)
    resp = _synthesize(client, ["d1", "s1"])
    assert resp.headers["Location"].endswith("/batches/1/haplotypes/1")

    haps = _haps(app)
    assert len(haps) == 1 and haps[0]["status"] == "failed"
    assert haps[0]["sequence"] is None          # 失败不留部分序列
    assert haps[0]["seq_length"] is None
    assert "d1" in haps[0]["error"] and "s1" in haps[0]["error"]

    page = client.get("/batches/1/haplotypes/1").data.decode()
    assert "失败" in page and "d1" in page and "s1" in page
    assert "整次拒绝" in page
    assert client.get("/batches/1/haplotypes/1/fasta").status_code == 404
    batch_page = client.get("/batches/1").data.decode()
    assert "失败" in batch_page


def test_stale_after_correction_but_still_viewable(app, client):
    _upload(client)
    _run(app)
    _synthesize(client, ["s1", "d1", "i1"])
    assert "当前" in client.get("/batches/1/haplotypes/1").data.decode()

    _correct_csv(client, "ID,POS,REF,ALT\ns9,4,T,A\n")
    _run(app)                                    # 报告 #2 成为当前结论

    page = client.get("/batches/1/haplotypes/1").data.decode()
    assert "已过期" in page and "当前</span>" not in page   # 旧依据结果不得标为当前
    assert APP_HAP in page                                # 但仍可查
    batch_page = client.get("/batches/1").data.decode()
    assert "已过期" in batch_page

    # 基于新报告合成的新单倍型才是当前
    _synthesize(client, ["s9"], rid=2)
    page = client.get("/batches/1/haplotypes/2").data.decode()
    assert "当前" in page


def test_interleaved_correction_old_basis_saved_but_not_current(app, client):
    """更正与生成交错：按旧报告生成的结果照常保存可查，但立即标为已过期。"""
    _upload(client)
    _run(app)
    _correct_csv(client, "ID,POS,REF,ALT\ns9,4,T,A\n")   # 更正先到
    _synthesize(client, ["s1"], rid=1)                   # 旧报告上的生成后至

    haps = _haps(app)
    assert len(haps) == 1 and haps[0]["status"] == "ok"
    assert (haps[0]["fasta_version"], haps[0]["csv_version"]) == (1, 1)
    page = client.get("/batches/1/haplotypes/1").data.decode()
    assert "已过期" in page and "当前</span>" not in page

    _run(app)                                            # 新版本报告完成
    page = client.get("/batches/1/haplotypes/1").data.decode()
    assert "已过期" in page                              # 依然不是当前


def test_empty_and_unknown_selection_create_no_record(app, client):
    _upload(client)
    _run(app)
    resp = _synthesize(client, [])
    assert resp.headers["Location"].endswith("/batches/1/reports/1")
    resp = _synthesize(client, ["ghost"])
    assert resp.headers["Location"].endswith("/batches/1/reports/1")
    assert _haps(app) == []                              # 表单错误不产生记录


def test_synthesis_form_only_on_current_report(app, client):
    _upload(client)
    _run(app)
    _correct_csv(client, "ID,POS,REF,ALT\ns9,4,T,A\n")
    _run(app)
    old_page = client.get("/batches/1/reports/1").data.decode()
    assert 'name="ids"' not in old_page                  # 过期报告不开放勾选
    new_page = client.get("/batches/1/reports/2").data.decode()
    assert 'name="ids"' in new_page


def test_synthesis_rejects_failed_or_missing_report(app, client):
    _upload(client, csv_text="ID,POS,REF,ALT\nbad,5,C,G\n")   # 整批失败
    _run(app)
    assert client.post("/batches/1/reports/1/haplotypes",
                       data={"ids": ["bad"]}).status_code == 404
    assert client.post("/batches/1/reports/99/haplotypes",
                       data={"ids": ["x"]}).status_code == 404
    assert client.get("/batches/1/haplotypes/99").status_code == 404
    assert _haps(app) == []
