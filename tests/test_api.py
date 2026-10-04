"""HTTP-level tests for the decode API."""

from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


def test_health():
    resp = client.get("/health")
    assert resp.status_code == 200
    assert resp.json()["status"] == "ok"


def test_decode_unique():
    ref = "ACGTACGAT"
    body = {
        "reference": ref,
        "read": ref * 4,
        "copies": 4,
        "max_edits": 0,
    }
    resp = client.post("/api/concatemers/decode", json=body)
    assert resp.status_code == 200
    data = resp.json()
    assert data["status"] == "unique"
    assert data["objective"] == {"total_edits": 0, "max_segment_edits": 0}
    assert data["rotated_reference"] == ref
    seg = data["witness"]["segments"][0]
    assert seg["cigar"] == "9M"
    # replay rows are present and consistent
    assert seg["aligned_reference"] == seg["aligned_read"] == ref


def test_decode_with_indel_and_substitution_smoke():
    # substitution + deletion + insertion across three copies (30 nt total)
    ref = "ACGTACGATC"  # 10 nt
    read = (
        "ATGTACGATC"    # C -> T substitution, 10 nt
        + "ACGTACGAT"   # trailing C deleted, 9 nt
        + "ACGTACGATCA"  # A inserted, 11 nt
    )
    assert len(read) == 30
    body = {"reference": ref, "read": read, "copies": 3, "max_edits": 1}
    resp = client.post("/api/concatemers/decode", json=body)
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["objective"] == {"total_edits": 3, "max_segment_edits": 1}
    witness = data.get("witness") or data["witnesses"][0]
    lengths = [e - b for b, e in witness["boundaries"]]
    assert lengths == [10, 9, 11]
    kinds = {op for s in witness["segments"] for op in s["cigar"] if op.isalpha()}
    assert {"M", "D", "I"} <= kinds


def test_decode_ambiguous_returns_two_witnesses():
    body = {
        "reference": "AAAAAAAAAA",
        "read": "A" * 30,
        "copies": 3,
        "max_edits": 0,
    }
    resp = client.post("/api/concatemers/decode", json=body)
    assert resp.status_code == 200
    data = resp.json()
    assert data["status"] == "ambiguous"
    assert len(data["witnesses"]) == 2
    assert data["witnesses"][0]["shift"] < data["witnesses"][1]["shift"]


def test_decode_infeasible_422_with_location():
    body = {
        "reference": "ACGTACGATC",
        # first copy carries two substitutions (A->T, C->T)
        "read": "TTGTACGATC" + "ACGTACGATC" * 2,
        "copies": 3,
        "max_edits": 1,
    }
    resp = client.post("/api/concatemers/decode", json=body)
    assert resp.status_code == 422
    data = resp.json()
    assert data["status"] == "infeasible"
    assert data["error"] == "constraint_failed"
    bad = data["nearest"]["violating_segments"]
    assert bad and bad[0]["required_edits"] == 2


def test_validation_reference_length():
    body = {"reference": "ACGT", "read": "A" * 30, "copies": 3, "max_edits": 0}
    resp = client.post("/api/concatemers/decode", json=body)
    assert resp.status_code == 422


def test_validation_read_length():
    body = {"reference": "ACGTACGT", "read": "A" * 20, "copies": 3, "max_edits": 0}
    resp = client.post("/api/concatemers/decode", json=body)
    assert resp.status_code == 422


def test_validation_copies_range():
    body = {
        "reference": "ACGTACGT",
        "read": "A" * 30,
        "copies": 9,
        "max_edits": 0,
    }
    resp = client.post("/api/concatemers/decode", json=body)
    assert resp.status_code == 422


def test_validation_edits_range():
    body = {
        "reference": "ACGTACGT",
        "read": "A" * 30,
        "copies": 3,
        "max_edits": 4,
    }
    resp = client.post("/api/concatemers/decode", json=body)
    assert resp.status_code == 422


def test_validation_bad_characters():
    body = {
        "reference": "ACGTACGN",
        "read": "A" * 30,
        "copies": 3,
        "max_edits": 0,
    }
    resp = client.post("/api/concatemers/decode", json=body)
    assert resp.status_code == 422


def test_lowercase_is_normalized():
    ref = "acgtacgatc"
    body = {"reference": ref, "read": ref.upper() * 3, "copies": 3, "max_edits": 0}
    resp = client.post("/api/concatemers/decode", json=body)
    assert resp.status_code == 200
    assert resp.json()["status"] == "unique"


# -------------------------------------------------------------- strand modes


def test_omitted_strand_mode_is_legacy_response():
    ref = "ACGTACGATC"
    body = {"reference": ref, "read": ref * 3, "copies": 3, "max_edits": 0}
    resp = client.post("/api/concatemers/decode", json=body)
    assert resp.status_code == 200
    data = resp.json()
    # No strand fields leak into the legacy response.
    assert "strand" not in data["witness"]
    assert "strand_mode" not in data["request"]
    assert data["rotated_reference"] == ref


def test_explicit_forward_matches_omitted_response_minus_strand_fields():
    ref = "ACGTACGATC"
    payload = {"reference": ref, "read": ref * 3, "copies": 3, "max_edits": 0}
    legacy = client.post(
        "/api/concatemers/decode", json=payload
    ).json()
    explicit = client.post(
        "/api/concatemers/decode", json={**payload, "strand_mode": "forward"}
    ).json()
    assert explicit["request"]["strand_mode"] == "forward"
    assert explicit["witness"]["strand"] == "forward"
    # Strip the new fields and compare the envelopes for equality.
    explicit["request"].pop("strand_mode")
    explicit["witness"].pop("strand")
    explicit.pop("ordering")
    legacy.pop("ordering")
    assert explicit == legacy


def test_reverse_mode_api():
    from app.solver import reverse_complement, rotate

    ref = "ACGTACGATC"
    oriented = rotate(reverse_complement(ref), 4)
    body = {
        "reference": ref,
        "read": oriented * 3,
        "copies": 3,
        "max_edits": 0,
        "strand_mode": "reverse",
    }
    resp = client.post("/api/concatemers/decode", json=body)
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["status"] == "unique"
    witness = data["witness"]
    assert witness["strand"] == "reverse"
    assert witness["shift"] == 4
    # The response hands back the oriented, rotated reference so the CIGAR
    # replays directly against it and the submitted read fragment.
    assert data["rotated_reference"] == oriented
    seg = witness["segments"][0]
    assert seg["reference"] == oriented
    assert seg["aligned_reference"] == seg["aligned_read"] == oriented
    assert data["request"]["strand_mode"] == "reverse"


def test_auto_mode_resolves_direction():
    from app.solver import reverse_complement, rotate

    ref = "ACGTACGATC"
    read = rotate(reverse_complement(ref), 4) * 3
    body = {
        "reference": ref,
        "read": read,
        "copies": 3,
        "max_edits": 0,
        "strand_mode": "auto",
    }
    resp = client.post("/api/concatemers/decode", json=body)
    assert resp.status_code == 200
    data = resp.json()
    assert data["status"] == "unique"
    assert data["witness"]["strand"] == "reverse"


def test_auto_cross_direction_ambiguity_smoke():
    # Reference equal to its own reverse complement: both directions tie.
    pal = "ACGATATCGT"
    body = {
        "reference": pal,
        "read": pal * 3,
        "copies": 3,
        "max_edits": 0,
        "strand_mode": "auto",
    }
    resp = client.post("/api/concatemers/decode", json=body)
    assert resp.status_code == 200
    data = resp.json()
    assert data["status"] == "ambiguous"
    w0, w1 = data["witnesses"]
    assert (w0["strand"], w1["strand"]) == ("forward", "reverse")
    assert data["request"]["strand_mode"] == "auto"
    # Stable (strand, shift, boundaries, CIGAR) ordering.
    assert (w0["strand"], w0["shift"], w0["boundaries"]) <= (
        w1["strand"],
        w1["shift"],
        w1["boundaries"],
    )


def test_auto_both_infeasible_is_422_constraint_failure():
    body = {
        "reference": "ACGTACGATC",
        "read": "TTGTACGATC" + "G" * 10 + "C" * 10,
        "copies": 3,
        "max_edits": 1,
        "strand_mode": "auto",
    }
    resp = client.post("/api/concatemers/decode", json=body)
    assert resp.status_code == 422
    data = resp.json()
    assert data["status"] == "infeasible"
    assert data["error"] == "constraint_failed"
    assert data["constraint"]["name"] == "per_segment_edit_budget"
    assert data["nearest"] is None


def test_invalid_strand_mode_rejected_by_field():
    body = {
        "reference": "ACGTACGATC",
        "read": "A" * 30,
        "copies": 3,
        "max_edits": 0,
        "strand_mode": "backwards",
    }
    resp = client.post("/api/concatemers/decode", json=body)
    assert resp.status_code == 422
    # Pydantic locates the offending field.
    detail = resp.json()["detail"]
    locs = {tuple(err["loc"]) for err in detail}
    assert ("body", "strand_mode") in locs
