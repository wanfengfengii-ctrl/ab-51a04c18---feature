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


# --------------------------------------------------------------- strand_mode


def test_omitted_strand_mode_response_has_no_strand_fields():
    ref = "ACGTACGATC"
    body = {"reference": ref, "read": ref * 3, "copies": 3, "max_edits": 0}
    omitted = client.post("/api/concatemers/decode", json=body).json()
    explicit = client.post(
        "/api/concatemers/decode", json={**body, "strand_mode": "forward"}
    ).json()
    # Omitted request keeps the historical request echo (no strand_mode).
    assert "strand_mode" not in omitted["request"]
    assert "strand" not in omitted["witness"]
    # Explicit forward echoes the mode and labels the strand.
    assert explicit["request"]["strand_mode"] == "forward"
    assert explicit["witness"]["strand"] == "forward"
    # The decoding itself is identical.
    assert omitted["objective"] == explicit["objective"]
    assert omitted["witness"]["boundaries"] == explicit["witness"]["boundaries"]


def test_reverse_strand_mode():
    ref = "ACGTACGATC"
    rref = "GATCGTACGT"  # reverse complement of ref
    body = {
        "reference": ref,
        "read": rref * 3,
        "copies": 3,
        "max_edits": 0,
        "strand_mode": "reverse",
    }
    resp = client.post("/api/concatemers/decode", json=body)
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["status"] == "unique"
    w = data["witness"]
    assert w["strand"] == "reverse"
    assert w["shift"] == 0
    # The oriented reference is the reverse complement, replayable directly.
    assert data["rotated_reference"] == rref
    seg = w["segments"][0]
    assert seg["reference"] == rref
    assert seg["aligned_reference"] == seg["aligned_read"] == rref
    assert data["request"]["strand_mode"] == "reverse"


def test_auto_resolves_to_reverse():
    ref = "ACGTACGATC"
    rref = "GATCGTACGT"
    body = {
        "reference": ref,
        "read": rref * 3,
        "copies": 3,
        "max_edits": 0,
        "strand_mode": "auto",
    }
    resp = client.post("/api/concatemers/decode", json=body)
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["status"] == "unique"
    assert data["witness"]["strand"] == "reverse"
    assert data["rotated_reference"] == rref


def test_auto_cross_strand_tie_is_ambiguous():
    pal = "AACGTACGTT"  # 10-nt non-periodic reverse-complement palindrome
    body = {
        "reference": pal,
        "read": pal * 3,
        "copies": 3,
        "max_edits": 0,
        "strand_mode": "auto",
    }
    resp = client.post("/api/concatemers/decode", json=body)
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["status"] == "ambiguous"
    strands = [w["strand"] for w in data["witnesses"]]
    assert strands == ["forward", "reverse"]
    assert data["more_witnesses"] is False
    for w in data["witnesses"]:
        assert w["rotated_reference"] == pal


def test_auto_both_infeasible_returns_locatable_422():
    ref = "ACGTACGATC"
    body = {
        "reference": ref,
        "read": "TTGTACGATC" + ref * 2,
        "copies": 3,
        "max_edits": 1,
        "strand_mode": "auto",
    }
    resp = client.post("/api/concatemers/decode", json=body)
    assert resp.status_code == 422
    data = resp.json()
    assert data["status"] == "infeasible"
    assert data["error"] == "constraint_failed"
    assert data["constraint"]["strand_mode"] == "auto"
    assert data["feasible_strands"] == {"forward": [], "reverse": []}
    assert data["nearest"]["strand"] in ("forward", "reverse")


def test_illegal_strand_mode_rejected_as_field_error():
    body = {
        "reference": "ACGTACGATC",
        "read": "ACGTACGATC" * 3,
        "copies": 3,
        "max_edits": 0,
        "strand_mode": "backwards",
    }
    resp = client.post("/api/concatemers/decode", json=body)
    assert resp.status_code == 422
    data = resp.json()
    # Pydantic field-validation envelope, not a constraint failure.
    assert data.get("status") != "infeasible"
    locs = [tuple(err["loc"]) for err in data["detail"]]
    assert ("body", "strand_mode") in locs


def test_omitted_mode_preserves_legacy_failure_body():
    base = {
        "reference": "ACGTACGATC",
        "read": "TTGTACGATC" + "ACGTACGATC" * 2,
        "copies": 3,
        "max_edits": 1,
    }
    # Omitted strand_mode keeps the historical 422 envelope exactly.
    resp = client.post("/api/concatemers/decode", json=base)
    assert resp.status_code == 422
    data = resp.json()
    assert data["status"] == "infeasible"
    assert "feasible_shifts" in data
    assert "feasible_strands" not in data
    assert "strand_mode" not in data["constraint"]
    assert "strand" not in data["nearest"]
    assert "strand_mode" not in data["request"]


def test_explicit_forward_failure_is_strand_aware():
    body = {
        "reference": "ACGTACGATC",
        "read": "TTGTACGATC" + "ACGTACGATC" * 2,
        "copies": 3,
        "max_edits": 1,
        "strand_mode": "forward",
    }
    resp = client.post("/api/concatemers/decode", json=body)
    assert resp.status_code == 422
    data = resp.json()
    assert data["status"] == "infeasible"
    assert "feasible_strands" in data
    assert "feasible_shifts" not in data
    assert data["constraint"]["strand_mode"] == "forward"
    assert data["nearest"]["strand"] == "forward"
