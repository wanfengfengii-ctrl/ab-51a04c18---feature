#!/usr/bin/env python3
"""One-shot verification run by the ``verify`` compose service.

It performs, in order, exiting non-zero on the first failure:

1. wait for the API ``/health`` endpoint to report healthy,
2. run the code test suite (pytest),
3. sanity-check the running build (importable app + dependency versions),
4. run a decode smoke test against the live API that exercises
   insertions, deletions and substitutions together, plus unique,
   ambiguous and infeasible responses;
5. exercise strand_mode: a reverse-complement read decoded in forward and
   reverse, an auto request whose direction is uniquely determined, and a
   cross-strand tie that must come back ambiguous.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
import urllib.error
import urllib.request

API = os.environ.get("API_BASE", "http://api:8000").rstrip("/")
HEALTH_URL = f"{API}/health"
DECODE_URL = f"{API}/api/concatemers/decode"

# Repository root (this file lives in <root>/scripts/verify.py).  In the
# image the tree is copied under /srv, so the same derivation holds there.
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def step(title: str) -> None:
    print(f"\n=== verify: {title} ===", flush=True)


def wait_for_health(timeout: float = 60.0) -> bool:
    step(f"waiting for API health at {HEALTH_URL}")
    deadline = time.time() + timeout
    last_error = None
    while time.time() < deadline:
        try:
            with urllib.request.urlopen(HEALTH_URL, timeout=2) as resp:
                if resp.status == 200:
                    body = json.loads(resp.read().decode())
                    if body.get("status") == "ok":
                        print("health ok:", body)
                        return True
        except (urllib.error.URLError, ConnectionError, OSError) as exc:
            last_error = exc
        time.sleep(1.0)
    print(f"health check timed out after {timeout:.0f}s: {last_error}")
    return False


def run_tests() -> bool:
    step("running code tests (pytest)")
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", "tests", "-q", "-p", "no:cacheprovider"],
        cwd=ROOT,
    )
    if proc.returncode != 0:
        print("pytest failed")
        return False
    print("pytest passed")
    return True


def check_build() -> bool:
    step("sanity-checking the running build")
    code = (
        "import fastapi, uvicorn, pydantic; "
        "from app.main import app; "
        "print('fastapi', fastapi.__version__, "
        "'uvicorn', uvicorn.__version__, "
        "'pydantic', pydantic.VERSION, "
        "'routes', len(app.routes))"
    )
    proc = subprocess.run([sys.executable, "-c", code], cwd=ROOT)
    if proc.returncode != 0:
        print("build sanity check failed")
        return False
    print("build sanity check passed")
    return True


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGT", "TGCA")
    return seq.translate(table)[::-1]


def post(payload: dict):
    req = urllib.request.Request(
        DECODE_URL,
        data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=10) as resp:
            return resp.status, json.loads(resp.read().decode())
    except urllib.error.HTTPError as exc:
        return exc.code, json.loads(exc.read().decode())


def smoke_decode() -> bool:
    step("decode smoke against live API (insertion + deletion + substitution)")

    # 10 nt reference; three copies with one edit of each kind.
    ref = "ACGTACGATC"
    read = (
        "ATGTACGATC"     # substitution (C -> T)
        + "ACGTACGAT"    # deletion of trailing C
        + "ACGTACGATCA"  # insertion of A
    )
    assert len(read) == 30
    status, data = post(
        {"reference": ref, "read": read, "copies": 3, "max_edits": 1}
    )
    print("mixed-noise status:", status)
    print(json.dumps(data, indent=2)[:1600])
    if status != 200:
        print("FAIL: expected HTTP 200")
        return False
    if data["objective"] != {"total_edits": 3, "max_segment_edits": 1}:
        print("FAIL: unexpected objective", data["objective"])
        return False
    witness = data.get("witness") or data["witnesses"][0]
    ops = {
        op
        for seg in witness["segments"]
        for op in seg["cigar"]
        if op.isalpha()
    }
    if not {"M", "D", "I"} <= ops:
        print("FAIL: CIGAR must contain M, D and I; got", ops)
        return False
    # The replay rows must reconstruct both original strings.
    for seg in witness["segments"]:
        if seg["aligned_reference"].replace("-", "") != seg["reference"]:
            print("FAIL: reference replay mismatch")
            return False
        if seg["aligned_read"].replace("-", "") != seg["read"]:
            print("FAIL: read replay mismatch")
            return False
    print("mixed-noise smoke passed (substitution + deletion + insertion)")

    # Unique clean case.
    status, data = post(
        {"reference": ref, "read": ref * 3, "copies": 3, "max_edits": 0}
    )
    if status != 200 or data["status"] != "unique" or data["objective"][
        "total_edits"
    ] != 0:
        print("FAIL: clean unique case", status, data.get("status"))
        return False
    print("unique case passed")

    # Ambiguous case (homopolymer shifts).
    status, data = post(
        {
            "reference": "AAAAAAAAAA",
            "read": "A" * 30,
            "copies": 3,
            "max_edits": 0,
        }
    )
    if (
        status != 200
        or data["status"] != "ambiguous"
        or len(data["witnesses"]) != 2
    ):
        print("FAIL: ambiguous case", status, data.get("status"))
        return False
    print("ambiguous case passed (two witnesses returned)")

    # Infeasible case must be locatable (422 + violating segment).
    status, data = post(
        {
            "reference": ref,
            "read": "TTGTACGATC" + ref * 2,
            "copies": 3,
            "max_edits": 1,
        }
    )
    if status != 422 or data.get("status") != "infeasible":
        print("FAIL: infeasible case", status, data.get("status"))
        return False
    bad = data["nearest"]["violating_segments"]
    if not bad or bad[0]["required_edits"] != 2:
        print("FAIL: infeasible case not locatable", data.get("nearest"))
        return False
    print("infeasible case passed (violating segment located)")
    return True


def smoke_strands() -> bool:
    step("strand_mode smoke (forward / reverse / auto + cross-strand tie)")

    ref = "ACGTACGATC"
    rref = reverse_complement(ref)
    assert rref == "GATCGTACGT"

    # Strict forward cannot explain a clean reverse-complement read at cap 0.
    status, data = post(
        {
            "reference": ref,
            "read": rref * 3,
            "copies": 3,
            "max_edits": 0,
            "strand_mode": "forward",
        }
    )
    if status != 422 or data.get("status") != "infeasible":
        print("FAIL: reverse-complement read must be infeasible in forward", status)
        return False
    print("forward rejects reverse-complement read (422 infeasible)")

    # Reverse mode decodes it uniquely with strand == reverse.
    status, data = post(
        {
            "reference": ref,
            "read": rref * 3,
            "copies": 3,
            "max_edits": 0,
            "strand_mode": "reverse",
        }
    )
    if status != 200 or data["status"] != "unique":
        print("FAIL: reverse mode should decode", status, data.get("status"))
        return False
    witness = data["witness"]
    if witness["strand"] != "reverse" or witness["shift"] != 0:
        print("FAIL: unexpected reverse witness", witness.get("strand"), witness.get("shift"))
        return False
    if data["rotated_reference"] != rref:
        print("FAIL: oriented reference mismatch", data["rotated_reference"])
        return False
    # The CIGAR must replay directly against the oriented reference/read.
    seg = witness["segments"][0]
    if seg["aligned_reference"] != seg["aligned_read"] or seg["aligned_reference"].replace("-", "") != rref:
        print("FAIL: reverse replay mismatch")
        return False
    print("reverse decodes uniquely; CIGAR replays against oriented reference")

    # Auto must determine the direction as reverse (never prefer the
    # first-computed forward direction).
    status, data = post(
        {
            "reference": ref,
            "read": rref * 3,
            "copies": 3,
            "max_edits": 0,
            "strand_mode": "auto",
        }
    )
    if status != 200 or data["status"] != "unique":
        print("FAIL: auto should resolve uniquely", status, data.get("status"))
        return False
    if data["witness"]["strand"] != "reverse":
        print("FAIL: auto picked the wrong direction", data["witness"].get("strand"))
        return False
    print("auto determines direction = reverse")

    # Cross-strand tie: a non-periodic reverse-complement palindrome gives an
    # equal optimum in both directions, which must be reported ambiguous.
    pal = "AACGTACGTT"
    assert reverse_complement(pal) == pal and len(pal) == 10
    status, data = post(
        {
            "reference": pal,
            "read": pal * 3,
            "copies": 3,
            "max_edits": 0,
            "strand_mode": "auto",
        }
    )
    if status != 200 or data["status"] != "ambiguous":
        print("FAIL: cross-strand tie must be ambiguous", status, data.get("status"))
        return False
    strands = [w["strand"] for w in data["witnesses"]]
    if strands != ["forward", "reverse"]:
        print("FAIL: witnesses must be ordered forward then reverse", strands)
        return False
    if data.get("more_witnesses") is not False:
        print("FAIL: exactly two cross-strand witnesses expected")
        return False
    print("cross-strand tie reported ambiguous, stably ordered forward/reverse")

    # An illegal strand_mode is a field-validation rejection, not a decode.
    status, data = post(
        {
            "reference": ref,
            "read": rref * 3,
            "copies": 3,
            "max_edits": 0,
            "strand_mode": "backwards",
        }
    )
    if status != 422 or data.get("status") == "infeasible":
        print("FAIL: illegal strand_mode must be a field 422", status)
        return False
    locs = [tuple(err.get("loc", [])) for err in data.get("detail", [])]
    if ("body", "strand_mode") not in locs:
        print("FAIL: field error must point at strand_mode", locs)
        return False
    print("illegal strand_mode rejected as a field error")
    return True


def main() -> int:
    if not wait_for_health():
        return 1
    if not run_tests():
        return 1
    if not check_build():
        return 1
    if not smoke_decode():
        return 1
    if not smoke_strands():
        return 1
    print("\n=== verify: ALL CHECKS PASSED ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
