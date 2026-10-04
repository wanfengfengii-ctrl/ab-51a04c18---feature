"""HTTP API for concatemer decoding.

POST /api/concatemers/decode
    {
      "reference": "8-20 nt circular reference",
      "read": "30-160 nt tandem read",
      "copies": 3-8,                 // expected tandem copy count
      "max_edits": 0-3,              // per-copy edit budget
      "strand_mode": "forward"       // optional: forward | reverse | auto
    }

``strand_mode`` (omit or "forward" by behaviour):
    forward - decode against the submitted reference (historical behaviour),
    reverse - decode against the reverse complement of the reference,
    auto    - consider both orientations jointly; a tie across directions is
              reported as ambiguous, never resolved in favour of the direction
              computed first.

Responses (HTTP 200):
    status == "unique"     -> single optimal explanation
    status == "ambiguous"  -> multiple optima, the first two witnesses
                              (sorted by strand, shift, boundaries, CIGAR)
                              are shown
Constraint failure:
    HTTP 422 with status == "infeasible" and a locatable ``nearest`` block.
An illegal ``strand_mode`` is rejected as a field-validation error (422).
"""

from __future__ import annotations

import os
from typing import Literal, Optional

from fastapi import FastAPI
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field, field_validator

from .solver import reverse_complement, rotate, solve

app = FastAPI(
    title="Concatemer Decode API",
    version="1.1.0",
    description="Recover a common cut point from noisy tandem barcode reads.",
)

_DNA = set("ACGT")


class DecodeRequest(BaseModel):
    reference: str = Field(..., description="8-20 nt circular reference")
    read: str = Field(..., alias="read", description="30-160 nt tandem read")
    copies: int = Field(..., ge=3, le=8)
    max_edits: int = Field(..., ge=0, le=3)
    strand_mode: Optional[Literal["forward", "reverse", "auto"]] = Field(
        default=None,
        description=(
            "Reference orientation: forward (default), reverse complement, "
            "or auto to consider both and report a cross-strand tie as "
            "ambiguous."
        ),
    )

    model_config = {"populate_by_name": True, "extra": "ignore"}

    @field_validator("reference", "read")
    @classmethod
    def _validate_dna(cls, value: str) -> str:
        seq = value.strip().upper()
        if not seq:
            raise ValueError("sequence must be non-empty")
        bad = sorted({ch for ch in seq if ch not in _DNA})
        if bad:
            raise ValueError(
                f"invalid DNA character(s): {''.join(bad)!r}; only A/C/G/T allowed"
            )
        return seq

    @field_validator("reference")
    @classmethod
    def _check_reference_length(cls, value: str) -> str:
        if not 8 <= len(value) <= 20:
            raise ValueError("reference must be 8 to 20 nt long")
        return value

    @field_validator("read")
    @classmethod
    def _check_read_length(cls, value: str) -> str:
        if not 30 <= len(value) <= 160:
            raise ValueError("read must be 30 to 160 nt long")
        return value


class Health(BaseModel):
    status: str
    service: str


@app.get("/health", response_model=Health)
def health() -> Health:
    return Health(status="ok", service="concatemer-decode")


@app.post("/api/concatemers/decode")
def decode(request: DecodeRequest):
    # An omitted strand_mode preserves the original forward request, response
    # and failure behaviour exactly; an explicit "forward" is strand-aware.
    omitted = request.strand_mode is None
    mode = "forward" if omitted else request.strand_mode
    result = solve(
        request.reference,
        request.read,
        request.copies,
        request.max_edits,
        mode,
    )

    echo = {
        "reference": request.reference,
        "read": request.read,
        "copies": request.copies,
        "max_edits_per_segment": request.max_edits,
    }
    if not omitted:
        echo["strand_mode"] = request.strand_mode
    result["request"] = echo

    if result["status"] == "infeasible":
        if omitted:
            _legacy_infeasible(result)
        return JSONResponse(status_code=422, content=result)

    if omitted:
        # Drop the strand label so the historical success body is unchanged.
        if result["status"] == "unique":
            result["witness"].pop("strand", None)
        else:
            for witness in result["witnesses"]:
                witness.pop("strand", None)

    def oriented_reference(witness: dict) -> str:
        strand = witness.get("strand", "forward")
        base = (
            reverse_complement(request.reference)
            if strand == "reverse"
            else request.reference
        )
        return rotate(base, witness["shift"])

    # The strand field only appears for an explicitly supplied strand_mode, so
    # omitted requests keep ordering text without a strand component.
    strand_aware = not omitted

    # Attach the oriented reference used by the optimal explanation(s).
    if result["status"] == "unique":
        witness = result["witness"]
        result["rotated_reference"] = oriented_reference(witness)
        result["ordering"] = (
            "objective lexicographically minimizes (total_edits, "
            "max_segment_edits); ties ordered by (strand, shift, "
            "boundaries, CIGAR)"
            if strand_aware
            else "objective lexicographically minimizes (total_edits, "
            "max_segment_edits); ties ordered by (shift, boundaries, CIGAR)"
        )
    else:
        for witness in result["witnesses"]:
            witness["rotated_reference"] = oriented_reference(witness)
        result["ordering"] = (
            "witnesses sorted by (strand, shift, boundaries, CIGAR); "
            "only the first two are returned"
            if strand_aware
            else "witnesses sorted by (shift, boundaries, CIGAR); "
            "only the first two are returned"
        )
    return result


def _legacy_infeasible(result: dict) -> None:
    """Collapse a strand-aware infeasible envelope to the original shape."""
    constraint = result["constraint"]
    constraint.pop("strand_mode", None)
    feasible_strands = result.pop("feasible_strands", {})
    request_echo = result.get("request")
    nearest = result.get("nearest")
    if nearest is not None:
        # Rebuild to keep the historical key order (shift first).
        nearest = {
            "shift": nearest["shift"],
            "total_edits": nearest["total_edits"],
            "max_segment_edits": nearest["max_segment_edits"],
            "boundaries": nearest["boundaries"],
            "violating_segments": nearest["violating_segments"],
        }
    # Rebuild the top level to preserve the historical key order; ``request``
    # was appended after the envelope in the original code path as well.
    rebuilt = {
        "status": result["status"],
        "error": result["error"],
        "message": result["message"],
        "constraint": constraint,
        "feasible_shifts": feasible_strands.get("forward", []),
        "nearest": nearest,
    }
    if request_echo is not None:
        rebuilt["request"] = request_echo
    result.clear()
    result.update(rebuilt)
