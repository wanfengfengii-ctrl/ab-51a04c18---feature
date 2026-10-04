"""HTTP API for concatemer decoding.

POST /api/concatemers/decode
    {
      "reference": "8-20 nt circular reference",
      "read": "30-160 nt tandem read",
      "copies": 3-8,                 // expected tandem copy count
      "max_edits": 0-3,              // per-copy edit budget
      "strand_mode": "forward" | "reverse" | "auto"  // optional
    }

``strand_mode`` selects which reference orientation participates in the
common-cut-point decode: ``forward`` (the submitted reference), ``reverse``
(its reverse complement) or ``auto`` (both jointly).  When omitted the
request behaves exactly as the forward-only API.

Responses (HTTP 200):
    status == "unique"     -> single optimal explanation
    status == "ambiguous"  -> multiple optima, the first two witnesses
                              (sorted by strand, shift, boundaries, CIGAR)
                              are shown
Constraint failure:
    HTTP 422 with status == "infeasible" and a locatable ``nearest`` block.
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
            "reference orientation: 'forward' as submitted, 'reverse' uses "
            "the reverse complement, 'auto' adjudicates both jointly; "
            "omitted means forward with the legacy response shape"
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
    # Omitting strand_mode keeps the legacy forward-only contract exactly,
    # including the response shape (no strand fields).
    legacy = request.strand_mode is None
    strand_mode = "forward" if legacy else request.strand_mode
    result = solve(
        request.reference,
        request.read,
        request.copies,
        request.max_edits,
        strand_mode,
    )

    echo = {
        "reference": request.reference,
        "read": request.read,
        "copies": request.copies,
        "max_edits_per_segment": request.max_edits,
    }
    if not legacy:
        echo["strand_mode"] = strand_mode
    result["request"] = echo

    if result["status"] == "infeasible":
        if legacy and result.get("nearest") is not None:
            result["nearest"].pop("strand", None)
        return JSONResponse(status_code=422, content=result)

    def oriented_reference(strand: str) -> str:
        if strand == "reverse":
            return reverse_complement(request.reference)
        return request.reference

    # Attach the rotated, orientation-aware reference used by each witness.
    if result["status"] == "unique":
        witness = result["witness"]
        if legacy:
            witness.pop("strand", None)
            result["rotated_reference"] = rotate(
                request.reference, witness["shift"]
            )
            result["ordering"] = (
                "objective lexicographically minimizes (total_edits, "
                "max_segment_edits); ties ordered by (shift, boundaries, "
                "CIGAR)"
            )
        else:
            result["rotated_reference"] = rotate(
                oriented_reference(witness["strand"]), witness["shift"]
            )
            result["ordering"] = (
                "objective lexicographically minimizes (total_edits, "
                "max_segment_edits); ties ordered by (strand, shift, "
                "boundaries, CIGAR)"
            )
    else:
        for witness in result["witnesses"]:
            if legacy:
                witness.pop("strand", None)
                base = request.reference
            else:
                base = oriented_reference(witness["strand"])
            witness["rotated_reference"] = rotate(base, witness["shift"])
        result["ordering"] = (
            "witnesses sorted by ("
            + ("" if legacy else "strand, ")
            + "shift, boundaries, CIGAR); only the first two are returned"
        )
    return result
