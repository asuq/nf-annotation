#!/usr/bin/env python3
"""Parse and classify metadata atypical-warning strings."""

from __future__ import annotations

import re

MISSING_VALUE_TOKENS = {"", "na", "n/a", "null", "none"}
UNVERIFIED_SOURCE_ORGANISM = "unverified source organism"
ATYPICAL_REASON_SPLIT_RE = re.compile(r"[;,]")


def is_missing(value: str | None) -> bool:
    """Return True when one atypical-warning value should be treated as missing."""
    if value is None:
        return True
    return value.strip().lower() in MISSING_VALUE_TOKENS


def parse_atypical_reasons(value: str | None) -> tuple[str, ...]:
    """Return normalised atypical reasons extracted from one metadata value."""
    if is_missing(value):
        return ()
    assert value is not None
    return tuple(
        reason.casefold()
        for raw_reason in ATYPICAL_REASON_SPLIT_RE.split(value)
        if (reason := raw_reason.strip())
    )


def classify_atypical_warnings(value: str | None) -> tuple[bool, bool]:
    """Return atypical status and the sole unverified-source exception flag."""
    reasons = parse_atypical_reasons(value)
    if not reasons:
        return False, False
    return True, all(reason == UNVERIFIED_SOURCE_ORGANISM for reason in reasons)


def is_unverified_source_only(value: str | None) -> bool:
    """Return True when the only atypical reason is unverified source organism."""
    _is_atypical, is_exception = classify_atypical_warnings(value)
    return is_exception
