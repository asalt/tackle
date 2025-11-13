# formula_utils.py
"""
Utilities for extracting explicit GeneID references from tokens and strings.

Supported patterns (case-insensitive):

- Single token (e.g. from a parsed formula):
    geneid:ABC123
    GeneID_ABC123
    geneid-ABC123
    `geneid:ABC123`   # backtick-quoted
    # TODO also support case insensitive GID and GENE
    # TODO also test for if multiple are present e.g. gene_123-gene456 (with or without space hetween the operator) (we should probably use operators as part of our boundary txt extraction? )

- Inline in a larger string:
    "foo GeneID:ABC123 bar"
    "prefix geneid_ABC123.suffix"
"""

from __future__ import annotations

import re
from typing import Optional, List, Iterable


# Allowable characters in a GeneID key (tweak if needed)
_GENEID_KEY_CHARS = r"[A-Za-z0-9_.-]+"

# Case-insensitive matcher for explicit GeneID references as a *single token*
_GENEID_TOKEN_RE = re.compile(
    rf"^geneid_?({_GENEID_KEY_CHARS})$",
    re.IGNORECASE,
)

# Case-insensitive finder for explicit GeneID references *inside* larger strings
_GENEID_FINDER_RE = re.compile(
    rf"(?i)\bgeneid_?({_GENEID_KEY_CHARS})\b"
)


def _coerce_token(token: object | None) -> str:
    """Best-effort conversion to a normalized string token."""
    if token is None:
        return ""
    # Normalize whitespace around the token
    return str(token).strip()


def unquote_backticks(token: object | None) -> str:
    """
    Remove a single pair of surrounding backticks from a token, if present.

    Examples:
        `geneid:ABC`  -> geneid:ABC
        geneid:ABC    -> geneid:ABC
        ``weird``     -> `weird`   (only outermost pair is removed)
        None          -> ""
    """
    tok = _coerce_token(token)
    if len(tok) >= 2 and tok[0] == tok[-1] == "`":
        return tok[1:-1]
    return tok


def extract_geneid_key(token: object | None) -> Optional[str]:
    """
    Return the GeneID key if the token matches the explicit GeneID pattern, else None.

    Accepts:
        - Raw strings like "geneid:ABC123"
        - Backtick-quoted tokens like "`geneid_ABC123`"
        - Any object; will be coerced via str() (with None -> "")

    Returns:
        The extracted key (e.g. "ABC123") or None if no valid GeneID pattern is found.
    """
    tok = unquote_backticks(token)
    if not tok:
        return None

    m = _GENEID_TOKEN_RE.match(tok)
    if not m:
        return None

    key = m.group(1).strip()
    return key or None


def find_geneid_keys_in_string(s: object | None, *, unique: bool = False) -> List[str]:
    """
    Find all explicit GeneID keys in a larger string.

    Examples:
        "foo GeneID:ABC bar geneid_DEF baz" -> ["ABC", "DEF"]

    Args:
        s: Any object; None is treated as no matches.
        unique: If True, deduplicate while preserving the first occurrence order.

    Returns:
        List of matched keys (possibly empty).
    """
    if s is None:
        return []

    text = str(s)
    matches: List[str] = _GENEID_FINDER_RE.findall(text)

    if not unique:
        return matches

    # Deduplicate while preserving order
    seen: set[str] = set()
    uniq: List[str] = []
    for key in matches:
        if key not in seen:
            seen.add(key)
            uniq.append(key)
    return uniq
