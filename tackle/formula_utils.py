import re
from typing import Optional, List, Tuple, Iterable

# Allowable characters in a GeneID key (tune as needed)
_GENEID_KEY_CHARS = r"[A-Za-z0-9_.-]+"

# Case-insensitive matcher for explicit GeneID references
# Only support underscore or direct adjacency after 'geneid' to avoid ':' and '-' clashes in formulas.
_GENEID_TOKEN_RE = re.compile(rf"^geneid_?({_GENEID_KEY_CHARS})$", re.IGNORECASE)
_GID_TOKEN_RE = re.compile(rf"^gid_?({_GENEID_KEY_CHARS})$", re.IGNORECASE)
_GENEID_FINDER_RE = re.compile(rf"(?i)\bgeneid_?({_GENEID_KEY_CHARS})\b")
_GID_FINDER_RE = re.compile(rf"(?i)\bgid_?({_GENEID_KEY_CHARS})\b")


def unquote_backticks(token: str) -> str:
    """Remove a single pair of outer backticks from a token, if present."""
    if not isinstance(token, str):
        token = str(token)
    token = token.strip()
    if len(token) >= 2 and token[0] == token[-1] == "`":
        return token[1:-1]
    return token


def extract_geneid_key(token: str) -> Optional[str]:
    """Return the GeneID key if token matches the explicit GeneID pattern, else None."""
    tok = unquote_backticks(token)
    m = _GENEID_TOKEN_RE.match(tok)
    if not m:
        return None
    key = m.group(1).strip()
    return key or None


def extract_gid_key(token: str) -> Optional[str]:
    """Return the key if token matches an explicit GID pattern, else None."""
    tok = unquote_backticks(token)
    m = _GID_TOKEN_RE.match(tok)
    if not m:
        return None
    key = m.group(1).strip()
    return key or None


def extract_gene_key_any(token: object | None) -> Optional[Tuple[str, str]]:
    """Return (prefix, key) for explicit gene tokens: ('geneid'|'gid', key) or None."""
    tok = unquote_backticks(token)
    if not tok:
        return None
    m = _GENEID_TOKEN_RE.match(tok)
    if m:
        key = m.group(1).strip()
        return ("geneid", key) if key else None
    m = _GID_TOKEN_RE.match(tok)
    if m:
        key = m.group(1).strip()
        return ("gid", key) if key else None
    return None


def find_geneid_keys_in_string(s: str) -> List[str]:
    """Find all explicit GeneID keys in a larger string."""
    return _GENEID_FINDER_RE.findall(str(s))


def find_gid_keys_in_string(s: str) -> List[str]:
    """Find all explicit GID keys in a larger string."""
    return _GID_FINDER_RE.findall(str(s))


def find_gene_keys_in_string_any(s: str, *, unique: bool = False) -> List[str]:
    """Find keys from either GeneID_* or GID_* tokens in text."""
    text = str(s)
    keys = _GENEID_FINDER_RE.findall(text) + _GID_FINDER_RE.findall(text)
    if not unique:
        return keys
    seen = set()
    out: List[str] = []
    for k in keys:
        if k not in seen:
            seen.add(k)
            out.append(k)
    return out
