import re
from typing import Optional, List

# Allowable characters in a GeneID key (tune as needed)
_GENEID_KEY_CHARS = r"[A-Za-z0-9_.-]+"

# Case-insensitive matcher for explicit GeneID references
# Only support underscore or direct adjacency after 'geneid' to avoid ':' and '-' clashes in formulas.
_GENEID_TOKEN_RE = re.compile(rf"^geneid_?({_GENEID_KEY_CHARS})$", re.IGNORECASE)
_GENEID_FINDER_RE = re.compile(rf"(?i)\bgeneid_?({_GENEID_KEY_CHARS})\b")


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


def find_geneid_keys_in_string(s: str) -> List[str]:
    """Find all explicit GeneID keys in a larger string."""
    return _GENEID_FINDER_RE.findall(str(s))
