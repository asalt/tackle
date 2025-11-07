import re
from typing import Optional, List

# Case-insensitive matcher for explicit GeneID references
_GENEID_TOKEN_RE = re.compile(r"^geneid(?:[_:\-])(.+)$", re.IGNORECASE)
_GENEID_FINDER_RE = re.compile(r"(?i)\bgeneid(?:[_:\-])([A-Za-z0-9_.-]+)\b")


def unquote_backticks(token: str) -> str:
    if len(token) >= 2 and token[0] == "`" and token[-1] == "`":
        return token[1:-1]
    return token


def extract_geneid_key(token: str) -> Optional[str]:
    """Return the GeneID key if token matches the explicit GeneID pattern, else None."""
    tok = unquote_backticks(str(token))
    m = _GENEID_TOKEN_RE.match(tok)
    if not m:
        return None
    key = m.group(1)
    return key if key else None


def find_geneid_keys_in_string(s: str) -> List[str]:
    """Find all explicit GeneID keys in a larger string."""
    return _GENEID_FINDER_RE.findall(str(s))

